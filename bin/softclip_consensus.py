#!/usr/bin/env python

# Written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.

"""Extend the ends of a consensus with the consensus of the reads soft-clipped at those ends.

Pileup-based consensus callers (iVar, bcftools, ViralConsensus) only see aligned bases, so a
consensus can never grow past the reference it was mapped to: reads that run over an end are
soft-clipped and their overhang is dropped. Iterative refinement therefore cannot recover a
segment end that the assembler left off.

This script adds that overhang back. For each end of each contig it takes the clipped reads
anchored at the tip, groups them by where their alignment stops (reads shifted by an indel in the
contig's last few bases disagree position by position, so only the dominant group votes), walks
outward one position at a time and appends the majority base while
  - only reads still consistent with the extension so far vote (<= 1 mismatch per 20 bases),
    so a minority carrying an indel drops out instead of blurring every later position,
  - the majority base is carried by >= --min-depth unique fragments (not reads: a PCR duplicate
    family agreeing with itself is not independent evidence), and
  - >= --min-agreement of the reads at that position agree on it.
It then refuses any extension that creates a --kmer-size k-mer already present in the contig, in
either orientation. Assemblers fold segment ends back onto themselves (the extension is the
reverse complement of the true end); this keeps that artefact from being grown back in.

Meant to run between refinement iterations: the next iteration maps the reads onto the extended
contig and its consensus caller re-calls the new bases like any others.
"""

import argparse
import csv
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = set("ACGT")
BAM_CSOFT_CLIP, BAM_CHARD_CLIP = 4, 5


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extend consensus ends with the consensus of soft-clipped read overhangs.",
        epilog="Example: softclip_consensus.py --bam it1.bam --fasta it1.consensus.fa --prefix it1",
    )
    parser.add_argument("--bam", type=Path, required=True, help="Reads mapped to the reference the consensus was called from")
    parser.add_argument("--fasta", type=Path, required=True, help="Consensus fasta called from --bam")
    parser.add_argument("--prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("--min-depth", type=int, default=5, help="Minimum unique fragments carrying the majority base (default: 5)")
    parser.add_argument("--min-agreement", type=float, default=0.7, help="Minimum fraction of reads agreeing on the majority base (default: 0.7)")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality of a read to be used (default: 20)")
    parser.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality of an overhang base to be counted (default: 20)")
    parser.add_argument("--max-extension", type=int, default=150, help="Maximum bases added per end per run (default: 150)")
    parser.add_argument("--end-tolerance", type=int, default=5, help="Read alignment must start/stop within this many bases of the tip (default: 5)")
    parser.add_argument("--kmer-size", type=int, default=25, help="k-mer size of the self-repeat guard, 0 disables it (default: 25)")
    parser.add_argument(
        "-l", "--log-level", choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"), default="INFO",
        help="The desired log level (default INFO).",
    )
    return parser.parse_args(argv)


def fragment_key(read):
    """Identify the sequenced molecule, so duplicates of one fragment count once."""
    return (read.reference_start, read.reference_end, read.next_reference_start, read.template_length, read.is_reverse)


def usable(read, min_mapq):
    return not (
        read.is_unmapped or read.is_secondary or read.is_supplementary
        or read.is_qcfail or read.is_duplicate or read.mapping_quality < min_mapq
    )


def clip_length(cigar, first):
    """Soft-clip length at the start (first=True) or end of a CIGAR, looking past a hard clip."""
    ops = cigar if first else cigar[::-1]
    for op, length in ops[:2]:
        if op == BAM_CSOFT_CLIP:
            return length
        if op != BAM_CHARD_CLIP:
            return 0
    return 0


def overhang(reads, ref_len, side, args):
    """Build the sequence the clipped reads place beyond one end.

    Reads are grouped by how many contig bases lie beyond their alignment (the 'anchor offset').
    Reads in one group share a coordinate frame; reads in different groups can be shifted against
    each other when the contig's last few bases carry an indel, so only the dominant group votes.
    Its clip also re-calls the contig bases it overlaps, which the aligner refused to align.

    Returns (offset, outward sequence, reads at tip, reads clipped, dominant group share, reason).
    The outward sequence starts at the first replaced base (offset bases inside the contig) and
    runs away from the contig: for the 5' end, base 0 is at coordinate offset-1.
    """
    groups = defaultdict(list)  # offset -> [(fragment key, outward bases, None where low quality)]
    reaching = clipped = 0
    for read in reads:
        if not usable(read, args.min_mapq) or read.query_sequence is None:
            continue
        seq, qual = read.query_sequence, read.query_qualities
        if side == 5:
            offset = read.reference_start
            if offset > args.end_tolerance:
                continue
            reaching += 1
            clip = clip_length(read.cigartuples, first=True)
            # clip base j sits at coordinate offset - clip + j; walk outward from offset - 1
            idx = range(clip - 1, -1, -1)
        else:
            offset = ref_len - read.reference_end
            if offset > args.end_tolerance:
                continue
            reaching += 1
            clip = clip_length(read.cigartuples, first=False)
            # clip base j sits at coordinate reference_end + j
            idx = range(len(seq) - clip, len(seq))
        if clip <= offset:
            continue  # the clip does not reach past the tip
        clipped += 1
        bases = [seq[q].upper() if seq[q].upper() in BASES and (qual is None or qual[q] >= args.min_baseq) else None
                 for q in idx]
        groups[offset].append((fragment_key(read), bases))

    if not groups:
        return 0, "", reaching, clipped, 0.0, "no_reads"
    n_frags = {g: len({key for key, _ in reads_}) for g, reads_ in groups.items()}
    offset = max(n_frags, key=n_frags.get)
    share = n_frags[offset] / len({key for reads_ in groups.values() for key, _ in reads_})

    # Thread the reads: a read votes only while it agrees with the extension built so far, so reads
    # that carry an indel or a different haplotype drop out instead of blurring every later column.
    members = groups[offset]
    mismatches = [0] * len(members)
    out, reason = [], "max_extension"
    for i in range(offset + args.max_extension):
        tally, frags = Counter(), defaultdict(set)
        for m, (key, bases) in enumerate(members):
            if i < len(bases) and bases[i] and mismatches[m] <= max(1, i // 20):
                tally[bases[i]] += 1
                frags[bases[i]].add(key)
        if not tally:
            reason = "no_reads"
            break
        base, n = tally.most_common(1)[0]
        if n / sum(tally.values()) < args.min_agreement:
            reason = "disagreement"
            break
        if len(frags[base]) < args.min_depth:
            reason = "low_depth"
            break
        out.append(base)
        for m, (_, bases) in enumerate(members):
            if i < len(bases) and bases[i] and bases[i] != base:
                mismatches[m] += 1
    return offset, "".join(out), reaching, clipped, round(share, 3), reason


def kmers(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i + k]}


def guard(tail, ext, known, k):
    """Length of `ext` (appended after `tail`) kept before it creates a k-mer already in `known`.

    `tail` is the last k-1 contig bases, so every k-mer checked holds at least one new base.
    """
    if not ext or k <= 0:
        return len(ext)
    combined = tail + ext
    for i in range(len(combined) - k + 1):
        if combined[i:i + k] in known:
            return max(0, i - len(tail))
    return len(ext)


def extend(name, seq, bam, contig, args):
    """Return (new sequence, report rows) for one contig."""
    ref_len = bam.get_reference_length(contig)
    if ref_len != len(seq):
        # Indels in the consensus shift only interior coordinates; the ends still line up.
        logger.debug("%s: consensus %d bp, mapped reference %d bp", name, len(seq), ref_len)
    k = args.kmer_size
    known = kmers(seq, k) | kmers(revcomp(seq), k) if k > 0 else set()
    rows, start, stop, prefix, suffix = [], 0, len(seq), "", ""
    for side in (5, 3):
        tip = seq[0] if side == 5 else seq[-1]
        row = dict(contig=name, end=f"{side}'", reads_at_tip=0, reads_clipped=0, anchor_share=0.0,
                   replaced=0, extension=0, guard_trimmed=0)
        if tip == "N":
            rows.append({**row, "stop": "tip_masked"})
            continue
        lo, hi = (0, args.end_tolerance + 1) if side == 5 else (max(0, ref_len - args.end_tolerance - 1), ref_len)
        offset, out, reaching, clipped, share, reason = overhang(bam.fetch(contig, lo, hi), ref_len, side, args)
        row.update(reads_at_tip=reaching, reads_clipped=clipped, anchor_share=share)
        if side == 3:
            kept_seq = seq[:len(seq) - offset]
            keep = guard(kept_seq[-(k - 1):] if k > 1 else "", out, known, k)
            new = out[:keep]
        else:
            kept_seq = seq[offset:]
            # the 5' end is the 3' end of the reverse complement; the outward bases are its extension
            comp = out.translate(COMPLEMENT)
            keep = guard(revcomp(kept_seq[:k - 1]) if k > 1 else "", comp, known, k)
            new = revcomp(comp[:keep])
        if keep < len(out):
            reason = "self_repeat"
        if keep <= offset:
            # nothing beyond the tip survives: leave this end exactly as it was
            rows.append({**row, "guard_trimmed": len(out) - keep, "stop": reason})
            continue
        if side == 3:
            stop, suffix = len(seq) - offset, new
        else:
            start, prefix = offset, new
        rows.append({**row, "replaced": offset, "extension": keep - offset,
                     "guard_trimmed": len(out) - keep, "stop": reason})
    return prefix + seq[start:stop] + suffix, rows


def open_indexed(path):
    bam = pysam.AlignmentFile(str(path), "rb")
    if bam.has_index():
        return bam
    bam.close()
    try:
        logger.info("Indexing %s", path)
        pysam.index(str(path))
    except pysam.utils.SamtoolsError:
        logger.info("%s is not coordinate sorted, sorting a copy", path)
        sorted_path = Path(f"{path.stem}.sorted.bam")
        pysam.sort("-o", str(sorted_path), str(path))
        pysam.index(str(sorted_path))
        path = sorted_path
    return pysam.AlignmentFile(str(path), "rb")


def main(argv=None):
    args = parse_args(argv)
    logger.setLevel(args.log_level)
    records = [(r.name, r.comment, r.sequence.upper()) for r in pysam.FastxFile(str(args.fasta))]
    bam = open_indexed(args.bam)
    references = list(bam.references)

    # The consensus is renamed after calling, so pair records with BAM references by name, else by order.
    if all(name in references for name, _, _ in records):
        pairs = [(rec, rec[0]) for rec in records]
    elif len(records) == len(references):
        pairs = list(zip(records, references))
    else:
        logger.warning("Cannot pair %d consensus records with %d BAM references, leaving the consensus unchanged",
                       len(records), len(references))
        pairs = [(rec, None) for rec in records]

    report = []
    with open(f"{args.prefix}.fasta", "w") as out:
        for (name, comment, seq), contig in pairs:
            if contig is not None and seq:
                new, rows = extend(name, seq, bam, contig, args)
                report.extend(rows)
                logger.info("%s: %d -> %d bp", name, len(seq), len(new))
                seq = new
            out.write(f">{name}{' ' + comment if comment else ''}\n{seq}\n")

    fields = ["contig", "end", "reads_at_tip", "reads_clipped", "anchor_share", "replaced", "extension", "guard_trimmed", "stop"]
    with open(f"{args.prefix}.softclip.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
