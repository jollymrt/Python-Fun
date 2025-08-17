{\rtf1\ansi\ansicpg1252\cocoartf2820
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 import pysam\
from Bio import pairwise2\
from collections import defaultdict, Counter\
import sys\
\
def primer_match(read_seq, primer, min_identity=0.9):\
    """Check if read has >=90% identity to primer (local alignment)."""\
    aln = pairwise2.align.localms(read_seq, primer, 2, -1, -2, -0.5, one_alignment_only=True)[0]\
    matches = sum([1 for a, b in zip(aln.seqA, aln.seqB) if a == b])\
    identity = matches / len(primer)\
    return identity >= min_identity\
\
def filter_read(read, primer):\
    """Filter reads based on primer, clipping, and alignment flags."""\
    if read.is_secondary or read.is_supplementary:\
        return False\
    seq = read.query_sequence\
    # primer check\
    if not primer_match(seq, primer):\
        return False\
    # soft/hard clipping filter\
    clip_bases = 0\
    read_len = read.query_length\
    for (cigar_type, length) in read.cigartuples:\
        if cigar_type in (4, 5):  # 4 = soft clip, 5 = hard clip\
            clip_bases += length\
    if clip_bases / read_len > 0.05:  # >5% clipped\
        return False\
    return True\
\
def get_variants(read, region_start):\
    """Extract SNVs and indels adjusted for position in reference coordinates."""\
    variants = []\
    ref_pos = read.reference_start\
    query_pos = 0\
\
    for (cigar_type, length) in read.cigartuples:\
        if cigar_type == 0:  # match/mismatch\
            for i in range(length):\
                ref_base = read.get_reference_sequence()[i + query_pos]\
                read_base = read.query_sequence[i + query_pos]\
                if ref_base != read_base:\
                    variants.append(("SNV", ref_pos + i, ref_base, read_base))\
            ref_pos += length\
            query_pos += length\
\
        elif cigar_type == 1:  # insertion\
            ins_seq = read.query_sequence[query_pos: query_pos + length]\
            variants.append(("INS", ref_pos, "-", ins_seq))\
            query_pos += length\
\
        elif cigar_type == 2:  # deletion\
            variants.append(("DEL", ref_pos, read.get_reference_sequence()[query_pos: query_pos + length], "-"))\
            ref_pos += length\
\
        else:\
            # skip softclip, hardclip, padding etc.\
            if cigar_type in (4, 5):\
                query_pos += length\
\
    return variants\
\
def main():\
    if len(sys.argv) < 4:\
        print("Usage: python umi_variant_caller.py input.bam chr:start-end primer_seq")\
        sys.exit(1)\
\
    bam_file = sys.argv[1]\
    region = sys.argv[2]\
    primer = sys.argv[3]\
\
    samfile = pysam.AlignmentFile(bam_file, "rb")\
    chrom, coords = region.split(":")\
    start, end = map(int, coords.split("-"))\
\
    umi_dict = defaultdict(list)\
\
    for read in samfile.fetch(chrom, start, end):\
        if not filter_read(read, primer):\
            continue\
        umi = read.get_tag("RX") if read.has_tag("RX") else read.query_name\
        umi_dict[umi].append(read)\
\
    print("UMI\\tVariantType\\tPosition\\tRef\\tAlt\\tReadCount")\
\
    for umi, reads in umi_dict.items():\
        all_variants = []\
        for read in reads:\
            variants = get_variants(read, start)\
            all_variants.extend(variants)\
\
        variant_counts = Counter(all_variants)\
        for var, count in variant_counts.items():\
            vtype, pos, ref, alt = var\
            print(f"\{umi\}\\t\{vtype\}\\t\{pos\}\\t\{ref\}\\t\{alt\}\\t\{count\}")\
\
    samfile.close()\
\
if __name__ == "__main__":\
    main()}