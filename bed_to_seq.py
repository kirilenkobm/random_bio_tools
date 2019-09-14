#!/usr/bin/env python3
"""Convert bed lines to sequences."""
import argparse
import os
import sys
from twobitreader import TwoBitFile

__author__ = "Bogdan Kirilenko, 2018."
complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("bed_source", help="bed 12 file or stdin")
    app.add_argument("db", help="2 bit file or alias")
    app.add_argument("--utr", "-u", help="Load UTR sequences too", action="store_true", dest="utr")
    app.add_argument("--output", default="stdout", help="Output, stdout as default")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def get_2bit_path(db_opt):
    """Check if alias and return a path to 2bit file."""
    if os.path.isfile(db_opt):  # not an alias
        return db_opt  # there is nothing to do
    else:  # not a file
        raise FileNotFoundError(f"{db_opt} not a regular file!")


def revert_compl(string):
    """Revert string; changes bases with complementary."""
    string = string[::-1]
    new_str = ''
    for c in string:
        new_c = complement.get(c)
        new_str += new_c
    return new_str


def main():
    """Entry point."""
    args = parse_args()
    source = open(args.bed_source, "r") if args.bed_source != "stdin" else sys.stdin
    two_bit_data = TwoBitFile(get_2bit_path(args.db))
    # so let's read input
    for num, line in enumerate(source):
        bed_info = line[:-1].split("\t")
        # parse bed info
        chrom = bed_info[0]
        chrom_seq = two_bit_data[chrom]
        gene_seq = ""
        chromStart = int(bed_info[1])
        # chromEnd = int(bed_info[2])
        name = bed_info[3]  # gene_name usually
        # bed_score = int(bed_info[4])  # never used
        # strand = bed_info[5]  # otherwise:
        strand = True if bed_info[5] == '+' else False
        thickStart = int(bed_info[6])
        thickEnd = int(bed_info[7])
        # itemRgb = bed_info[8]  # never used
        blockCount = int(bed_info[9])
        blockSizes = [int(x) for x in bed_info[10].split(',') if x != '']
        blockStarts = [int(x) for x in bed_info[11].split(',') if x != '']
        # not-in-file info
        blockEnds = [blockStarts[i] + blockSizes[i] for i in range(blockCount)]
        blockAbsStarts = [blockStarts[i] + chromStart for i in range(blockCount)]
        blockAbsEnds = [blockEnds[i] + chromStart for i in range(blockCount)]
        # block-by-block
        for block_num in range(blockCount):
            if not args.utr:
                blockStart = blockAbsStarts[block_num]
                blockEnd = blockAbsEnds[block_num]
                # skip the block if it is entirely UTR
                if blockEnd <= thickStart:
                    continue
                elif blockStart >= thickEnd:
                    continue
                blockNewStart = blockStart if blockStart >= thickStart else thickStart
                blockNewEnd = blockEnd if blockEnd <= thickEnd else thickEnd
                exon_seq = chrom_seq[blockNewStart: blockNewEnd].upper()
            else:
                exon_seq = chrom_seq[blockAbsStarts[block_num]: blockAbsEnds[block_num]]
            gene_seq += exon_seq
        if len(gene_seq) == 0:
            continue
        gene_seq = gene_seq if strand else revert_compl(gene_seq)
        sys.stdout.write(">{}\n{}\n".format(name, gene_seq))
    source.close() if args.bed_source != "stdin" else None
    sys.exit(0)


if __name__ == "__main__":
    main()
