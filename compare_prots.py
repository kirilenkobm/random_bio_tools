#!/usr/bin/env python3
"""Compare two protein sequences."""
import argparse
import os
import sys
from collections import defaultdict

__author__ = "Bogdan Kirilenko, 2018."
EPSTEINS_MATRIX_PATH = os.path.join(os.path.dirname(__file__), "data", "Epsteins_difference.txt")
BLOSUM62_MATRIX_PATH = os.path.join(os.path.dirname(__file__), "data", "BLOSUM62.txt")
GAP_SCORE = 0  # score for -
MASK_SCORE = 0  # score for X
STOP_SCORE = 0  # score for *


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def read_fasta(fasta_stream):
    """Read fasta, return dict and type."""
    # open the file
    f = open(fasta_stream, "r") if fasta_stream != "stdin" else sys.stdin
    fasta_data = f.read().split(">")
    f.close()
    # assert fasta_data[0] == ""  # if a file starts with > it should be empty
    del fasta_data[0]  # remove it "" we don't need that
    sequences = {}  # accumulate data here
    order = []  # to have ordered list
    # read line by line

    for elem in fasta_data:
        raw_lines = elem.split("\n")
        header = raw_lines[0]  # it must be first ['capHir1', 'ATGCCGCGCCAATTCCCCAAGCTGA... ]

        lines = [x for x in raw_lines[1:] if x != ""]  # separate nucleotide-containing lines
        if len(lines) == 0:  # it is a mistake - empty sequene --> get rid of
            continue

        fasta_content = "".join(lines)
        sequences[header] = fasta_content
        order.append(header)
    if len(sequences) == 0:
        die("There are not fasta-formatted sequences in {0}!".format(fasta_stream))
    if len(sequences.keys()) != len(order):  # it is possible in case of non-unique headers
        err = "Error! Sequences names must be unique! There are" \
              " {0} sequences and {1} unique names!".format(len(sequences.keys()), len(order))
        die(err)
    return sequences


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("fasta_1", help="Fasta file containing the first sequence")
    app.add_argument("seq_1", help="Sequence 1 identifyer")
    app.add_argument("fasta_2", help="Fasta file containing the second sequence. "
                     "Write - if the first file is same")
    app.add_argument("seq_2", help="Second sequence identifyer, - if the same with seq_1")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def make_epsteins_matrix(matrix_path=EPSTEINS_MATRIX_PATH):
    """Read epsteins difference matrix."""
    # TODO: different matrixes
    f = open(matrix_path, "r")
    MATRIX, num = defaultdict(dict), 0
    for line_num, line in enumerate(f):
        if line_num == 0:
            alphabet = line.split()
            continue
        # alphabet is given
        assert len(alphabet) > 0  # should be obtained before we reach this path
        line_data = line.split()
        del line_data[0]
        scores = [float(x) for x in line_data]
        row_char = alphabet[num]
        for subnum, score in enumerate(scores):
            col_char = alphabet[subnum]
            MATRIX[row_char][col_char] = score
        num += 1
    f.close()
    return MATRIX


def make_blosum_matrix(matrix_path=BLOSUM62_MATRIX_PATH):
    """Read BLOSUM62 matrix."""
    # TODO: different matrixes
    f = open(matrix_path, "r")
    MATRIX, num = defaultdict(dict), 0
    for line in (f):
        if line.startswith("#"):
            continue
        elif line.startswith(" "):
            alphabet = line.split()
            continue
        # alphabet is given
        assert len(alphabet) > 0  # should be obtained before we reach this path
        line_data = line.split()
        del line_data[0]
        scores = [int(x) for x in line_data]
        row_char = alphabet[num]
        for subnum, score in enumerate(scores):
            col_char = alphabet[subnum]
            MATRIX[row_char][col_char] = score
        num += 1
    f.close()
    return MATRIX


def main():
    """Entry point."""
    args = parse_args()
    # read data and check if it's correct
    fasta_1_data = read_fasta(args.fasta_1)
    fasta_2_data = read_fasta(args.fasta_2) if args.fasta_2 != "-" else fasta_1_data
    seq_1 = fasta_1_data.get(args.seq_1)
    seq_2 = fasta_2_data.get(args.seq_2)
    die("Error! Sequence {} not found in {}!".format(args.seq_1, args.fasta_1)) if not seq_1 else None
    die("Error! Sequence {} not found in {}!".format(args.seq_2, args.fasta_2)) if not seq_2 else None
    die("Error! Aligned sequences required! (seq_1 and seq_2 have different lenghts)")\
        if len(seq_1) != len(seq_2) else None
    # so let's get started
    epsteins_matrix = make_epsteins_matrix()
    BLOSUM62_matxix = make_blosum_matrix()
    identical_aa = 0
    seq_len = comp_seq_len = len(seq_1)
    eps_sum_score, blosum_sum_score = 0, 0
    for ch1, ch2 in zip(seq_1, seq_2):
        # ignored cases
        if ch1 == ch2 == "X":
            # both masked
            comp_seq_len -= 1
            continue
        elif ch1 == ch2 == "-":
            # both gaps
            comp_seq_len -= 1
            continue
        elif ch1 == ch2 == "*":
            # both stops
            comp_seq_len -= 1
            continue
        # scores to be computed
        elif ch1 == ch2:
            ep_score = 1
            bl_score = BLOSUM62_matxix[ch1][ch2]
            identical_aa += 1
        elif ch1 == "-" or ch2 == "-":
            ep_score = GAP_SCORE
            bl_score = GAP_SCORE
        elif ch1 == "X" or ch2 == "X":
            ep_score = MASK_SCORE
            bl_score = BLOSUM62_matxix[ch1][ch2]
        elif ch1 == "*" or ch2 == "*":
            ep_score == STOP_SCORE
            bl_score = BLOSUM62_matxix[ch1][ch2]
        else:
            bl_score = BLOSUM62_matxix[ch1][ch2]
            ep_score = epsteins_matrix[ch1][ch2]
        eps_sum_score += ep_score
        blosum_sum_score += bl_score
    # output
    ep_perc_sim = eps_sum_score / comp_seq_len * 100
    perc_id = identical_aa / comp_seq_len * 100
    print("Percent similarity according the Epsteins matrix: {}".format(ep_perc_sim))
    print("Percent identity: {}".format(perc_id))
    print("Alignment score according BLOSUM62: {}".format(blosum_sum_score))
    print("Identical amino acids: {}".format(identical_aa))
    print("Sequence length: {}".format(seq_len))
    print("Effective sequence length: {}".format(comp_seq_len))
    sys.exit(0)


if __name__ == "__main__":
    main()
