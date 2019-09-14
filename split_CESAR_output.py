#!/usr/bin/env python3
"""Split CESAR output in separated exons data."""
import argparse
import os
import sys
import subprocess

__author__ = "Bogdan Kirilenko, 2019."


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
    app.add_argument("cesar_output")
    app.add_argument("--flank_size", type=int, default=10)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def read_cesar_out(cesar_line):
    """Return ref and query sequence."""
    cesar_content = cesar_line.split("\n")
    # del cesar_content[0]
    fractions = parts(cesar_content, 4)
    cesar_fractions = []
    for fraction in fractions:
        if len(fraction) == 1:
            continue
        ref_seq = fraction[1]
        query_name = fraction[2][1:]
        query_seq = fraction[3]
        fraction_ans = (query_name, ref_seq, query_seq)
        cesar_fractions.append(fraction_ans)
    return cesar_fractions


def get_starts_ends(r_seq):
    """Return start and end indexes."""
    r_len = len(r_seq)
    start_end = []
    for i in range(1, r_len):
        prev = r_seq[i - 1]
        cur = r_seq[i]
        if prev == " " and cur != " ":
            # exon starts
            start_end.append(i)
            continue
        elif prev != " " and cur == " ":
            start_end.append(i)
            continue
    return start_end


def main():
    """Entry point."""
    args = parse_args()
    with open(args.cesar_output, "r") as f:
        cesar_line = f.read()
    cesar_fractions = read_cesar_out(cesar_line)
    place_holder = " " * (args.flank_size + 1)
    # main loop
    for fraction in cesar_fractions:
        num_to_exons = {}
        q_name = fraction[0]
        # add placeholder to avoid index error
        r_seq = place_holder + fraction[1] + place_holder
        q_seq = place_holder + fraction[2] + place_holder
        starts_ends = get_starts_ends(r_seq)
        for ex_num, start_end in enumerate(parts(starts_ends, 2), 1):
            start, end = start_end
            flanked_r_exon = r_seq[start - args.flank_size: end + args.flank_size]
            flanked_q_exon = q_seq[start - args.flank_size: end + args.flank_size]
            num_to_exons[ex_num] = (flanked_r_exon, flanked_q_exon)
        exon_nums = sorted(num_to_exons.keys())
        print("# query ID == {}".format(q_name))
        for exon_num in exon_nums:
            r_seq, q_seq = num_to_exons[exon_num]
            r_codons, q_codons = [], []
            exon_len = len(r_seq)
            u_in_a_row = 0
            pointer = 0
            r_first_codon = ""
            q_first_codon = ""
            for i in range(1, exon_len):
                prev_ = r_seq[i - 1]
                cur_ = r_seq[i]
                if not prev_.isupper() and cur_.isupper():
                    r_first_codon = r_seq[: i]
                    q_first_codon = q_seq[: i]
                    pointer = i
                    break
            r_codons.append(r_first_codon)
            q_codons.append(q_first_codon)

            for i in range(exon_len):
                r_char = r_seq[i]
                r_upp = r_char.isupper()
                u_in_a_row = u_in_a_row if r_upp is False else u_in_a_row + 1
                if u_in_a_row == 3:
                    r_codon = r_seq[pointer: i + 1]
                    q_codon = q_seq[pointer: i + 1]
                    r_codons.append(r_codon)
                    q_codons.append(q_codon)
                    pointer = i + 1
                    u_in_a_row = 0
            r_last_codon = r_seq[pointer: ]
            q_last_codon = q_seq[pointer: ]
            r_codons.append(r_last_codon)
            q_codons.append(q_last_codon)
            print(">ref_exon_{}\n{}".format(exon_num, " ".join(r_codons)))
            print(">que_exon_{}\n{}".format(exon_num, " ".join(q_codons)))

    sys.exit(0)


if __name__ == "__main__":
    main()