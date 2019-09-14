#!/usr/bin/env python3
"""Reorder muscle HTML file."""
import argparse
import sys

__author__ = "Bogdan Kirilenko, 2019."
HTML_START = '<HTML>\n<BODY BGCOLOR="#FFEEE0">\n<PRE>'
HTML_END = '</SPAN>\n</PRE>\n</BODY>\n</HTML>'


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
    app.add_argument("html_in", help="Html file produced by MUSCLE")
    app.add_argument("order_key", help="Ordered list of species or fasta")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def make_sort_key(order_file):
    """Generate key for sorting."""
    with open(order_file, "r") as f:
        order_file_content = f.readlines()
    # check if it's fasta
    starts_with_more = sum([1 for l in order_file_content if l.startswith(">")])
    prop_of_more = starts_with_more / len(order_file_content)
    # half or less of lines start with >
    is_fasta = True if prop_of_more <= 0.51 else False
    if not is_fasta:
        species = [l[:-1] for l in order_file_content]
    else:
        species = [l[1:-1] for l in order_file_content if l.startswith(">")]
    spec_order = {s: n for n, s in enumerate(species)}
    return spec_order


def main():
    """Entry point."""
    args = parse_args()
    sort_key = make_sort_key(args.order_key)
    with open(args.html_in, "r") as f:
        html_data = f.read()
    html_data = html_data.replace(HTML_START, "").replace(HTML_END, "")
    html_blocks = html_data.split("\n\n")
    reordered_blocks = []
    for block in html_blocks:
        block_lines = [x for x in block.split("\n") if x != ""]
        spec_line = []
        for line in block_lines:
            species = line.split()[1].split(">")[-1]
            spec_line.append((species, line))
        sorted_lines = sorted(spec_line, key=lambda x: sort_key.get(x[0]))
        sorted_block = "\n".join([x[1] for x in sorted_lines])
        reordered_blocks.append(sorted_block)
    reordered_blocks_str = "\n\n".join(reordered_blocks)
    reordered_html_line = HTML_START + reordered_blocks_str + HTML_END
    print(reordered_html_line)
    sys.exit(0)


if __name__ == "__main__":
    main()
