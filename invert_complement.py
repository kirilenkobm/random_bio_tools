#!/usr/bin/env python3
"""Make reverse complement sequence."""
import sys

try:
    seq = sys.argv[1]
except IndexError:
    sys.exit("Usage: {} [DNA sequence]".format(sys.argv[0]))

complemenletters = {"A": "T", "T": "A", "G": "C", "C": "G", "-": "-", "N": "N"}
reverse = seq[::-1]
reverse_complement = "".join([complemenletters.get(c) if complemenletters.get(c) else "X" for c in reverse])
print(reverse_complement)

