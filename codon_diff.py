#!/sw/bin/python3
"""Difference between two sequences."""
import argparse
import sys
from evolve import get_alts

__author__ = "Bogdan Kirilenko, 2018."

# genetic code
AA_CODE = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
           "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
           "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
           "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
           "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
           "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
           "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
           "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
           "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
           "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
           "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
           "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
           "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
           "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
           "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
           "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
           "---": "-", "NNN": "X"}
# combinations to take care of
MASKED = ["NNN", "---", "TAG", "TGA", "TAA"]


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
    app.add_argument("first_fasta", type=str, help="First fasta file.")
    app.add_argument("first_sp", type=str, help="Species from the first file.")
    app.add_argument("second_fasta", type=str, help="First fasta file. Write - if the same with the first.")
    app.add_argument("second_sp", type=str, help="Species from the first file. "
                     "Write - to use the same species with the first.")
    app.add_argument("--no_mask", action="store_true", dest="no_mask",
                     help="Do not ignore NNN's and ---'s.")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_fasta(fasta_file):
    """Read fasta, return dict and type."""
    # open the file
    with open(fasta_file, "r") as f:
        fasta_data = f.read().split(">")
    assert fasta_data[0] == ""  # if a file starts with > it should be empty
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
        die("There are not fasta-formatted sequences in {0}!".format(fasta_file))
    if len(sequences.keys()) != len(order):  # it is possible in case of non-unique headers
        err = "Error! Sequences names must be unique! There are" \
              " {0} sequences and {1} unique names!".format(len(sequences.keys()), len(order))
        die(err)
    return sequences, order


def parts(lst, n=25):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def codons_dist(codon1, codon2):
    """Number of different bases in codon."""
    return sum([1 for i in range(3) if codon1[i] != codon2[i]])


def main():
    """Entry point."""
    args = parse_args()
    # read fastas, the first one
    first_seqs, first_species = read_fasta(args.first_fasta)
    if args.first_sp not in first_species:
        die("Error! There is no {0} in the {1}".format(args.first_sp, args.first_fasta))
    # and the second one
    second_fasta = args.second_fasta if args.second_fasta != "-" else args.first_fasta
    second_sp = args.second_sp if args.second_sp != "-" else args.first_sp
    second_seqs, second_species = read_fasta(second_fasta)
    if second_sp not in second_species:
        die("Error! There is no {0} in the {1}".format(args.second_sp, second_fasta))
    # get sequences
    first_seq = first_seqs[args.first_sp]
    second_seq = second_seqs[second_sp]
    # split into codons and check
    first_codons = parts(first_seq, n=3)
    second_codons = parts(second_seq, n=3)
    if len(first_codons[-1]) != len(second_codons[-1]) != 3:
        die("Error! Codon alignment required!")
    if len(first_codons) != len(second_codons):
        die("Error! Sequences of the same length are required!")
    codons_num = len(first_codons)
    # initial values
    diff_number, changes = 1, 0
    syns, non_syns = 0, 0
    output_line = "{0}| codon num {1}; nucl: {2} -> {3}; AA: {4} -> {5}; {6}\n"
    syn_codons, nsyn_codons = 0, 0

    # loop by codons
    for codon_num in range(codons_num):
        first_codon = first_codons[codon_num]
        second_codon = second_codons[codon_num]

        syn_codons_seqs, nsyn_codons_seqs = get_alts(first_codon)
        syn_num, nsyn_num = len(syn_codons_seqs), len(nsyn_codons_seqs)
        syn_codons += syn_num
        nsyn_codons += nsyn_num
        # same codons, ignore
        if first_codon == second_codon:
            continue
        # define what's with masking
        if first_codon in MASKED and not args.no_mask:
            continue  # default value, it not --no_mask
        elif second_codon in MASKED and not args.no_mask:
            continue
        # write difference
        f_AA = AA_CODE.get(first_codon)
        s_AA = AA_CODE.get(second_codon)
        changes += codons_dist(first_codon, second_codon)
        if f_AA == s_AA:
            syns += 1
            synline = "syn"
        else:
            non_syns += 1
            synline = "non-syn"
        # output it
        sys.stdout.write(output_line.format(diff_number, codon_num,
                                            first_codon, second_codon,
                                            f_AA, s_AA, synline))
        diff_number += 1  # number of differencies
    omega_base = nsyn_codons / syn_codons if syn_codons != 0 else 9999
    non_syn_to_syn = non_syns / syns if syns != 0 else 9999
    omega = non_syn_to_syn / omega_base if syns != 0 else 9999
    sys.stdout.write("Overall {0} different codons and {1} changes; {2} synonymous and {3} non-synonymous;\n" \
                     "There are {4} synonymous and {5} non-synonymous sites, ratio is {6}\n" \
                     "Omega: {7}\n".format(diff_number - 1, changes, syns, non_syns, syn_codons,
                                           nsyn_codons, omega_base, omega))
    sys.exit(0)


if __name__ == "__main__":
    main()
