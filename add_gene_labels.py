#!/usr/bin/env python3
"""Post processing tool.

Replaces all the EnsemblIDs with EnsemblID<dot>Gene_name
in any file you want."""
import argparse
import re
import os
import sys

__author__ = "Bogdan Kirilenko, 2018."
PATTERN = r"ENS\w[\d]{11}"  # like ENST00000000000
LOCATION = os.path.dirname(__file__)


def eprint(msg):
    """Like print but for stderr."""
    sys.stderr.write(msg + "\n")


def die(msg, rc=1):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input_file", type=str, help="Input file containing ensembl ids.")
    app.add_argument("--gene_names_table", type=str, default="{}/data/ensGeneIdToName.txt".format(LOCATION),
                     help="Biomart table containing gene names and ensembl ids.")
    app.add_argument("--output_file", default="stdout", help="Output file, default stdout.")
    app.add_argument("--inplace", "-i", action="store_true", dest="inplace",
                     help="Save output in the input file.")
    app.add_argument("--remove_ens_ids", "-r", action="store_true", dest="remove_ens_ids",
                     help="Replace ensemblIDs with gene names instead of appending")
    app.add_argument("--show_none", "-n", action="store_true", dest="show_none",
                     help="Instead of skipping such the cases, assign a name None to the genes "
                     "for what the gene name was not found.")
    app.add_argument("--sep", "-s", default=".", help="Separator between gene ID and name")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_ensembl_data(gene_names_table):
    """Return dict ensemblID to gene_name."""
    with open(gene_names_table, "r") as f:
        ens_data = [x[:-1].split("\t") for x in f.readlines()]
    header, fields = ['Gene stable ID', 'Transcript stable ID', 'Gene name'], ens_data[0]
    try:
        gene_stable_id_index = header.index(fields[0])
        trans_stable_id_index = header.index(fields[1])
        gene_name_index = header.index(fields[2])
    except ValueError:  # absence on one of the necessary fields
        err_msg = "Error! Please make sure that gene_names_table contains the following fields:\n{0}\n"\
                  "You can download the proper input from the EnsemblBiomart".format("\t".join(fields))
        die(err_msg)
    # the file is proper, let's fill the dict
    gene_stable_ids = [x[gene_stable_id_index] for x in ens_data[1:]]
    trans_stable_ids = [x[trans_stable_id_index] for x in ens_data[1:]]
    gene_names = [x[gene_name_index] for x in ens_data[1:]]
    gene_id_to_name, trans_id_to_name = {}, {}
    for i in range(len(gene_names)):
        gene_id_to_name[gene_stable_ids[i]] = gene_names[i]
        trans_id_to_name[trans_stable_ids[i]] = gene_names[i]
    return gene_id_to_name, trans_id_to_name


def main():
    """Entry point."""
    args = parse_args()
    # read Ensembl data
    gene_id_to_name, trans_id_to_name = read_ensembl_data(args.gene_names_table)
    # replace ENST we found with our ids
    output_buffer = []

    f = open(args.input_file, "r") if args.input_file != "stdin" else sys.stdin
    for line in f:
        # split in tokens and find those with ENSX[numbers] inside
        tokens, enss = line.split(), []
        for token in tokens:
            ens_match_span = re.search(PATTERN, token)
            if not ens_match_span:
                continue
            ens_match = token[ens_match_span.start(): ens_match_span.end()]
            enss.append(ens_match)

        for ens in enss:
            enst = trans_id_to_name.get(ens)
            ensg = gene_id_to_name.get(ens)
            not_found = enst is None and ensg is None
            # the main switch
            if not_found and args.show_none and args.remove_ens_ids:
                line = line.replace(ens, "None")
            elif not_found and args.show_none:
                line = line.replace(ens, "{0}{1}{2}".format(ens, args.sep, "None"))
            elif not_found:
                continue
            elif enst and args.remove_ens_ids:
                line = line.replace(ens, enst)
            elif enst:
                line = line.replace(ens, "{0}{1}{2}".format(ens, args.sep, enst))
            elif ensg and args.remove_ens_ids:
                line = line.replace(ens, ensg)
            elif ensg:
                line = line.replace(ens, "{0}{1}{2}".format(ens, args.sep, ensg))
            else:
                continue
            output_buffer.append(line)
    f.close()

    # save the output
    output_file = args.output_file if not args.inplace else args.input_file
    if args.inplace and args.input_file == "stdin":
        # cannot write in stdin!
        output_file = "stdout"
    f = open(output_file, "w") if output_file != "stdout" else sys.stdout
    f.write("".join(output_buffer))

    f.close() if output_file != "stdout" else None
    sys.exit(0)


if __name__ == "__main__":
    main()
