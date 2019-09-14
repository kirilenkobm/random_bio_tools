#!/sw/bin/python3
"""Find all intersections for chain X bed.

Writes to stdout the following table:
chain_id<tab>comma-separated list of ovrelapped genes.
"""
import sys
from collections import defaultdict
import subprocess

__author__ = "Bogdan Kirilenko, 2018."


def parse_chain(chain):
    """Return chrom: ranges from chain file."""
    # I need the headers only
    chrom_range = defaultdict(list)
    cmd = "cat {0} | grep chain".format(chain)
    headers = subprocess.check_output(cmd, shell=True).decode("utf-8")
    for header in headers.split("\n"):
        if len(header) == 0:
            continue
        header_info = header.split()
        chrom = header_info[2]
        start = int(header_info[5])
        end = int(header_info[6])
        chain_id = header_info[12]
        chrom_range[chrom].append((chain_id, start, end))
    return chrom_range


def parse_bed(bed):
    """Return chrom: ranges from bed file."""
    chrom_range = defaultdict(list)
    f = open(bed, "r")
    for line in f:
        line_info = line.split("\t")
        chrom = line_info[0]
        start = int(line_info[1])
        end = int(line_info[2])
        gene = line_info[3]
        chrom_range[chrom].append((gene, start, end))
    f.close()
    return chrom_range


def intersect(range_1, range_2):
    """Return intersection size."""
    return min(range_1[2], range_2[2]) - max(range_1[1], range_2[1])


def find_first(chains, beds):
    """Find indexes for the first intersection."""
    if intersect(chains[0], beds[0]) > 0:  # no need to search
        return 0, 0
    first_bed_start, first_bed_end = beds[0][1], beds[0][2]
    first_chain_start, first_chain_end = chains[0][1], chains[0][2]
    if first_chain_end < first_bed_start:
        # we have lots of chains in the beginning not intersecting beds
        for i in range(len(chains)):
            if intersect(chains[i], beds[0]) > 0:
                return i, 0
            elif chains[i][1] > beds[0][2]:
                return i, 1
    else:  # lots of beds in the beginning not covered by chains
        for i in range(len(beds)):
            if intersect(chains[i], beds[0]) > 0:
                return 0, i
            elif beds[i][1] > chains[0][2]:
                return 1, i


def overlap(chains, beds):
    """Return intersections for chain: bed."""
    # init state, find FIRST bed intersecting the FIRST chain
    chain_beds = defaultdict(list)
    chain_init, start_with = find_first(chains, beds)
    chains_num = len(chains)
    bed_len = len(beds)
    for i in range(chain_init, chains_num):
        FLAG = False  # was intersection or not?
        FIRST = True
        chain = chains[i]
        while True:
            if FIRST:  # start with previous start, first iteration
                bed_num = start_with
                FIRST = False  # guarantee that this condition works ONCE per loop
            else:  # just increase the pointer
                bed_num += 1  # to avoid inf loop

            if bed_num >= bed_len:
                break  # beds are over
            # pick the bed range
            bed = beds[bed_num]

            if chain[2] < bed[1]:  # too late
                break  # means that bed is "righter" than chain

            if intersect(chain, bed) > 0:
                if not FLAG:  # the FIRST intersection of this chain
                    start_with = bed_num  # guarantee that I will assign to starts with
                    # only the FIRST intersection (if it took place)
                    FLAG = True  # otherwise starts with will be preserved
                # save the intersection
                chain_beds[chain[0]].append(bed[0])

            else:  # we recorded all the region with intersections
                if chain[1] > bed[2]:  # too early
                    # in case like:
                    # gene A EEEEE----------------------------------------EEEEEE #
                    # gene B               EEEEEEEEEE                            #
                    # gene C                               EEEEEEEEE             #
                    # chain                                    ccccc             #
                    # at gene A I will get FLAG = True and NO intersection with gene B
                    # --> I will miss gene C in this case without this condition.
                    continue

                elif FLAG:  # this is not a nested gene
                    break  # and all intersections are saved --> proceed to the next chain

    return chain_beds


def chain_bed_intersect(chain, bed):
    """Entry point."""
    # get list of chrom: ranges for both
    chain_data = parse_chain(chain)
    bed_data = parse_bed(bed)
    chroms = list(bed_data.keys())
    chain_bed_dict = {}  # out answer
    # main loop
    for chrom in chroms:
        # sort the ranges
        bed_ranges = sorted(bed_data[chrom], key=lambda x: x[1])
        chain_ranges = sorted(chain_data[chrom], key=lambda x: x[1])
        chrom_chain_beds = overlap(chain_ranges, bed_ranges)
        chain_bed_dict.update(chrom_chain_beds)
    return chain_bed_dict


def save(dct, output="stdout"):
    """Save output in the file given."""
    f = open(output, "w") if output != "stdout" else sys.stdout
    for k, v in dct.items():
        f.write("{0}\t{1}\n".format(k, ",".join(v) + ","))
    f.close()


if __name__ == "__main__":
    try:  # read args
        chain_file = sys.argv[1]
        bed_file = sys.argv[2]
    except IndexError:
        sys.stderr.write("Usage: {} [chain_file] [bed file]\n".format(sys.argv[0]))
        sys.stderr.write("Output goes to stdout.\n")
        sys.exit(0)
    chain_bed_dict = chain_bed_intersect(chain_file, bed_file)
    save(chain_bed_dict)
