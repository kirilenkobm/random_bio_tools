#!/sw/bin/python3
"""Check codon alignment quality."""
import argparse
import sys
from collections import Counter, defaultdict
import operator

__author__ = "Bogdan Kirilenko, 2018"
# defaults
LETTERS = ["A", "T", "G", "C", "N", "-"]
STATES = [-2, -1, 0, 1, 2]


def eprint(msg):
    """Write to stderr."""
    sys.stderr.write(msg + "\n")


def die(msg, rc=1):
    """Interrupt program."""
    eprint(msg)
    eprint("Program finished with exit code {}".format(rc))
    sys.exit(rc)


def parse_args():
    """Read args, check (if possible)."""
    # read arguments
    app = argparse.ArgumentParser("A tool for codon alignment quality estimation.")
    app.add_argument("input", type=str, help="Input fasta file containg aligned sequences.")
    app.add_argument("species", type=str, help="Species of interest, like 'mm10'. Single or comma-separated list.")
    app.add_argument("output", type=str, help="Save result in...")
    app.add_argument("--window_size", "-w", type=int, default=7, help="Window size, 7"
                     " is recommended (default), must be > 2.")
    app.add_argument("--threshold", type=float, default=2.0, help="Score considered as significant.")
    app.add_argument("-a", "--append", action="store_true", dest="append",
                     help="Use to append the results in a file that already exists.")
    args = app.parse_args()
    # check if everything is alright
    if len(sys.argv) < 3:  # there are no arguments
        app.print_help()
        sys.exit(0)
    if args.window_size < 3:  # it makes no sence to use window of this size
        die("Error! Expected window size  > 2, {} given.".format(args.window_size))
    if args.window_size <= args.threshold:
        die("Error! Requested threshold {0} is unreachable! "
            "It should be less than window size {1}".format(args.threshold, args.window_size))
    elif args.threshold < 0:
        die("Error! Threshold should be > 0, {0} given.".format(args.threshold))
    return args


def read_fasta(fasta_file):
    """Read fasta, return sequences."""
    # open the file
    with open(fasta_file, "r") as f:
        fasta_data = f.read().split(">")

    assert fasta_data[0] == ""  # if assertion is failed - it is not a fasta
    sequences = {}  # save results here
    order = []  # in case if we are interested in the order

    # read line by line
    for elem in fasta_data[1:]:  # the first elem is ""
        raw_lines = elem.split("\n")  # might be useful if there is fasta-80 or fasta-60
        header = raw_lines[0]  # >header\nATGAGAGAGTTAC --> header, ATGACACATACGA
        lines = [x for x in raw_lines[1:] if x != ""]  # the last elem is "" and the first is the header
        sequence = "".join(lines)
        # add data to collectors
        sequences[header] = sequence
        order.append(header)

    # check if data is correct
    if len(sequences) == 0:  # either empty or the file is corrupted
        die("Error, there are no fasta-formatted sequences in {}!".format(fasta_file))

    if len(sequences.keys()) != len(order):  # means there are non-unique headers
        err_msg = "Error! Sequence names must be unique! There are " \
                  "{0} sequences with {1} unique names!".format(len(sequences.keys()), len(order))
        die(err_msg)

    return sequences, order


def compute_frequences(sequences, seq_len):
    """Return dict position: frequences for bases."""
    frequences = {}
    sequences_num = len(sequences.keys())
    columns = defaultdict(list)  # revert sp: sequence to number - bases

    for k, seq in sequences.items():
        if len(seq) != seq_len:      # check that sequences are equal in length
            err_msg = "Error! Sequences must have the equal lenght! Make sure that you " \
                      "operate with aligned sequences! {0} has {1} bases but your sp has {2}".format(k, len(seq), seq_len)
            die(err_msg)

        # they are equal --> let's work
        for num, base in enumerate(seq):  # num of base and base
            columns[num].append(base)  # num: [base base base]

    # convert to frequences
    for num, column in columns.items():
        quantities = Counter(column)  # count frequences of As, Ts etc
        frequences[num] = {k: v / sequences_num for k, v in quantities.items()}  # absolute numbers to relatives
    # consider special cases
    for cond in [-2, -1, seq_len, seq_len + 1]:
        frequences[cond] = {}  # in case we will call frequences [-2]
    return frequences


def split_windows(sp_seq, window_size):
    """Split a sequence in a set of subsequences according the window size."""
    if window_size > len(sp_seq):  # in case if sequence is too short (window is too huge) it makes no sence
        die("Error! Window size {0} is bigger than sequence len {1}!".format(window_size, len(sp_seq)))
    return [sp_seq[i:i+window_size] for i in range(len(sp_seq) - window_size + 1)]


def get_scores(frequences, windows):
    """Return score for each shift for each window."""
    window_scores = {}
    # main loop
    for w_num, window in enumerate(windows):  # num of the first base
        # define start point in absolute coordinates
        window_scores[w_num] = {}  # collect shift: score for each
        # test each possible shift between -2 -1 0 1 2
        for shift in STATES:
            shift_score = 0  # each shift has the own score
            # check score of each position
            for b_num, base in enumerate(window):
                freq_pos = w_num + b_num + shift  # start of the window + number of base in window + shift
                freq_data = frequences[freq_pos]  # get corresponding dict
                base_score = freq_data.get(base) if freq_data.get(base) else -0.25  # in case of None
                shift_score += base_score  # add base score to shift score
            window_scores[w_num][shift] = shift_score
            window_scores[w_num]["seq"] = window  # save sequence just in case
    return window_scores


def check_scores(window_scores, threshold):
    """Return positions: shift in case if shift is not the best."""
    suspect = {}
    for w_num, scores in window_scores.items():
        # get sequence and remove it
        sequence = scores["seq"]
        del scores["seq"]
        # get the maximal score
        shift, score = max(scores.items(), key=operator.itemgetter(1))
        # if shift != zero --> alignment is broken here
        # and also there are some conditions to check
        zero_condition = shift != 0  # it should be not shift 0
        thr_condirion = score > threshold  # better then defined threshold
        sign_condition = score > scores[0] * 1.5  # and bigger than 0
        # are all these conditions true?
        all_conditions = zero_condition and thr_condirion and sign_condition
        if all_conditions:  # and if yes append this stuff
            suspect[w_num] = (shift, score, sequence, scores[0])
    return suspect


def save_result(broken_places, window_size, output, append, sp):
    """Save the result."""
    output_line = "sp\tpositions\tseq\tshift\tscore\tzero\n"  # collect result there
    for start_point, shift_score_seq_zer in broken_places.items():
        shift, score, seq, zero = shift_score_seq_zer  # it is a tuple (shift, score, sequence, zero)
        end_point = start_point + window_size  # user-friendly output
        position = "{0}-{1}".format(start_point + 1, end_point)
        new_line = "{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5:.4f}\n".format(sp, position, seq, shift, score, zero)
        output_line += new_line  # append this line
    # define the mode of file object according -a param
    m = "a" if append else "w"  # a - append to existing file, w - create a new file
    if m == "a" and output == "stdout":  # if stdout there isn't file to append
        eprint("Warning! Param -a makes no sence if write to stdout.")
    # in case of stdout we should initiate object f in an alternative way
    f = sys.stdout if output == "stdout" else open(output, m)
    # write it
    f.write(output_line)
    f.close()


def main():
    """Entry point."""
    args = parse_args()  # load args
    sequences, order = read_fasta(args.input)  # read fasta
    if args.species not in order:  # check if we can use these species
        err_msg = "Error! There is no sequence with name {0} " \
                  "The possible options are:\n{1}".format(args.species, " ".join(sorted(order)))
        die(err_msg)
    sp_seq = sequences[args.species]  # pick the sequence of interest
    del sequences[args.species]  # don't need it there
    # the main part of program
    frequences = compute_frequences(sequences, len(sp_seq))  # compute frequences for each column
    windows = split_windows(sp_seq, args.window_size)  # get all the sequences with the certain window size
    # compute all the possible scores
    window_scores = get_scores(frequences, windows)
    # filter scores if zero-shift is not the best one
    broken_places = check_scores(window_scores, args.threshold)
    # save it wherever we need it
    save_result(broken_places, args.window_size, args.output, args.append, args.species)
    sys.exit(0)  # say good bye


if __name__ == "__main__":
    main()
