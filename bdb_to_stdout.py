#!/usr/bin/env python3
"""Sample."""
import argparse
import sys
import bsddb3

__author__ = "Bogdan Kirilenko, 2018."


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
    app.add_argument("bdb_file")
    app.add_argument("query") if len(sys.argv) > 2 else app.add_argument("--query")
    # app.add_argument("--k_num", "-k", type=int, default=0)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def handle_db(db_file):
    """Load db."""
    try:
        db = bsddb3.btopen(db_file, "r")
    except Exception:
        die("Cannon open {0}".format(db_file))
    return db


def get_value(db, query_str):
    """Load a value according the key."""
    key = query_str.encode()
    try:
        value = db[key].decode("utf-8")
        db.close()
        return value
    except KeyError:
        db.close()
        die("Cannon find {0} in the file.".format(query_str))
        return None


def db_keys(db):
    """Show keys."""
    keys = [k.decode("utf-8") for k in db.keys()]
    db.close()
    output = "In the db keys are:\n"
    output += ",".join(keys)
    return output


def main():
    """Entry point."""
    args = parse_args()
    db = handle_db(args.bdb_file)
    result = get_value(db, args.query) if args.query else db_keys(db)
    sys.stdout.write(result + "\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
