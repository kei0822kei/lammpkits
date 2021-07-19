#!/usr/bin/env python


"""
This script deals with lammps output.
"""

from pprint import pprint
import argparse
from matplotlib import pyplot as plt
from lammpkits.lammps.output import LammpsOutput


# argparse
def get_argparse():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-l', '--logfile', type=str, default=None,
        help="Output log.lammps file.")

    parser.add_argument(
        '-k', '--key', type=str, default=None,
        help="Lammps key shown in log.lammps.")

    parser.add_argument(
        '-p', '--plot', action='store_true',
        help="Plot fig.")

    return parser.parse_args()


def main(logfile,
         key,
         is_plot,
         ):
    lmp_out = LammpsOutput(logfile=logfile)
    print("Load: %s" % logfile)
    print("Available keys are: {}".format(lmp_out.get_keys()))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    if is_plot:
        ax = lmp_out.plot_transition(ax=ax, key=key)
        plt.show()


if __name__ == '__main__':
    args = get_argparse()
    main(logfile=args.logfile,
         key=args.key,
         is_plot=args.plot,
         )
