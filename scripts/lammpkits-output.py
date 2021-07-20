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
        '-k', '--keys', type=str, default=None,
        help="Lammps keys shown in log.lammps.")

    parser.add_argument(
        '-p', '--plot', action='store_true',
        help="Plot fig.")

    return parser.parse_args()


def main(logfile,
         keys,
         is_plot,
         ):
    lmp_out = LammpsOutput(logfile=logfile)
    print("Load: %s" % logfile)
    print("Available keys are: {}".format(lmp_out.get_keys()))

    num_of_plots = len(keys)
    if num_of_plots == 1:
        l, k = 1, 1
    elif num_of_plots == 2:
        l, k = 2, 1
    elif num_of_plots <= 4:
        l, k = 2, 2
    elif num_of_plots <= 6:
        l, k = 3, 2
    elif num_of_plots <= 9:
        l, k = 3, 3
    elif num_of_plots <= 12:
        l, k = 4, 3
    elif num_of_plots <= 16:
        l, k = 4, 4
    else:
        raise RuntimeError("Too many keys (more than 16).")

    fig = plt.figure(figsize=(5*l,5*k))
    if is_plot:
        for i, key in enumerate(keys):
            ax = fig.add_subplot(l,k,i+1)
            ax = lmp_out.plot_transition(ax=ax, key=key)
        plt.show()


if __name__ == '__main__':
    args = get_argparse()
    if args.keys is not None:
        keys = list(map(str, args.keys.split()))
    else:
        keys = None

    main(logfile=args.logfile,
         keys=keys,
         is_plot=args.plot,
         )
