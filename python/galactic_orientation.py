#!/usr/bin/env python3
"""galactic_orientation
Uses cumulative distribution functions (CDF) for portions of a galaxy to
determine the orientation of the galaxy.
"""
from sys import stdin
from argparse import ArgumentError, ArgumentParser, FileType
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import spline

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'), default=stdin,
        help='input file containing x and y coordinates, '
             'as well as magnitude (default is stdin)')
    parser.add_argument('-o', '--output', type=str,
        help='file to save plot to')
    parser.add_argument('-n', '--north', type=float, default=0.0,
        help='lower bound on northern region')
    parser.add_argument('-s', '--south', type=float, default=0.0,
        help='upper bound on southern region')
    parser.add_argument('-e', '--east', type=float, default=0.0,
        help='lower bound on eastern region')
    parser.add_argument('-w', '--west', type=float, default=0.0,
        help='upper bound on western region')
    parser.add_argument('regions', metavar='R', nargs='+',
        help='list of regions to plot cdf for')

    args = parser.parse_args()

    if args.south > args.north:
        raise ArgumentError('upper bound on south exceeds lower bound on north')
    if args.east > args.west:
        raise ArgumentError('upper bound on east exceeds lower bound on west')

    return args

def cdf(region):
    return lambda mag: region[region <= mag].size / region.size

def plot_cdf(axes, region, label):
    if region.size == 0: return

    sorted_region = np.sort(region)
    region_cdf = cdf(sorted_region)
    region_cd  = list(map(region_cdf, sorted_region))
#    smooth_x   = np.linspace(sorted_region.min(), sorted_region.max(), 300)
#    smooth_y   = spline(sorted_region, region_cd, smooth_x)
#    axes.plot(smooth_x, smooth_y, label=label)
    axes.plot(sorted_region, region_cd, label=label)

def main():
    args = get_args()

    x, y, mag = np.loadtxt(args.input, unpack=True)

    west    =               x >= args.west
    vcenter = np.logical_and(x < args.west,
                             x >= args.east)
    east    =               x <  args.east
    north   =               y >= args.north
    hcenter = np.logical_and(y < args.north,
                             y >= args.south)
    south   =               y <  args.south

    regions = dict(
        C  = mag[np.logical_and(vcenter, hcenter)],
        N  = mag[north],
        S  = mag[south],
        E  = mag[east],
        W  = mag[west],
        NC = mag[np.logical_and(vcenter, north)],
        SC = mag[np.logical_and(vcenter, south)],
        EC = mag[np.logical_and(hcenter, east)],
        WC = mag[np.logical_and(hcenter, west)],
        NE = mag[np.logical_and(north,   east)],
        NW = mag[np.logical_and(north,   west)],
        SE = mag[np.logical_and(south,   east)],
        SW = mag[np.logical_and(south,   west)]
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    ax.invert_xaxis()

    for region in args.regions:
        plot_cdf(ax, regions[region], region)

    ax.set_xlabel('M')
    ax.set_ylabel('P(M >= m)')

    ax.legend()
    fig.savefig(args.output)

if __name__ == '__main__':
    exit(main())
