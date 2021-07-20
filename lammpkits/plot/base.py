#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module provides various kinds of tools for plotting.
"""

import numpy as np
from matplotlib import pyplot as plt

DEFAULT_COLORS = ['r', 'b', 'm', 'y', 'g', 'c']
DEFAULT_COLORS.extend(plt.rcParams['axes.prop_cycle'].by_key()['color'])
DEFAULT_MARKERS = ['o', 'v', ',', '^', 'h', 'D', '<', '*', '>', 'd']


def line_chart(ax,
               xdata:np.array,
               ydata:np.array,
               xlabel:str,
               ylabel:str,
               label:str=None,
               sort_by='x',
               **kwargs):
    """
    Plot line chart in ax.

    Args:
        ax: subplot of matplotlib
        xdata (np.array): Input xdata.
        ydata (np.array): Input ydata.
        xlabel (str): x label.
        ylabel (str): y label.
        label (str): Label for ax.
        sort_by (str): if sort_by == 'y', sort by y data.
        kwargs: c, marker, facecolor, s, alpha.

    Notes:
        'kwargs' is parsed to ax.scatter. for more detailed information,
        see documentation for ax.scatter.
    """
    if 'c' in kwargs.keys():
        c = kwargs['c']
    else:
        c_num = len(ax.get_lines()) % len(DEFAULT_COLORS)
        c = DEFAULT_COLORS[c_num]

    if 'marker' in kwargs.keys():
        marker = kwargs['marker']
    else:
        marker_num = len(ax.get_lines()) % len(DEFAULT_MARKERS)
        marker = DEFAULT_MARKERS[marker_num]

    if 's' in kwargs.keys():
        s = kwargs['s']
    else:
        s = None

    if 'facecolor' in kwargs.keys():
        facecolor = kwargs['facecolor']
    else:
        facecolor_num = len(ax.get_lines()) % 2
        if facecolor_num == 0:
            facecolor = 'None'
        else:
            facecolor = c

    if 'alpha' in kwargs.keys():
        alpha = kwargs['alpha']
    else:
        alpha = 1.

    raw = np.array([xdata, ydata])
    if sort_by == 'y':
        idx = np.array(ydata).argsort()
    else:
        idx = np.array(xdata).argsort()
    sort = raw[:,idx]
    ax.plot(sort[0,:], sort[1,:], linestyle='--', linewidth=0.5, c=c,
            alpha=alpha)
    ax.scatter(sort[0,:], sort[1,:], facecolor=facecolor, marker=marker,
               edgecolor=c, alpha=alpha, label=label, s=s)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
