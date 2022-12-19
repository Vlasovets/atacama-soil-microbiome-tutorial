import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def plot_heatmap(cov, precision, mask, low_rank=None, low=False):
    
    meta_ticks = np.array(cov.columns[-14:]) # TO DO: change according to the data
    bug_ticks = np.arange(len(cov.columns[:-14]))
    ticks = np.hstack((bug_ticks, meta_ticks))

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9  # the right side of the subplots of the figure
    bottom = 0.1  # the bottom of the subplots of the figure
    top = 0.9  # the top of the subplots of the figure
    wspace = -0.6  # the amount of width reserved for blank space between subplots,
    hspace = 0.5  # the amount of height reserved for white space between subplots,
    fontsize = 56
    cmap = "coolwarm"
    vmin = -0.5
    vmax = 0.5
    linewidth = .5
    square = True
    cbar = False

    if low:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(90, 35))

            plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

            ax1.get_shared_y_axes().join(ax2)
            ax3.get_shared_y_axes().join(ax4)

            g1 = sns.heatmap(cov, cmap=cmap, cbar=cbar, ax=ax1, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g1.set_ylabel('')
            g1.set_xlabel('Covariance', fontsize=fontsize)

            g2 = sns.heatmap(precision, cmap=cmap, cbar=cbar, ax=ax2, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g2.set_ylabel('')
            g2.set_xlabel('Inverse covariance', fontsize=fontsize)
            g2.set_yticks([])

            g3 = sns.heatmap(low_rank, cmap=cmap, ax=ax3, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g3.set_ylabel('')
            g3.set_xlabel('Low-rank solution', fontsize=fontsize)
            g3.set_yticks([])

            g4 = sns.heatmap(mask, cmap=cmap, ax=ax4, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g4.set_ylabel('')
            g4.set_xlabel('Mask', fontsize=fontsize)
            g4.set_yticks([])
    else:

        wspace = 0.5  # the amount of width reserved for blank space between subplots,
        hspace = 0.5

        fig, (ax1, ax2 ,ax3) = plt.subplots(1, 3, figsize=(90, 35))

        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        ax1.get_shared_y_axes().join(ax2, ax3)

        g1 = sns.heatmap(cov, cmap=cmap, cbar=cbar, ax=ax1, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g1.set_ylabel('')
        g1.set_xlabel('Covariance', fontsize=fontsize)

        g2 = sns.heatmap(precision, cmap=cmap, cbar=cbar, ax=ax2, vmin=vmin, vmax=vmax, linewidth=linewidth,
                         square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g2.set_ylabel('')
        g2.set_xlabel('Inverse covariance', fontsize=fontsize)
        g2.set_yticks([])

        g3 = sns.heatmap(mask, cmap=cmap, ax=ax3, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g3.set_ylabel('')
        g3.set_xlabel('Mask', fontsize=fontsize)
        g3.set_yticks([])

    return fig