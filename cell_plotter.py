#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes two text files as inputs: one (speck) contains a list of the
fluorescence intensity values in cells with SPOP in speckles and the other
(nospeck) contains a list of the fluorescence intensity values in cells with
diffuse SPOP. Each data point represents one cell.

The script performs a two-tailed student's t-test (built into scipy) to
determine statistical significance and plots asterisks and raw p values on the
plot. The script also plots the data as a violin plot to show data distribution
and as individual data points for transparency.

@author: emeryusher
gau1@psu.edu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

speck = np.loadtxt('wt_speck.txt')
nospeck = np.loadtxt('wt_nospeck.txt')


data = [speck, nospeck]

#number of cells in each data set
s = len(speck)
n = len(nospeck)

#perform student's t-test
#scipy-calculated p value, is very similar to longhand
t2, p2 = stats.ttest_ind(speck,nospeck)


print('')
print('p = ' +str(p2))

print('number of speckled cells = ' +str(s))
print('number of non-speckled cells = ' +str(n))

#creates violin plot
pos = [1, 1.3]
plt.violinplot(data, pos, widths=0.2, showmeans=True, showextrema=True, showmedians=False, bw_method=0.3)

plt.ylabel('GFP intensity (au)')
plt.title('WT Pdx1 expression and SPOP speckling')
plt.ylim(-1000, 40000)
plt.rcParams["axes.grid"] = False

#conditional to annotate the p value concisely on the plot
if p2 > 0.05:
    plt.annotate('no significance', xy = (1.15, 37500), horizontalalignment='center')

elif 0.05 >= p2 > 0.01:
    plt.annotate('p < 0.05', xy = (1.15, 37500), horizontalalignment='center')

elif 0.01 >= p2 > 0.001:
    plt.annotate('p < 0.01', xy = (1.15, 37500), horizontalalignment='center')

elif 0.001 >= p2 > 0.0001:
    plt.annotate('p < 0.001', xy = (1.15, 37500), horizontalalignment='center')

else:
    plt.annotate('p < 0.0001',  xy = (1.15, 37500), horizontalalignment='center')

#hard coded labels for each data set on plot
plt.annotate('cells w/ speckled SPOP', xy = (1, 10000), horizontalalignment='center')
plt.annotate('cells w/ diffuse SPOP', xy = (1.3, 34000), horizontalalignment='center')

#optional add individual data points into plot
ys = speck
yn = nospeck
xs = np.random.normal(pos[0], 0.01, s)
xn = np.random.normal(pos[1], 0.01, n)

plt.scatter(xn, yn, color = 'navy', s = 24)
plt.scatter(xs, ys, color = 'navy', s = 24)

#save figure as appropriate
#plt.savefig('wt-pdx1_analysis2.png',  format = 'png', dpi = 600)
#plt.savefig('190719wt-pdx1_analysis.ps',  format = 'ps', dpi = 300)
