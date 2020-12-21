#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluorescence Anisotropy Direct Fitting Script

This script takes an input array containing macromolecule concentration
values (x) and fluorescence polarization values (y) and generates a nonlinear
least squares fit according to Roehrl et al, 2005 to extract a dissociation
constant.

NB: this script uses anisotropy as a direct input. To convert polarization to
anisotropy, use the following relation:

A = 2P / (3-P)

User inputs:
-data (text file) (line 37)
-Puc  (fluorescent probe concentration) (line 4)
-Kd1 (Kd of macromolecule and probe, measured directly) (line 46)
-Q (quantum yield) (line 48)

Feel free to use the provided text file (fa_direct_example.txt) to try out
this script.

@author: emeryusher
gau1@psu.edu
"""

#imports relevant modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

#indicate the text file containing data
#column 1 must contain the concetration values in uM
#column 2 must contain the ANISOTROPY values
data = np.loadtxt('fa_direct_example.txt')

#create separate arrays for the X and Y data points
pep = data[:,0]
I = data[:,1]

#pull the min and max anisotropy values from dataset
Af = np.min(data[:,1])
Ab = np.max(data[:,1])

#Puc represents the total concentration of fluorophore in uM
Puc = 0.010
#Q is instrument/wavelength-specific; quantum yield
Q = 1

#fitting function 'kdfit' is defined with input parameters:
    #SPOP = Total concentration of SPOP in uM
    #Kd = dissociation constant based on midpoint of fit (in uM)
    #pep = list of concentration x values
#fitting optimizes Kd (direct dissociation constant fitted value) and G (scaling factor)
def kdfit(pep, Kd, G):

    FB = (((Kd + pep + Puc) - np.sqrt(((Kd + pep + Puc)**2) - 4 * pep * Puc)) / (2 * Puc))

    Aobs = ((Q * FB * Ab) + ((Af *(1 - FB))/ (1 - (FB * (1 - Q)))))

    return G * Aobs

#user-inputted guesses for Kd and G
guess1 = [1, 8]

#p1 is where the fit is stored
p1, pcov1 = optimize.curve_fit(kdfit, pep, I, guess1)

fit1 = np.array(kdfit(pep, p1[0], p1[1]))
Kdround = round(p1[0], 2)

##PLOTTING##
plt.figure(figsize = (4,4))

#plotting
plt.ylabel('anisotropy (mA)')
plt.xlabel('SPOP-MATH concentration (uM)')

plt.xscale('log')
plt.title('fPuc-spop direct binding')

#plots Aobs anisotropy data
plt.scatter(pep, I, s = 48, facecolors='none', edgecolors='r' )

#plots fit from Kd function
plt.plot(pep, fit1 , 'k')
plt.annotate('Kd = ' +str(Kdround)+ ' uM', xy = (5, 150))

#plt.savefig('spop-fPUcFA.ps',  format = 'ps', dpi = 600)
