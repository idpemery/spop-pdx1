#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluorescence Anisotropy Competition Fitting Script

This script takes an input array containing competitor molecular concentration
values (x) and fluorescence anisotropy values (y) and generates a nonlinear
least squares fit according to Roehrl et al, 2005 and Bouchard et al, 2018 to
extract the dissociation constant of the competitor toward the macromolecule.

NB: this script uses anisotropy as a direct input. To convert polarization to
anisotropy, use the following relation:

A = 2P / (3-P)

User inputs:
-data (text file) (line 38)
-SPOP (macromolecule concentration) (line 45)
-Puc (probe concentration) (line 47)
-Kd1 (Kd of macromolecule and probe, measured directly) (line 47)
-Q (quantum yield) (line 49)

Feel free to use the provided text file (fa_competition_example.txt) to try out
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
#column 2 must contain the POLARIZATION values
data = np.loadtxt('fa_competition_example.txt')

#create array for the x values (peptide concentration) and y values (anisotropy)
pep = data[:,0]
Aobs = data[:,1]

#pull the min and max anisotropy values from dataset
Af = np.min(data[:,1])
Ab = np.max(data[:,1])

#SPOP represents the total concentration of macromolecule in uM assuming it is unchanging over the duration of the experiment
SPOP = 6.0
#Puc is the concentration of fluorophore-labeled probe, which is unchanged over the course of the experiment
Puc = 0.040
#Kd1 is the measured dissociation constant from direct binding assays between SPOP and probe in uM
Kd1 = 2.55
#Q is instrument/wavelength-specific; quantum yield
Q = 3.0

#fitting function 'kd2fit' is defined with input parameters:
    #SPOP = Total concentration of SPOP in uM
    #Kd = dissociation constant based on midpoint of fit (in uM)
    #pep = list of concentration x values
#fitting optimizes Kd2 (competition dissociation constant fitted value) and G (scaling factor)
def kd2fit(pep, Kd2, G):
    d = Kd1 + Kd2 + Puc + pep - SPOP
    e = ((pep - SPOP) * Kd1) + ((Puc - SPOP) * Kd2) + (Kd1 * Kd2)
    f =  -1 * Kd1 * Kd2 * SPOP

    theta_top = (-2* d**3) + (9 * d * e) - (27 * f)
    theta_bottom = 2 * np.sqrt(((d**2)-3*e)**3)
    th = np.arccos(theta_top/theta_bottom)

    top = (2 * np.sqrt(d**2 - 3 * e) * np.cos(th / 3)) - d
    bottom = (3 * Kd1) + (2 * np.sqrt((d**2 - 3 * e)) * np.cos(th / 3)) - d

    FB = (top/bottom)

    Aobs = (((Q * FB * Ab) + (Af *(1 - FB)))/ (1 - (FB * (1 - Q))))

    return G * Aobs

#user-inputted guesses for Kd2 and G
#G is arbitrary scaling factor
guess1 = [8, 1]

#p1 is where the fit is stored
p1, pcov1 = optimize.curve_fit(kd2fit, pep, Aobs, guess1)

fit1 = np.array(kd2fit(pep, p1[0], p1[1]))
Kd2round = round(p1[0], 2)

##PLOTTING##
plt.figure(figsize = (4,4))

plt.ylabel('anisotropy (mA)')
plt.xlabel('competitor concentration (uM)')
#xlim may be changed dependent on concentration regime
plt.xlim(0.01, 10000)
#ylim may be changed based on anisptropy values
plt.ylim(110, 170)

plt.xscale('log')
plt.title('puc-spop competition binding')
#plots Aobs anisotropy data
plt.scatter(pep, Aobs, s = 48, facecolors='none', edgecolors='b', marker = '^' )
plt.annotate('Kd2 = ' +str(Kd2round)+ ' uM', xy = (5, 150))

#plots fit from Kd2 function
plt.plot(pep, fit1 , 'k')


#plt.savefig('puc-spopFA.ps',  format = 'ps', dpi = 600)
