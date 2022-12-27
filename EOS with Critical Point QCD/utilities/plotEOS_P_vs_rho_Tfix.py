#/usr/bin/env python3
# Copyright Chun Shen @ 2018

from numpy import *
import sys
from os import path
from scipy import interpolate
import matplotlib.pyplot as plt

try:
    EOS_table_folder = str(sys.argv[1])
except:
    print("Usage: python3 %s EOS_table_folder" % str(sys.argv[0]))
    exit(1)

hbarC    = 0.19733

filename_string = EOS_table_folder.split("Files_")[1].split("/")[0]

ed_table = loadtxt(path.join(EOS_table_folder,
                             "EnerDens_Final_%s_3D.dat" % filename_string))
nB_table = loadtxt(path.join(EOS_table_folder,
                             "BarDens_Final_%s_3D.dat" % filename_string))
P_table  = loadtxt(path.join(EOS_table_folder,
                             "Press_Final_%s_3D.dat" % filename_string))

n_muB = 451
n_T   = 771
muB = linspace(0.000, 0.450, n_muB)
T   = linspace(0.030, 0.800, n_T)

muB_table, T_table = meshgrid(muB, T)
ed_table = ed_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)
nB_table = nB_table[:, 2].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)
P_table  = P_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)

color_list = ['k', 'r', 'g', 'b', 'm']
fig = plt.figure()
ax  = plt.axes([0.13, 0.13, 0.82, 0.82])
for idx, T_idx in enumerate([120, 118, 115, 113, 110]):
    print("T = %g GeV" % T[T_idx])
    plt.plot(nB_table[T_idx, :]**2., P_table[T_idx, :], linestyle='-',
             color='%s' % color_list[idx], label = "T = %g GeV" % T[T_idx])
plt.legend(loc=0, fontsize=15) 
plt.ylabel("P (GeV/fm$^4$)", fontsize=15)
plt.xlabel("$\\rho_B^2$ (1/fm$^6$)", fontsize=15)
plt.savefig("P_vs_rhobsq_Tfix.pdf", fmt='pdf')

fig = plt.figure()
ax  = plt.axes([0.13, 0.13, 0.82, 0.82])
for idx, T_idx in enumerate([120, 118, 115, 113, 110]):
    print("T = %g GeV" % T[T_idx])
    plt.plot(nB_table[T_idx, :]**2., ed_table[T_idx, :], linestyle='-',
             color='%s' % color_list[idx], label = "T = %g GeV" % T[T_idx])
plt.legend(loc=0, fontsize=15) 
plt.ylabel("e (GeV/fm$^4$)", fontsize=15)
plt.xlabel("$\\rho_B^2$ (1/fm$^6$)", fontsize=15)
plt.savefig("e_vs_rhobsq_Tfix.pdf", fmt='pdf')
