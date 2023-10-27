#!/usr/bin/env python3
import numpy as np
import matplotlib.cm as cm
import re
import matplotlib.pyplot as plt

# function to compute normalized Legendre function
def plm(theta, degree):
    p_lm = np.zeros((degree + 1, degree + 1))
    u = np.sqrt(1 - np.cos(theta)**2)
    for l in range(degree +1):
        if l == 0:
            p_lm[l, l] = 1
        elif l == 1:
            p_lm[l, l] = u * np.sqrt(3)
        elif l > 1:
            p_lm[l, l] = u * np.sqrt((2*l + 1) / (2*l)) * p_lm[l-1, l-1]

    for l in range(degree + 1):
        for m in range(0, l+1):
            if m == l-1:
                if m == 0:
                    delta = 1
                else:
                    delta = 0
                a_lm = (2/np.sqrt(1+delta)) * ((m+1)/np.sqrt((l-m)*(l+m+1)))
                p_lm[l, m] = a_lm * (np.cos(theta)/u) * p_lm[l, m+1]

    for l in range(degree + 1):
        for m in range(l-2, 0 - 1, -1):
            if m == 0:
                delta = 1
            else:
                delta = 0
            a_lm = (2/np.sqrt(1+delta)) * ((m+1)/np.sqrt((l-m)*(l+m+1)))
            b_lm = (1 / np.sqrt(1 + delta)) * ((np.sqrt((l + m + 2) * (l - m - 1))) / (np.sqrt((l - m) * (l + m + 1))))
            p_lm[l, m] = a_lm * (np.cos(theta)/u) * p_lm[l, m + 1] - b_lm * p_lm[l, m+2]
    return p_lm

# function to read the Stokes coefficients
def readstrokescoefficients(filename):
    with open(filename) as f:
        reach_end_of_head = 0
        for line in f:
            line = line.strip()
            if reach_end_of_head == 0:
                if line.startswith("earth_gravity_constant"):
                    GM = float(line[line.index(" ") + 1:])
                elif "radius" in line:
                    line = re.sub(' +', ' ', line)
                    R = float(line.split(" ")[1])
                elif line.startswith("max_degree"):
                    line = re.sub(' +', ' ', line)
                    max_degree = int(line.split(" ")[1])
                    C = np.zeros((max_degree + 1, max_degree + 1))
                    S = np.zeros((max_degree + 1, max_degree + 1))
                else:
                    if line.startswith("end_of_head"):
                        reach_end_of_head = 1
            else:
                line = re.sub(' +', ' ', line)
                line = line.split()
                L = int(line[1])
                M = int(line[2])
                C[L, M] = float(line[3])
                S[L, M] = float(line[4])
    return C, S, R, GM


