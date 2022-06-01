"""Homework 1: Distance functions on vectors.

Homework 1: Distance functions on vectors
Course    : Data Mining (636-0018-00L)

Auxiliary functions.

This file implements the distance functions that are invoked from the main
program.
"""
# Author: Damian Roqueiro <damian.roqueiro@bsse.ethz.ch>
# Author: Bastian Rieck <bastian.rieck@bsse.ethz.ch>

#import numpy as np
import math
from math import sqrt

def manhattan_dist(v1, v2):
    d = 0
    length = len(v1)
    for i in range(length):
        d = float(d+abs(v1[i]-v2[i]))
    return d

def hamming_dist(v1, v2):
    dist = 0
    length = len(v1)
    for i in range (length):
        number1 = ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', v1[i]))
        number2 = ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', v2[i]))
        count = 0
        for x in range(32):
            if number1[x] != number2[x]:
                count = count+1
        if count>0:
            dist=dist+1
    return dist


def euclidean_dist(v1, v2):
    dist = 0
    length = len(v1)
    for i in range(length):
        dist = dist+(v1[i]-v2[i])**2
    dist = sqrt(dist)
    return dist


def chebyshev_dist(v1, v2):
    dist = 0
    length = len(v1)
    for i in range(length):
        dist_new = sqrt((v1[i]-v2[i])**2)
        if dist_new > dist:
            dist = dist_new
    return dist

def minkowski_dist(v1, v2, d):
    dist = 0
    length = len(v1)
    for i in range(length):
        dist = dist+abs(v1[i]-v2[i])**d
    dist = dist**(1/d)
    return dist
