#!/usr/bin/env python
# an oversimplified (QED-like) parton shower
# for Zuoz lectures (2016) by Gavin P. Salam
from random import random
from math import pi, exp, log, sqrt

ptHigh = 100.0
ptCut  = 1.0
alphas = 0.12
CA=3

def main():
    for iev in range(0,10):
        print "\nEvent", iev
        event()

def event():
    # start with maximum possible value of Sudakov
    sudakov  = 1
    while (True):
        # scale it by a random number 
        sudakov *= random()
        # deduce the corresponding pt
        pt = ptFromSudakov(sudakov)
        # if pt falls below the cutoff, event is finished
        if (pt < ptCut): break
        print "  primary emission with pt = ", pt

def ptFromSudakov(sudakovValue):
    """Returns the pt value that solves the relation 
       Sudakov = sudakovValue (for 0 < sudakovValue < 1)
    """
    norm = (2*CA/pi)
    # r = Sudakov = exp(-alphas * norm * L^2)
    # --> log(r) = -alphas * norm * L^2
    # --> L^2 = log(r)/(-alphas*norm)
    L2 = log(sudakovValue)/(-alphas * norm)
    pt = ptHigh * exp(-sqrt(L2))
    return pt
    
main()

