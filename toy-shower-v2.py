#!/usr/bin/env python
# an oversimplified (QED-like) parton shower, 
# for Zuoz lectures (2016) by Gavin P. Salam
# with additional modifications (2017) by Matteo Cacciati and Gavin Salam
#
# - inclusion of running coupling effects, using the rejection method
#   (you can get back to fixed coupling by modifying the alphas initialisation
#   line)
# - printout of additional observables (scalar pt sum, "Higgs" pt == vector pt sum)
from random import random, seed
from math import pi, exp, log, sqrt, cos, sin

ptHigh = 100
ptCut  = 1.0
alphas_mz = 0.12
CA=3

# fix the seed for reproducibility between runs
seed(1)

# some handy global variables to carry around the place
alphas = None
alphas_ptCut = -1 


def main():
    # set up the coupling, including information about its maximum value
    global alphas, alphas_ptCut
    # you can choose whether you want to use a fixed coupling or a running one
    alphas = StrongCoupling(alphas_mz, fixed=False)
    alphas_ptCut = alphas(ptCut)

    # loop over events
    for iev in range(0,100000):
        print "\nEvent", iev
        event()

def event():
    # start with maximum possible value of Sudakov
    sudakov  = 1
    scalar_pt_sum = 0
    higgs_vector_pt = Perp(0,0)
    while (True):
        # scale it by a random number 
        sudakov *= random()
        # deduce the corresponding pt
        pt = ptFromSudakov(sudakov)
        # if pt falls below the cutoff, event is finished
        if (pt < ptCut): break
        # account for running of the coupling with the rejection method 
        acceptance_prob = alphas(pt) / alphas_ptCut
        # if random is larger than the acceptance probability, then we
        # ignore this emission and continue on down
        if (random() > acceptance_prob): continue

        print "  primary emission with pt = ", pt

        scalar_pt_sum += pt
        # now get azimuthal angle
        phi = random() * 2*pi
        # so as to calculate the vector pt sum
        higgs_vector_pt += Perp(pt*cos(phi), pt*sin(phi))
    print "scalar_pt_sum, higgs_vector_pt = ", scalar_pt_sum, higgs_vector_pt.abs()
    
def ptFromSudakov(sudakovValue):
    """Returns the pt value that solves the relation 
       Sudakov = sudakovValue (for 0 < sudakovValue < 1).
       (This is a fixed-coupling Sudakov -- use the rejection method to get a running coupling)
    """
    norm = (2*CA/pi)
    # r = Sudakov = exp(-alphas * norm * L^2)
    # --> log(r) = -alphas * norm * L^2
    # --> L^2 = log(r)/(-alphas*norm)
    L2 = log(sudakovValue)/(-alphas_ptCut * norm)
    pt = ptHigh * exp(-sqrt(L2))
    return pt

#------------------------------------------------------------------------------
class StrongCoupling(object):
    """Contains a one-loop, 5-flavour running coupling (or alternatively a fixed coupling)
    """
    def __init__(self, alphas_mz, fixed=False):
        """Creates the class, initialised with the value of the coupling at the Z mass"""
        self.alphas_mz = alphas_mz
        self.mz = 91.2
        self.b0 = (11*3 - 2*5)/(12*pi)
        self.fixed = fixed
   
    def __call__(self, mu):
        "Returns the value of the coupling at scale mu"
        if (self.fixed): return self.alphas_mz
        else:            return self.alphas_mz / (1 + 2 * self.b0 * self.alphas_mz * log(mu/self.mz))

#------------------------------------------------------------------------------
class Perp(object): 
    """A minimal 2d vector class to store a transverse momentum
    """
    def __init__(self,px,py):
        self.px = px
        self.py = py
   
    def __add__(self, other):
        return Perp(self.px+other.px, self.py+other.py)

    def __iadd__(self, other):
        self.px += other.px
        self.py += other.py
        return self

    def abs(self):
        return sqrt(self.px**2 + self.py**2)
main()

