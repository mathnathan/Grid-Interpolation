from readF90binary import *
import numpy as np

f = FortranFile('data/SMinterp.dat')
x = f.readReals()

print "leaving python viz"
