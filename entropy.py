#!/usr/bin/env python

# Calculate the entropy given the base and the probabilities.
from math import log
import sys

usage = "Usage: python " + sys.argv[0] + " base p1 p2 p3...\n"

if len(sys.argv) <= 1:
    sys.stderr.write(usage)
    sys.exit(1)

base = float(sys.argv[1])
prob = [float(i) for i in sys.argv[2:]]

ent = 0.0
for p in prob:
    ent += p * log(p, base)
ent = -ent

print ent
