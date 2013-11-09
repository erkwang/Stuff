#!/usr/bin/env python
import sys
import math

#input comes from STDIN (standard input)
for line in sys.stdin:
    line = line.strip()
    x, y = line.split("\t", 1)
    x_low = math.floor(float(x) * 10) / 10
    x_hi = math.ceil(float(x) * 10) / 10
    y_low = math.floor(float(y) * 10) / 10
    y_hi = math.ceil(float(y) * 10) / 10
    print '%s\t%s\t%s\t%s' % (x_low, x_hi, y_low, y_hi)

