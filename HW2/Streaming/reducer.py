#!/usr/bin/env python
import sys

current_count = 1
first = True
result = []
for line in sys.stdin:
    line = line.strip()
    xy_new = line
    if first:
        xy = xy_new
        first = False
    if xy != xy_new:
        print '%s,%s' % (",".join(xy.split("\t")), str(current_count))
        xy = xy_new
        current_count = 1
    else:
        current_count += 1
print '%s,%s' % (",".join(xy.split("\t")), str(current_count))