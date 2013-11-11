#!/usr/bin/env python
import sys

first = True
for line in sys.stdin:
    line = line.strip()
    xy_new, current_count = line.split("\t", 1)
    current_count = int(current_count)
    if first:
        xy = xy_new
        count = 1
        first = False
    else:
        if xy != xy_new:
            print '%s,%s' % (xy, str(count))
            xy = xy_new
            count = 1
        else:
            count += current_count
print '%s,%s' % (xy, str(count))