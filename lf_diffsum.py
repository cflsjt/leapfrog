#!/usr/bin/env python

import sys
import collections

counts = collections.defaultdict(list)
total = 0
simmat = collections.defaultdict(int)
libs = set()
with open(sys.argv[1]) as F:
    for line in F:
        total += 1
        line = line.strip()
        ls = line.split("\t")

        elem = tuple(ls[:4])
        sets = tuple(sorted(ls[4:]))
        counts[sets].append(elem)
        for a in sets:
            libs.add(a)
            for b in sets:
                simmat[(a, b)] += 1


print('total', total)
print(libs)
libs = sorted(libs)
base = sys.argv[2]

with open(base + '.venn', 'w') as F:
    kys = sorted([(len(x), x) for x in counts.keys()])
    for ln, s in kys:
        F.write("%s\t%s\n" % (len(counts[s]), ",".join(s[1:])))

with open(base + '.sim', 'w') as F:
    for l in libs:
        F.write("\t%s" % l)
    F.write("\n")

    for l in libs:
        F.write("%s" % l)
        for m in libs:
            F.write("\t%s" % simmat[(l, m)])
        F.write("\n")

with open(base + '.simf', 'w') as F:
    for l in libs:
        F.write("\t%s" % l)
    F.write("\n")

    for l in libs:
        F.write("%s" % l)
        for m in libs:
            F.write("\t%f" % (float(simmat[(l, m)])/total))
        F.write("\n")
