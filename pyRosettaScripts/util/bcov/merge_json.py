#!/usr/bin/python


import sys
import os
import gzip



def gzopen(name, mode="r"):
    if (name.endswith(".gz")):
        return gzip.open(name, mode)
    else:
        return open(name, mode)


f = gzopen(sys.argv[1], "w")
f.write("[\n")



first = True

for file in sys.stdin:
    file = file.strip()
    if (not first):
        f.write(",\n")
    first = False

    g = gzopen(file)
    f.write(g.read())
    g.close()


f.write("\n]\n")
f.close()

