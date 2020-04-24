# This file reads a file in fasta format and extracts the sequences
import sys
from os.path import basename, splitext

ifile = open(sys.argv[1], 'r')
outfile = open(splitext(basename(sys.argv[1]))[0]+".backtxt", 'w')
count = 1
for l in ifile.readlines():
    if l[0] == '>':
        continue
    print(l.strip(), file=outfile)
    count = count + 1
ifile.close()
outfile.close()