''' this file takes a file containing list of epitopes and converts it to fasta
eg: convert epitopesnew.csv to .fasta

QAGTSLSLFGDS	0
LVGDRVRVASKS	1

to

>e1-0
QAGTSLSLFGDS

>e2-1
LVGDRVRVASKS
'''

import sys
from os.path import basename, splitext

ifile = open(sys.argv[1], 'r')
outfile = open(splitext(basename(sys.argv[1]))[0]+".fasta", 'w')
count = 1
for l in ifile.readlines():
    line = l.strip().split()
    seq = line[0].strip()
    # ep = line[1].strip()
    print(">e-" + str(count), file=outfile)
    print(seq, file=outfile)
    count = count + 1
ifile.close()
outfile.close()
    


