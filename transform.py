import sys
import numpy as np
readdata = open(sys.argv[1], 'r')
aapdata = {}
for i in readdata.readlines():
    if 'X' not in i.split()[0].strip():
        aapdata[i.split()[0]] = np.log(float(i.split()[1]))
normalstats = open(sys.argv[1]+".normal", 'w')
mx = max(aapdata.values())
mn = min(aapdata.values())
print(mx,mn)
for i in aapdata.keys():
    aapdata[i] = 2 * ((aapdata[i] - mn) / (mx - mn)) - 1
    print(i, aapdata[i], file=normalstats)
readdata.close()
normalstats.close()
