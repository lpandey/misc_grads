# -*- coding: utf-8 -*-
"""
Non-Mass-Weighted centroid from xyz coords (average coordinate)

Written by Laxman Pandey on 01.26.2011
"""
import os

def myCenter(natm,sumx,sumy,sumz):
    xc = (1.0/natm)*sumx
    yc = (1.0/natm)*sumy
    zc = (1.0/natm)*sumz
    return xc, yc, zc

def readxyz(xyzf):
    rf = open(xyzf, 'r')
    global topnatm
    topnatm = None
    natm = 0
    xl, yl, zl = [], [], []
    for line in rf:
        line = line.strip()
        if not line: continue
        line = line.split()

        # testing and getting to just the coordinate lines
        if (len(line) == 1  and line[0].isdigit()): # if num of atoms available
           topnatm = int(line[0])
           #print topnatm
        if len(line)  < 4: continue
        elm, x, y, z = line[0], line[1], line[2], line[3]
        if type(elm) != str or len(elm) > 2: continue

        try:
            xi, yi, zi = float(x), float(y), float(z)
            xl.append(xi), yl.append(yi), zl.append(zi)
            natm += 1
        except ValueError: continue # skip non-coord lines
    rf.close()    

    #print('no. of atoms: {}, collected: {}'.format(topnatm, len(xl)))

    if topnatm:
        assert (len(xl) == natm and len(xl) == topnatm), 'atom numbers issues !!'
    else:
        assert len(xl) == natm, 'atom numbers issues !!'
    sumx, sumy, sumz = sum(xl), sum(yl), sum(zl)
    return natm, sumx, sumy, sumz

# specify file to read
parentPath = os.getcwd()
xyzf = os.path.join(parentPath, 'TobyNelson', 'dP-BDT-tO-2.xyz')

# Call functions now
natm, sumx, sumy, sumz = readxyz(xyzf) 
xc, yc, zc = myCenter(natm,sumx,sumy,sumz)
print('centriod: {:15.8f} {:15.8f} {:15.8f}'.format(xc, yc, zc))
