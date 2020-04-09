#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())
    
from plotter import readfile
from r_sample import MWW
import numpy as np
print('lets rock!')    
def main(dk, dm, sbNMF, K0 = 1.744, Mu0 = 5.147,
         n1 = 17, n2 = 14, N = 10): 
    result = []
    step = 0
    for i in dk:
        result.append([])
        for p in dm:
            count = 0
            for j in range(N):
                step += 1
                if step%10000 == 0:
                    print(step) 
                sample1 = sbNMF#np.random.vonmises(Mu0, K0, n1)
                sample2 = np.random.vonmises(Mu0 - p, K0 + i, n2)
                pr = MWW(sample1, sample2, False)[1]
                if pr > 0.5:
                    count += 1
            result[-1].append(count/N)
    return np.array(result)

if __name__ == '__main__':
    dks = np.linspace(-0.9, 5, 100)
    dms = np.linspace(0, 150, 100) * np.pi/180
    sbNMF, ci1 = readfile('ErubeNMF.txt')
    sbNMF = list(sbNMF[:, 0] * np.pi/180)
    for i in range(len(sbNMF)):
        if sbNMF[i] > np.pi:
            sbNMF[i] -= 2*np.pi

    sbNMF = np.array(sbNMF)        
    ans = main(dks, dms, sbNMF)
    np.savetxt('porn.txt', ans, delimiter=',')
        