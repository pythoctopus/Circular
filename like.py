# -*- coding: utf-8 -*-
import numpy as np
import sys, os
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from plotter import md
def century(ci, ciR, n, mu, k, N = 10000):
    answer = []
    lci, uci = ci
    lr, ur = ciR
    for i in k:
        answer.append([])
        #print('and one more time')
        for p in mu:
            count = 0
            for j in range(N):
                sample = np.random.vonmises(p, i, n)
                r, A = md(sample, False)
                if lci <= A <= uci and lr <= r <= ur:
                    count += 1
            L = count/N
            answer[-1].append(L)
    return np.array(answer)

if __name__ == '__main__':
    cisbNMF = np.array([162.362, 233.974]) * np.pi/180
    cisbCMF = np.array([171.909, 222.676]) * np.pi/180
    ciRNF = [0.2762572521479965, 0.6164401082174988]
    ciRCF = [0.45649609273939273, 0.799670346646901]
    n1 = 25
    n2 = 20
    k = np.linspace(0.5, 5, 100)
    mu = np.linspace(150, 300, 150) * np.pi/180
    ans1 = century(cisbNMF, ciRNF, n1, mu, k)
    np.savetxt('ans1.txt', ans1, delimiter=',')
    del ans1
    ans2 = century(cisbCMF, ciRCF, n2, mu, k)
    np.savetxt('ans2.txt', ans2, delimiter=',')
