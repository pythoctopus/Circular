# -*- coding: utf-8 -*-
import sys, os
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())
    
from plotter import readfile
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp
'''
бесполезная функция: просто переписал вызов функции из numpy так,
чтобы было удобнее
возвращает случайную выборку по вон мизесу
'''
def randvm(n, mu = 197, k = 1.58):
    mu*= pi/180
    return(np.random.vonmises(mu, k, n))
    
'''
первые два аргумента - сравниваемые выборки
третий - логический. Если True - массив на входе двухмерный, каждый элемент:
[значение угла, длина вектора]
Просто для того чтобы можно было использовать ф-ю readfile из plotter
Возвращает значение критерия, p-val, длину вектора
'''
def MWW(arr1, arr2, length = True):
    sins, coss = [], []
    if length:
        sample1 = arr1[:,0]
        sample2 = arr2[:,0]
    else:
        sample1 = arr1
        sample2 = arr2
    n1 = len(sample1)
    n2 = len(sample2)
    total = np.sort(np.concatenate((sample1, sample2)))
    for i in range(len(total)):
        if total[i] in sample1:
            coss.append(cos(i*2*pi/(n1+n2)))
            sins.append(sin(i*2*pi/(n1+n2)))
    R2 = sum(sins)**2 + sum(coss)**2
    W = 2*(n1 + n2 - 1)*R2/(n1*n2)
    p = exp(-W/2)
    return(W, p, R2)

'''
несмотря на название, считает не мощность, а 1-мощность.
описание в доковском файле
'''    
def power_MWW_mu(effect, mu, pars2, sbNMF, N, K, pars1 = None):
    result = []
    for i in effect:
        count = 0
        for p in range(N):
            sample1 = sbNMF#randvm(pars1[0], mu, pars1[1])
            sample2 = randvm(pars2[0], mu-i, pars2[1]) #+ np.pi
            pr = MWW(sample1, sample2, False)[1]
            if pr > K:
                count += 1
        result.append(np.array([i, count/N]))
    return np.array(result)

'''
аналогично с концентрацией
'''
def power_MWW_k(effect, k, pars1, pars2, sbNMF, N, K):
    result = []
    for i in effect:
        count = 0
        for p in range(N):
            sample1 = sbNMF#randvm(pars1[0], pars1[1], k)
            sample2 = randvm(pars2[0], pars2[1], k + i)
            pr = MWW(sample1, sample2, False)[1]
            if pr > K:
                count += 1
        result.append(np.array([i, count/N]))
    return np.array(result)

if __name__ == '__main__':
    sbNMF, ci1 = readfile('SboriNMF.txt')
    sbCMF, ci2 = readfile('SboriCMF.txt')
    sbNMF, ci3 = readfile('ErubeNMF.txt')
    sbNMF = list(sbNMF[:, 0] * np.pi/180)
    for i in range(len(sbNMF)):
        if sbNMF[i] > pi:
            sbNMF[i] -= 2*pi

    sbNMF = np.array(sbNMF)
    N = 10000
    ks = np.linspace(-0.9, 6, 25)
    effect = [i for i in range(0, 151, 10)]
    powersMu = power_MWW_mu(effect, 5.147*180/pi,
                            (17, 1.744), (14, 1.292), sbNMF, N, 0.5)
    powersK = power_MWW_k(ks, 1.744, (17, 5.147*180/pi),
                          (14, 4.643*180/pi), sbNMF, N, 0.5)
 #   powersMu = power_MWW_mu(effect, 197, (25, 0.93), (20, 1.58), N)
 #   powersK = power_MWW_k(ks, 0.93, (25, 197), (20, 197), N)
    f, axes = plt.subplots(2, figsize = (8, 13))
    ax, ax2 = axes
    ax.plot(powersMu[:, 0], powersMu[:, 1], color = 'r', marker = 'o')
    ax.set_yticks(np.linspace(0, 1, 11))
    ax.set_ylabel('probability level')
    ax.set_xlabel('Difference between zero value (295°) and mean direction of theoretical distribution, degrees',
                  ha = 'left', x = 0)
    ax2.plot(powersK[:, 0], powersK[:, 1], color = 'r', marker = 'o')
    ax2.set_yticks(np.linspace(0, 1, 11))
    ax2.set_ylabel('probability level')
    ax2.set_xlabel('Difference between zero value (1.74) and concentration of theoretical distribution',
                   ha = 'left', x = 0) 
#    ax.set_title('E. rubecula\n')
    ax.set_title('E. rubecula\n')
    plt.savefig('powerERMWW2.png', dpi = 500, bbox_inches = 'tight') 
    plt.show()      
