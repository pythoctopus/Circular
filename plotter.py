# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from math import pi, sin, cos, log, atan, exp
import numpy as np
'''
Считывает текстовый файл. Если в файле указаны значения доверительного
интервала, параметру CI нужно присвоить значение True(стоит по умолчанию),
иначе - False
Возвращает два объекта
'''
def readfile(path, CI = True):
    with open(path, 'r') as inf:
        arr = inf.read().strip().split('\n')
        arr = [float(i) for i in arr]
    if CI == True:
        data = arr[:-2]
        CI95 = arr[-2:]
    else:
        data = arr
        CI95 = None
    for i in range(len(data)):
        data[i] = np.array([data[i], 1])
    return(np.array(data), CI95)

'''
Создает основу графика - объект Axes, с некоторыми надстройками, 
принимает два логических значения. Если на графике нужны доверительные
интервалы и внутренняя окружность с критическим значением - оба True.
Вспомогательная функция.
'''
def cor(CI = False, p = False):
    ax = plt.axes(projection = 'polar')
    ax.set_title('gN\n')
    ax.set_rorigin(0)
    ax.set_rlim(0, 1.06)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if not CI:
        ax.set_thetagrids((0, 90, 180, 270),
                          ('0', '90', '180', '270'))
    if not p:
        ax.set_rgrids((0, 0),labels = ('', ''))
    return(ax)
    
'''
Считает критическое значение для длины среднего вектора
'''
def crcal(a, n): #compute critical value of resultant r for p-lvel a with sample n
    x = log(a)
    K = ((- x - (2*x + x**2)/(4*n))/n)**0.5
    return(K)

'''
Рассчет среднего вектора и направления - нету в питоне модуля для
круговой статистики. Возвращает длину вектора и угол в радианах.
На вход принимает массив формата выходного из ф-ии readfile
'''    
def md(arr, l = True):
    S, C, n = 0, 0, len(arr)
    if l:
        for i in arr:
            S += sin(i[0]*pi/180) * i[1]
            C += cos(i[0]*pi/180) * i[1]
    else:
        for i in arr:
            S += sin(i)
            C += cos(i)
    R = (S**2 + C**2)**0.5
#    print(S, C, R)
    r = R/n
    A = atan(S/C)
    if C < 0:
        A = A + pi
    elif S < 0:
        A = 2*pi + A
    return(r, A)

'''
Считает концентрацию
'''    
def conc(r):
    if r < 0.53:
        c = 2*r + r**3 + 5*r**5/6
    elif r >= 0.85:
        c = 1/(3*r - 4*r**2 + r**3)
    else:
        c = -0.4 + 1.39*r + 0.43/(1 - r)
    return(c)

'''
критерий релея по Мардии 1972
'''    
def r_test(r, n):
    L = 0.5*exp(-n*r**2)
    return(L)
    
'''
Строит график на основе одной выборки
data - массив формата выходного из readfile: [[fi, 1]....[fiN, 1]]
CI, p - наличие на графике
ci - array-like, 2 значения: среднее направление +- CI
pval - 0.05 or 0.01
color - цвет, см color в matplotlib.pyplot
(examples 2-3)
'''
def single(data, CI = True, p = True, ci = None, pval = None, color = 'red',
           markersize = 7):
    r, fi = md(data)
    N = len(data)
    ax = cor(CI=CI, p=p)
    ax.plot(data[:,0]*pi/180, data[:,1], 'o', color=color,
            markersize = markersize)
    ax.text(0, -0.1, 'N = %d' %N, fontsize = 12, transform=ax.transAxes)
    ax.text(0.8, -0.1, 'r = %.2f' %r, fontsize = 12, transform=ax.transAxes)
    ax.quiver(0,0, sin(fi), cos(fi), scale = 2*1.06/r, zorder = 3)
    if p:
        K = crcal(pval, N)
        chapter = fi*2//pi
        angle = (chapter-1)*pi/2 + pi/4
        ax.set_rgrids((K, 0), labels = ('p<%.2f' %pval, ''),
                       angle = angle*180/pi)
    if CI:
        et = np.array([0, 1/2, 1, 3/2, 2])
        delta1 = np.fabs(et - ci[0]/180)
        delta2 = np.fabs(et - ci[0]/180)
        critical = 15/180
        labs = ['0', '90', '180', '270']
        for i in range(5):
            if delta1[i] <= critical or delta2[i] <= critical:
                if i == 4:
                    labs[0] = ''
                else:
                    labs[i] = ''
        labs.append('-CI95')
        labs.append('+CI95')
        ax.set_thetagrids((0, 90, 180, 270, ci[0], ci[1]),
                           tuple(labs))
    return(ax)

'''
arrs - array-like, несколько выборок
colors - array-like, цвет для каждой выборки
(example1)
'''
def multiple(arrs, colors):
    ax = cor()
    for i in range(len(arrs)):
        r, fi = md(arrs[i])
        ax.quiver(0,0, sin(fi) ,cos(fi), scale = 2*1.06/r, color=colors[i])
        ax.plot(arrs[i][:,0]*pi/180, arrs[i][:,1], 'o',
                color=colors[i], markersize = 7)
    return(ax)
        
'''
Параметры значат то же, что и в 2-х предыдущих функциях.
Это на случай, если нужно объединить две выборки в одну, но отобразить
их на одном графике как разные. На примере понятнее: я сделал так в
дипломе, объединив контроль и опыт у славок, и сделав график, на котором
птицы из контроля и опыта обозначены разными цветами, но один
средний вектор и доверительный интервал на всю выборку
(example4)
'''
def combine(arrs, colors, CI = True, p = True, ci = None, pval = None,
            markersize = 7, text = False):
    ax = cor(CI, p)
    data = np.concatenate(arrs, 0)
    r, fi = md(data)
    N = len(data)
    ax.quiver(0,0, sin(fi), cos(fi), scale = 2*1.06/r, zorder = 3)
    if text:
        ax.text(0, -0.1, 'N = %d' %len(arrs[0]), fontsize = 12,
                transform=ax.transAxes)
        ax.text(0.8, -0.1, 'r = %.2f' %r, fontsize = 12, transform=ax.transAxes)
    for i in range(len(arrs)):
        r, fi = md(arrs[i])
        ax.plot(arrs[i][:,0]*pi/180, arrs[i][:,1], 'o',
                color=colors[i], markersize = markersize)
    if p:
        K = crcal(pval, N)
        chapter = fi*2//pi
        angle = (chapter-1)*pi/2 + pi/4
        ax.set_rgrids((K, 0), labels = ('p<%.2f' %pval, ''),
                       angle = angle*180/pi)
    if CI:
        et = np.array([0, 1/2, 1, 3/2, 2])
        delta1 = np.fabs(et - ci[0]/180)
        delta2 = np.fabs(et - ci[0]/180)
        critical = 15/180
        labs = ['0', '90', '180', '270']
        for i in range(4):
            if delta1[i] <= critical or delta2[i] <= critical:
                if i == 4:
                    labs[0] = ''
                else:
                    labs[i] = ''
        labs.append('-CI95')
        labs.append('+CI95')
        ax.set_thetagrids((0, 90, 180, 270, ci[0], ci[1]),
                           tuple(labs))
    return(ax)
        
if __name__ == '__main__': #тут несколько примеров использования
    sbNMF, ci1 = readfile('SboriNMF.txt')
    sbCMF, ci2 = readfile('SboriCMF.txt')
    erNMF, n = readfile('ErubeNMF.txt')
    erCMF, n = readfile('ErubeCMF.txt', False)
#    multiple((sbNMF,sbCMF), ('red', 'blue'))
#    plt.savefig('example1.png', format = 'png', bbox_inches='tight', 
#                dpi = 500)
#    plt.show()
#    single(erNMF, CI= False, pval = 0.01)
#    plt.savefig('example2.png', format = 'png', bbox_inches='tight', 
#                dpi = 500)
#    plt.show()
    print(conc(md(erNMF)[0]), md(erNMF))
    print(conc(md(erCMF)[0]), md(erCMF))
    single(erNMF, pval = 0.01, ci = ci1)
    plt.savefig('example3.png', format = 'png', bbox_inches='tight', 
                dpi = 500)
    plt.show()
#    combine((sbNMF, sbCMF), ('red', 'blue'),
#                  ci = (176.092, 219.299), pval = 0.01)
#    plt.savefig('example4.png', format = 'png', bbox_inches='tight', 
#                dpi = 500)
#    plt.show()
