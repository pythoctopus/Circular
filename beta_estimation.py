import numpy as np
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
#from tqdm import tqdm # if running not in Jupyter

class BetaEstimator:
    def __init__(self, arr1, arr2, name1, name2, effect_mu, effect_k,
                 is_degrees=True, is3d=True):
        
        self.sample1 = CircSample(arr1, name1, is_degrees=is_degrees)
        self.sample2 = CircSample(arr2, name2, is_degrees=is_degrees)
        self.sample1.compute_params()
        self.sample2.compute_params()
        self.effect_mu = effect_mu
        self.effect_k = effect_k
        self.is3d = is3d
        self.MWW_p = None
        self.result = None
        
    def estimate(self, n_estim=1000, savefile=False):
        W, p = MWW(self.sample1.data, self.sample2.data, isreal=True)
        self.MWW_p = p
        if self.is3d:
            self.result = beta3d(self.effect_mu, self.effect_k, self.sample1.data, 
                                 self.sample1.data.shape[0], self.sample1.K,
                                 self.sample1.mean, critical_p=p,
                                 n_estim=n_estim, savefile=savefile,
                                 path='{}_{}_MWW_beta3d.txt'.format(self.sample1.name,
                                                                   self.sample2.name))
            
        else:
            mus, kas = beta2d(self.effect_mu, self.effect_k, self.sample1.data, 
                              self.sample1.data.shape[0], self.sample1.K, self.sample1.mean,
                              n_estim=n_estim, critical_p=p)
            self.result = [mus, kas]
    
    def show(self, savefile=False):
        if self.is3d:
            f, ax = plt.subplots(figsize=(10,8))
            ax.set_xlabel('effect on mean')
            ax.set_ylabel('effect on concentration')
            color = np.zeros((11,3))
            color[:, 0] = np.linspace(0, 1, 11)
            cs = ax.contour(self.effect_mu * 180/np.pi, self.effect_k, self.result,
                           levels = [0.2, 0.4, 0.6, 0.8], colors='#969696')
            css = ax.contourf(self.effect_mu * 180/np.pi, self.effect_k, self.result,
                             levels = np.linspace(0, 1, 11), colors=color)
            ax.clabel(cs)
            f.colorbar(css, label='probability level')
            f.suptitle('{0}\nMWW p-val: {1:.3f}'.format(self.sample1.name, self.MWW_p))
            if savefile:
                plt.savefig('{}.png'.format(self.sample1.name), dpi=500)
            plt.show()
        else:
            f, ax = plt.subplots(2, figsize=(8,13))
            ax[0].plot(self.result[0][:, 0] * 180/np.pi, self.result[0][:, 1], color='r', marker='o')
            ax[1].plot(self.result[1][:, 0], self.result[1][:, 1], color='r', marker='o')
            for i in 0,1:
                ax[i].set_yticks(np.linspace(0, 1, 11))
                ax[i].set_ylabel('probability level')
                ax[i].set_xlabel('effect size')
                ax[i].grid()
                
            ax[0].set_xticks(self.effect_mu * 180/np.pi)
            ax[0].tick_params(axis='x', rotation=45)
            ax[0].set_title('Mean')
            ax[1].set_title('Concentration')
            f.suptitle('{0}\nMWW p-val: {1:.3f}'.format(self.sample1.name, self.MWW_p))
            if savefile:
                plt.savefig('{}.png'.format(self.sample1.name), dpi=500, bbox_inches='tight')
            plt.show()
                
class CircSample:
    def __init__(self, arr, name, is_degrees=True):
        if is_degrees:
            self.data = arr * np.pi/180
        else:
            self.data = arr
            
        self.name = name
        self.mean = None
        self.mean_vector = None
        self.K = None
        self.rayleigh_p = None
    
    def compute_params(self):
        A, l = compute_mean(self.data)
        self.mean = A
        self.mean_vector = (A, l)
        self.K = concentration(l, self.data.shape[0])
        
    def rayleigh(self):
        n = self.data.shape[0]
        r = self.mean_vector[1]
        L = 0.5 * np.exp(-n * r**2)
        self.rayleigh_p = L
    
    def show(self):
        x = np.log(0.05)
        n = self.data.shape[0]
        critical_vlen = ((- x - (2*x + x**2)/(4*n))/n)**0.5
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes(projection='polar')
        ax.set_rorigin(0)
        ax.set_rlim(0, 1.06)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rgrids((critical_vlen, 0), labels=None)
        ax.set_yticklabels([])
        ax.set_thetagrids((0, 90, 180, 270), ('0', '90', '180', '270'))
        ax.plot(self.data, np.ones(self.data.shape), 'o', color='black',
               markersize=12)
        ax.quiver(0, 0, self.mean, self.mean_vector[1], scale_units='xy',
                  scale=1, angles='xy', zorder=3)
        fig.suptitle('{0}\nmd, deg: {1:.3f}, concentration: {2:.3f}\nRayleigh p: {3:.3f}, N: {4}'.\
                     format(self.name, self.mean*180/np.pi, self.K, self.rayleigh_p, self.data.shape[0]))
        plt.show()

####################################################################################
# beta estimation
def beta2d(effect_mu, effect_k, sample, sN, sK, sMu, n_estim=1000, critical_p=0.05):
    fin_k, fin_mu = np.zeros(2), np.zeros(2)
    for i in tqdm(effect_mu):
        count = 0
        for p in range(n_estim):
            sample1 = np.random.vonmises(sMu-i, sK, sN)
            sample1 += (sample1 < 0) * (2*np.pi)
            pr = MWW(sample, sample1)[1]
            if pr >= critical_p:
                count += 1
        fin_mu = np.vstack((fin_mu, np.array([i, count/n_estim])))
        
    for i in tqdm(effect_k):
        count = 0
        for p in range(n_estim):
            sample1 = np.random.vonmises(sMu, sK+i, sN)
            sample1 += (sample1 < 0) * (2*np.pi)
            pr = MWW(sample, sample1)[1]
            if pr >= critical_p:
                count += 1
        fin_k = np.vstack((fin_k, np.array([i, count/n_estim])))
        
    return fin_mu[1:], fin_k[1:]

def beta3d(effect_mu, effect_k, sample, sN, sK, sMu, n_estim=1000,
           critical_p=0.05, savefile=True, path='beta3d_result.txt'):
    
    result = np.zeros(effect_mu.shape[0])
    for i in tqdm(effect_k):
        this_k = []
        for j in effect_mu:
            count = 0
            for p in range(n_estim):
                sample1 = np.random.vonmises(sMu-j, sK+i, sN)
                sample1 += (sample1 < 0) * (2*np.pi)
                pr = MWW(sample, sample1)[1]
                if pr >= critical_p:
                    count += 1
            this_k.append(count/n_estim)
        result = np.vstack((result, np.array(this_k)))
    
    if savefile:
        np.savetxt(path, result[1:], delimiter=',')
    return result[1:]

def MWW(arr1, arr2, isreal=False):
    if isreal:
        arr1 += np.random.rand(arr1.shape[0])*0.01 - 0.005
    S1, C1, S2, C2 = 0, 0, 0, 0
    n1 = arr1.shape[0]
    n2 = arr2.shape[0]
    total = np.sort(np.concatenate((arr1, arr2)))
    for i in range(n1+n2):
        if total[i] in arr1:
            S1 += np.sin(i*2*np.pi/(n1+n2))
            C1 += np.cos(i*2*np.pi/(n1+n2))
        else:
            S2 += np.sin(i*2*np.pi/(n1+n2))
            C2 += np.cos(i*2*np.pi/(n1+n2))
    W = (S1**2 + C1**2)/n1 + (S2**2 + C2**2)/n2
    p = np.exp(-W)
    return W*2, p

###################################################################################
# sample properties
def compute_mean(arr):
    N = np.shape(arr)[0]
    S = np.sum(np.sin(arr))
    C = np.sum(np.cos(arr))
    tg = S/C
    l = (S**2 + C**2)**0.5 /N
    A = np.arctan(tg)
    if C < 0:
        A = A + np.pi
    elif S < 0:
        A = 2*np.pi + A
    return A, l

def concentration(r, n):
    if r < 0.53:
        c = 2*r + r**3 + 5*r**5/6
    elif r >= 0.85:
        c = 1/(3*r - 4*r**2 + r**3)
    else:
        c = -0.4 + 1.39*r + 0.43/(1 - r)
        
    if n < 15:
        if c < 2:
            c = max([c - 2/(n*c), 0])
        else:
            c = (n -1)**3 * c / (n * (n**2 + 1))
            
    return(c)
    
