# Circular
## Description
This is a simple set of functions for estimation of II-type error probability 
for Mardia-Watsom-Wheeler test.
The basic idea is multimle comparisons of zero
sample with random samples from theoretical Von-Mises population with parameters of 
mean and concentration, similar to estimated from zero sample, but shifted by specific
value. This scrypt estimates the probability of II-type error depend on shift size.

## Usage
Some comments are written directky in the code. The working class is BetaEstimator.
It's __init__ method takes 8 arguments:
* arr1 - first numpy array of data (np.ndarray)
* arr1 - second array (another sample)
* name1 & name2 - description of arr1 and arr2 (string)
* effect_mu - the biases for mean (in radians) (list-like)
* effect_k - for concentration
* is_degrees - is data in degrees or radians (by default True)
* is3d - do u need a 3d or 2d plot of estimation (by dafault True)

Class attributes:
* sample1 & sample2 - CircSample objects, our samples arr1 & arr2
* effect_mu
* effect_k
* is3d
* MWW_p - p-value of Mardia-Watson-Wheeler test, is None untill __estimate__ method is called
* result - the result of estimation. 
> if __is3d__ is True: result is a np.ndarray of size KxM, K: effect_k size, M: effect_mu size.
> in another case result is a list of two elements, each of which is a np.ndarry of this kind: [effects(shifts), estimations].

Methods:
* estimate() - performs whe estimation. Params: n_estim - a number of generations of random samples, savefile - True or False,
only used if is3d is True. If savefile is True - saves result in txt file.
* show() - plot results

And a few worlds about sample1, sample2 and CircSample class. This class has 5 attributes and 3 methods:
* data - array of sample
* name - description
* mean - mean angle, radians (default None)
* mean_vector - tuple (mean angle, mean vector length) (default None)
* K - concentration (default None)
* rayleigh_p - p-value for rayleigh test (default None)
* compute_params() - compute mean, mean_vector, K
* rayleigh() - compute rayleigh_p
* show() - plot this sample (inner circle is a 0.05 significance level for Rayleigh test)


