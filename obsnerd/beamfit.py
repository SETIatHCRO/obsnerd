import numpy
from scipy.optimize import curve_fit


def func_gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))


def gaussian(data, A, mu, sigma):
# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [A, mu, sigma]
    x = list(range(len(data)))

    coeff, var_matrix = curve_fit(func_gauss, x, data, p0=p0)

    # Get the fitted curve
    fitg = func_gauss(x, *coeff)

    return coeff, fitg

