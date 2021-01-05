import numpy as np
import matplotlib.pyplot as plt
from numba import njit, int32, float64

@njit(float64(float64, int32))
def log_poisson(mu, n0):
    sum_ = sum(np.log(n0+1))
    return n0 * np.log(mu) - mu - sum_

