#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


class hist():
    def __init__(self, cnp_hist):
        """
        Init with np.hist, which is 2 arrays,
        first one is weihgt,
        second one is bin edges.
        """
        self.weight = cnp_hist[0]
        self.bin_edge = cnp_hist[1]
    
    def show(self):
        plt.grid()
        plt.hist(self.bin_edge[:-1], self.bin_edge, weights=self.weight);
        
    def get_bin_middle(self):
        """
        return array-like object, mid-point of each bin.
        """
        return self.bin_edge[1:] + np.diff(self.bin_edge)/2
        
    def selective_sampling(self, fsize):
        """
        Give a histogram, sampling it.
        """
    
        fbin_middle = self.get_bin_middle()
    
        fweight = self.weight
    
        sampled_list = np.array([])
    
        while sampled_list.size < fsize:
            randy = np.random.uniform(0, fweight.max(), size = fsize)
            randx = np.random.uniform(fbin_middle[0], fbin_middle[-1], size = fsize)
            yy = np.interp(randx, fbin_middle, fweight)
            sampled_list = np.hstack((sampled_list, randx[randy<yy]))
        if sampled_list.size > fsize:
            return sampled_list[:fsize]
        else:
            return sampled_list
        
    def total(self):
        if self.weight.dtype == 'float64':
            print("Histogram has been normalized!!")
        else:
            return np.sum(self.weight)
        
    def mean(self):
        """
        有偏大的嫌疑。需要调整
        """
        fbin_middle = self.get_bin_middle()
        return np.sum(fbin_middle*self.weight)/self.total()


# In[3]:


def generate_distribution(mu, fsize):
    """
    randomly generate two gaussian distribution.
    """
    rng = np.random.RandomState(10)  # deterministic random data
    rand_dis = np.hstack((np.random.uniform(size=1)*rng.normal(loc=-mu, scale=1, size=fsize),
                   np.random.uniform(size=1)*rng.normal(loc=mu, scale=1, size=fsize)))
    return rand_dis # np.histogram(a, bins="auto", density=True)


# In[58]:


if __name__ == "__main__":
    h = hist(np.histogram(generate_distribution(3, 30000), bins="auto"))
    h2 = hist(np.histogram(generate_distribution(3, 30000)**2, bins="auto"))
    h3 = hist(np.histogram(h.selective_sampling(10000)))
    print(str(h3.mean()) + ' ' + str(h.mean()))

