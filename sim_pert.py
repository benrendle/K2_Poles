''' Perturbation Script for Trilegal using observations '''

import numpy as np
import pandas as pd

class PERTURBATION():
    def __init__(self,_sim,_err):
        self.alt_sim = _sim
        self.err = _err

    def __call__(self):
        self.pert()
        return self.alt_sim

    def pert(self):
        ''' Propagate out to functions for each of the parameters to be perturbed '''
        self.alt_sim['Mass'] = self.gauss('Mass',2)
        # self.alt_sim['radius'] = self.gauss('radius')
        # self.alt_sim['age'] = self.gauss('age')
        self.alt_sim['logAge'] = self.gauss('logAge',4)
        # self.alt_sim['M_H'] = self.gauss('M_H')
        # self.alt_sim['Teff'] = self.gauss('Teff')

    def gauss(self,param,j):
        ''' Gaussian perturber for specified parameter '''
        x = self.alt_sim[param]
        sigma = 0.0
        sig = []
        for i in self.err:
            if i[0] == param:
                sigma = i[1] # percentage uncert
        # print(sigma)
        sig = x*(sigma/j/100)
        x = np.random.normal(x,sig)
        return x
