import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

from PyCRO_SL import *
from PyCRO_SL.CRO_SL import CRO_SL
from PyCRO_SL.AbsObjetiveFunc import AbsObjetiveFunc
from PyCRO_SL.SubstrateInt import SubstrateInt
from PyCRO_SL.SubstrateReal import SubstrateReal
import random

class HydroProblem(AbsObjetiveFunc):
    def __init__(self, rscript_name='exec_problem.R', opt="min"):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)

        inf_lim = np.array([1e-4,10,0,1e-4,1e-4,2,-6])
        sup_lim = np.array([1,2000,1,100,1,11,1])
        self.size = 7

        # Loading the function we have defined in R.
        self.exec_function_r = robjects.globalenv['hydro_prob']
        super().__init__(self.size, opt, sup_lim, inf_lim)

    def objetive(self, solution):
        metrics = self.exec_function_r(0, "NSE", solution)
        return metrics[0]

    def random_solution(self):
        return (self.sup_lim-self.inf_lim) * np.random.random(self.size) + self.inf_lim
    
    def check_bounds(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)


objfunc = HydroProblem("exec_problem.R")


substrates_real = [
    SubstrateReal("Gauss", {"F": 12.16}),
    SubstrateReal("Cauchy", {"F": 0.02}),
    SubstrateReal("2point"),
    SubstrateReal("DE/best/2", {"F":1.4, "Cr":0.9}),
    SubstrateReal("BLXalpha", {"Cr": 0.35}),
    SubstrateReal("Firefly", {"a":0.7, "b":1, "d":0.95, "g":10})
]

params = {
    "popSize": 100,
    "rho": 0.8,
    "Fb": 0.98,
    "Fd": 0.15,
    "Pd": 0.99,
    "k": 3,
    "K": 10,
    "group_subs": True,

    "stop_cond": "ngen",
    "time_limit": 400.0,
    "Ngen": 2000,
    "Neval": 3e4,
    "fit_target": 1000,

    "verbose": True,
    "v_timer": 1,

    "dynamic": True,
    "dyn_method": "fitness",
    "dyn_metric": "best",
    "dyn_steps": 500,
    "prob_amp": 0.015
}
c = CRO_SL(objfunc, substrates_real, params)

c.safe_optimize()