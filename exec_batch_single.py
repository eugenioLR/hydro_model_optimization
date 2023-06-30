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
    def __init__(self, rscript_name='exec_optim.R', metric="MSE", model_used=0):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)

        inf_lim = np.array([1e-4,10,0,1e-4,1e-4,2,-6])
        sup_lim = np.array([1,2000,1,100,1,11,1])
        self.size = 7

        self.metric = metric
        self.model_used = model_used

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"
        else:
            opt = "min"

        # Loading the function we have defined in R.
        self.exec_function_r = robjects.globalenv['hydro_prob']
        super().__init__(self.size, opt, sup_lim, inf_lim)

    def objetive(self, solution):
        metrics = self.exec_function_r(self.model_used, solution)
        if self.metric == "MSE":
            return float(metrics[0])
        elif self.metric == "RMSE":
            return float(metrics[1])
        elif self.metric == "NSE":
            return float(metrics[3])
        elif self.metric == "R2":
            return float(metrics[4])
        elif self.metric == "KGE":
            return float(metrics[5])
        

    def random_solution(self):
        return (self.sup_lim-self.inf_lim) * np.random.random(self.size) + self.inf_lim
    
    def check_bounds(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)

substrates_real = [
    SubstrateReal("Gauss", {"F": 0.1}),
    SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.9}),
    SubstrateReal("BLXalpha", {"Cr": 0.35}),
    SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
    SubstrateReal("Perm", {"Cr": 4/7}),
]


params = {
    "popSize": 120,
    "rho": 0.6,
    "Fb": 0.98,
    "Fd": 0.15,
    "Pd": 0.99,
    "k": 3,
    "K": 10,
    "group_subs": True,

    "stop_cond": "ngen",
    "time_limit": 400.0,
    "Ngen": 100,
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

def execute_hydro_cro(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    objfunc = HydroProblem("exec_optim.R", metric, model)
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_simple_5043_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")


import multiprocessing

proc_pool = []

for i in ["MSE", "NSE", "R2", "KGE"]:
    for j in [0,1,2,3]:
        proc_pool.append(multiprocessing.Process(target=execute_hydro_cro, args=(i, j)))

proc_pool.append(multiprocessing.Process(target=execute_hydro_cro, args=("MSE", 0)))

[p.start() for p in proc_pool]

[p.join() for p in proc_pool]