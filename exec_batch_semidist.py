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

import multiprocessing
from itertools import product

class HydroProblemSemidist(AbsObjetiveFunc):
    def __init__(self, rscript_name='exec_optim_semidist.R', metric="MSE", model_used=0, basin_code=-1, prev_q=0):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)

        inf_lim = np.array([1e-4,10,0,1e-4,1e-4,2,-6])
        sup_lim = np.array([1,2000,1,100,1,11,1])
        self.size = 7

        self.metric = metric
        self.model_used = model_used
        self.basin_code = basin_code
        self.prev_q = prev_q

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"

        # Loading the function we have defined in R.
        self.exec_function_r = robjects.globalenv['eval_basin_param']
        super().__init__(self.size, opt, sup_lim, inf_lim)

    def objetive(self, solution):
        metrics = self.exec_function_r(self.model_used, solution, self.basin_code, self.prev_q)
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

rscript_name = "exec_optim_semidist.R"

r = robjects.r
r['source'](rscript_name)
get_basin_q = robjects.globalenv['get_basin_q']

basins = pd.read_csv("data/CHGbasins.txt")

substrates_real = [
    SubstrateReal("Gauss", {"F": 0.0001}),
    SubstrateReal("Cauchy", {"F": 0.01}),
    SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
    SubstrateReal("BLXalpha", {"Cr": 0.35}),
    # SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
    # SubstrateReal("Perm", {"Cr": 4/7}),
]

params = {
    "popSize": 125,
    "rho": 0.6,
    "Fb": 0.98,
    "Fd": 0.2,
    "Pd": 0.99,
    "k": 3,
    "K": 3,
    "group_subs": False,

    "stop_cond": "neval",
    "time_limit": 400.0,
    "Ngen": 100,
    "Neval": 12000,
    "fit_target": 1000,

    "verbose": True,
    "v_timer": 1,

    "dynamic": True,
    "dyn_method": "fitness",
    "dyn_metric": "best",
    "dyn_steps": 500,
    "prob_amp": 0.03
}

def execute_hydro_cro(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    agg_q = {}
    basin_params = {}


    prev_q = 0
    for idx in basins.index:
        basin_code = basins["code"][idx]
        codedown = basins["codedown"][idx]

        prev_q = 0
        if basin_code in agg_q:
            prev_q = agg_q[basin_code]        

        objfunc = HydroProblemSemidist(rscript_name, metric, model, basin_code, prev_q)
        c = CRO_SL(objfunc, substrates_real, params)
        
        print(f"Optimizing basin with code {basin_code}.\n\n")

        c.safe_optimize()

        output_name = f"config_semidist_nd_{model}_{metric}_{basin_code}"

        c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
        # c.save_solution(output_name+".csv")

        basin_params[basin_code] = c.best_solution()[0]

        if codedown not in agg_q:
            agg_q[codedown] = 0

        agg_q[codedown] += get_basin_q(model, basin_params[basin_code], basin_code, prev_q)
    
    
    result = []
    for i in basins["code"]:
        result += [basin_params[i]]
    result = np.array(result)

    np.savetxt(f"config_semidist_5043_{model}_{metric}.csv", result, delimiter=",")

    print(f"\n\n\nFINISHED for {metric} using model {model}")


def execute_hydro_cro_wrapper(x):
    execute_hydro_cro(*x)

def main():
    pool = multiprocessing.Pool(processes=16)

    args = product(["MSE", "NSE", "R2", "KGE"], [0,1,2,3])

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
