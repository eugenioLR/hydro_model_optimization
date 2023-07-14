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

class HydroProblemFull(AbsObjetiveFunc):
    def __init__(self, rscript_name='exec_optim_semidist.R', basin_data="basinsABCD.txt", metric="MSE", model_used=0, weights=None, basin_code=3005):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)
        
        basin_df = pd.read_csv(basin_data)
        basin_df.sort_values(by=["order"])

        inf_lim = np.tile(np.array([1e-4,10,0,1e-4,1e-4,2,-6]), len(basin_df))
        sup_lim = np.tile(np.array([1,2000,1,100,1,11,1]), len(basin_df))
        self.size = 7*len(basin_df)

        self.metric = metric
        self.model_used = model_used
        self.basin_df = basin_df
        self.basin_code = basin_code

        if weights is None:
            weights = np.full(len(basin_df), 1/len(basin_df))
        self.weights = weights

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"

        # Loading the function we have defined in R.
        self.get_basin_q = robjects.globalenv['get_basin_q']
        self.exec_function_r = robjects.globalenv['eval_basin']
        super().__init__(self.size, opt, sup_lim, inf_lim)

    def objetive(self, solution):
        agg_q = {}
        result_q = {}
        
        result = 0

        for i, idx in enumerate(self.basin_df.index):
            basin_code = self.basin_df["code"][idx]
            codedown = self.basin_df["codedown"][idx]

            prev_q = 0
            if basin_code in agg_q:
                prev_q = agg_q[basin_code]
            
            result_q[basin_code] = self.get_basin_q(self.model_used, solution, basin_code, prev_q)
            if codedown not in agg_q:
                agg_q[codedown] = 0
            agg_q[codedown] += result_q[basin_code]

            metrics = self.exec_function_r(result_q[basin_code], basin_code)
            if self.metric == "MSE":
                result_aux = float(metrics[0])
            elif self.metric == "RMSE":
                result_aux = float(metrics[1]) 
            elif self.metric == "NSE":
                result_aux = float(metrics[3])
            elif self.metric == "R2":
                result_aux = float(metrics[4])
            elif self.metric == "KGE":
                result_aux = float(metrics[5])

            result += result_aux * self.weights[i]
        return result
        

    def random_solution(self):
        return (self.sup_lim-self.inf_lim) * np.random.random(self.size) + self.inf_lim
    
    def check_bounds(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)

substrates_real = [
    SubstrateReal("Gauss", {"F": 0.0001}),
    SubstrateReal("Cauchy", {"F": 0.01}),
    SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
    SubstrateReal("BLXalpha", {"Cr": 0.35}),
    # SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
    # SubstrateReal("Perm", {"Cr": 4/7}),
]

params = {
    "popSize": 120,
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



weights = []
basin_df = pd.read_csv("data/CHGbasins.txt")
basin_df.sort_values(by=["order"])

full_data = pd.read_csv("data/CHGdataSIMPA.txt")
for i, idx in enumerate(basin_df.index):
    basin_code = basin_df["code"][idx]
    total_caudal = full_data[full_data["qhmobs"] != -100][full_data["code"] == basin_code]["qhmobs"].mean()
    weights.append(total_caudal)

weights = np.asarray(weights)/sum(weights)
print(weights)
    

def execute_hydro_cro(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    #objfunc = HydroProblemFull("exec_optim_semidist.R", "data/CHGbasins.txt", metric, model, basin_code=5043)


    objfunc = HydroProblemFull("exec_optim_semidist.R", "data/CHGbasins.txt", metric, model, weights=weights, basin_code=5043)
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_full_qmm_5043_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")


def execute_hydro_cro_wrapper(x):
    execute_hydro_cro(*x)

def main():
    pool = multiprocessing.Pool(processes=17)

    args = product(["MSE", "NSE", "R2", "KGE"], [0,1,2,3])

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    """
    """
    main()
