import multiprocessing
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

from PyCROSL import SubstrateReal, AbsObjectiveFunc, CRO_SL
import random
from objective_functions import *


substrates_real = [
    # SubstrateReal("Gauss", {"F": 0.1}),
    # SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.9}),
    # SubstrateReal("BLXalpha", {"F": 0.35}),
    # SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
    # SubstrateReal("Perm", {"N": 2}),

    SubstrateReal("Cauchy", {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])}),
    SubstrateReal("MutNoise", {"method": "Gauss", "F": 1e-4, "N": 1}),
    SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
    SubstrateReal("BLXalpha", {"F": 0.35}),
    SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
]

params = {
    "popSize": 100,
    "rho": 0.6,
    "Fb": 0.98,
    "Fd": 0.15,
    "Pd": 0.99,
    "k": 3,
    "K": 10,
    "group_subs": True,

    "stop_cond": "Neval",
    "time_limit": 400.0,
    "Ngen": 100,
    "Neval": 3e4,
    "fit_target": 1000,

    "verbose": True,
    "v_timer": 1,

    "Njobs": 3,

    "dynamic": True,
    "dyn_method": "fitness",
    "dyn_metric": "best",
    "dyn_steps": 500,
    "prob_amp": 0.015
}

def execute_hydro_cro(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    objfunc = HydroSimpleModelGOF("exec_optim.R", "data/CHGdataSIMPA5043AG.txt", "data/CHGbasins5043AG.txt", metric, model)
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_simple_5043_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")

def execute_hydro_cro_wrapper(x):
    execute_hydro_cro(*x)

def main(args):
    pool = multiprocessing.Pool(processes=16)

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    # args = product(["MSE", "NSE", "R2", "KGE"], [0,1,2,3])
    # args = product(["MSE", "NSE", "KGE"], [0,1,2,3])

    args = [('NSE', 1), ('MSE', 1), ('KGE', 1)]
    
    main(args)
