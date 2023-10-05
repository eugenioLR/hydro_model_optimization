import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

from PyCROSL import SubstrateReal, AbsObjectiveFunc, CRO_SL
import random
from objective_functions import *

import multiprocessing
from itertools import product

rscript_name = "exec_optim_semidist.R"
#data_file = "data/CHGdataSIMPA.txt"
#basin_file = "data/CHGbasins.txt"
data_file = "data/CHTdataSIMPAcal.txt"
basin_file = "data/basinsSimpa.txt"

r = robjects.r
r['source'](rscript_name)
get_basin_q = robjects.globalenv['get_basin_q']

basins = pd.read_csv(basin_file)

substrates_real = [
    SubstrateReal("Gauss", {"F": 0.0001}),
    SubstrateReal("Cauchy", {"F": 0.01}),
    SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
    SubstrateReal("BLXalpha", {"F": 0.35}),
    SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
    # SubstrateReal("Perm", {"N": 2}),
]

params = {
    # "popSize": 4,
    "popSize": 100,
    "rho": 0.6,
    "Fb": 0.98,
    "Fd": 0.15,
    "Pd": 0.99,
    "k": 3,
    "K": 10,
    "group_subs": True,
    # "group_subs": False,

    "stop_cond": "Ngen",
    "time_limit": 400.0,
    "Ngen": 100,
    # "Ngen": 1,
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

    agg_q = {}
    basin_params = {}


    prev_q = 0
    for idx in basins.index:
        basin_code = basins["code"][idx]
        codedown = basins["codedown"][idx]

        prev_q = 0
        if basin_code in agg_q:
            prev_q = agg_q[basin_code]        

        objfunc = HydroSemidistModelGOF(rscript_name, data_file, basin_file, metric, model, basin_code, prev_q)
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

    # np.savetxt(f"config_semidist_5043_{model}_{metric}.csv", result, delimiter=",")
    np.savetxt(f"test.csv", result, delimiter=",")

    print(f"\n\n\nFINISHED for {metric} using model {model}")


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
    args = [('KGE', 0)]
    main(args)
