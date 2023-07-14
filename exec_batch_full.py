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
data_file = "data/CHGdataSIMPA.txt"
basin_file = "data/CHGbasins.txt"

# substrates_real = [
#     SubstrateReal("Gauss", {"F": 0.01}),
#     SubstrateReal("Cauchy", {"F": 0.01}),
#     SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
#     SubstrateReal("BLXalpha", {"F": 0.35}),
#     SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
#     SubstrateReal("Perm", {"N": 2}),
# ]

substrates_real = [
    SubstrateReal("Gauss", {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])}),
    SubstrateReal("Cauchy", {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])}),
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
    # "Ngen": 2,
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



weights = []
basin_df = pd.read_csv(basin_file)
basin_df.sort_values(by=["order"])

full_data = pd.read_csv(data_file)
for i, idx in enumerate(basin_df.index):
    basin_code = basin_df["code"][idx]
    total_caudal = full_data[full_data["qhmobs"] != -100][full_data["code"] == basin_code]["qhmobs"].mean()
    weights.append(total_caudal)

weights = np.asarray(weights)/sum(weights)
print(weights)
    

def execute_hydro_cro(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    #objfunc = HydroProblemFull("exec_optim_semidist.R", "data/CHGbasins.txt", metric, model, basin_code=5043)


    objfunc = HydroFullModelGOF("exec_optim_semidist.R", data_file, basin_file, metric, model, weights=weights)#, basin_code=5043)
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_full_qmm_5043_{model}_{metric}"
    # output_name = f"test"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")


def execute_hydro_cro_wrapper(x):
    execute_hydro_cro(*x)

def main():
    pool = multiprocessing.Pool(processes=17)

    # args = product(["MSE", "NSE", "R2", "KGE"], [0,1,2,3])
    args = product(["MSE", "NSE", "KGE"], [0,1,2,3])

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    """
    """
    main()
