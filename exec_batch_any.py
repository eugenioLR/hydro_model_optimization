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

basin_n = 5043

rscript_name_single = "exec_optim.R"
data_file_single = "data/CHGdataSIMPA5043AG.txt"
basin_file_single = "data/CHGbasins5043AG.txt"

rscript_name_cascade = "exec_optim_semidist.R"
data_file_cascade = "data/CHGdataSIMPA.txt"
basin_file_cascade = "data/CHGbasins.txt"

r = robjects.r
r['source'](rscript_name_single)
r['source'](rscript_name_cascade)
get_basin_q = robjects.globalenv['get_basin_q']

weights = []
basin_df = pd.read_csv(basin_file_cascade)
basin_df.sort_values(by=["order"])

full_data = pd.read_csv(data_file_cascade)
for i, idx in enumerate(basin_df.index):
    basin_code = basin_df["code"][idx]
    total_caudal = full_data[full_data["qhmobs"] != -100][full_data["code"] == basin_code]["qhmobs"].sum()
    weights.append(total_caudal)

weights = np.asarray(weights)/sum(weights)
print(weights)

def execute_hydro_cro_single(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    objfunc = HydroSimpleModelGOF(rscript_name_single, data_file_single, basin_file_single, metric, model)
    substrates_real = [
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

        "Njobs": 1,

        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 500,
        "prob_amp": 0.015
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_single_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")

def execute_hydro_cro_cascade(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    substrates_real = [
        SubstrateReal("Gauss", {"F": 0.0001}),
        SubstrateReal("Cauchy", {"F": 0.01}),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
        SubstrateReal("BLXalpha", {"F": 0.35}),
        SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
        # SubstrateReal("Perm", {"N": 2}),
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

        "stop_cond": "Ngen",
        "time_limit": 400.0,
        "Ngen": 100,
        "Neval": 4e4,
        "fit_target": 1000,

        "verbose": True,
        "v_timer": 1,

        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 500,
        "prob_amp": 0.015
    }

    agg_q = {}
    basin_params = {}


    prev_q = 0
    for idx in basin_df.index:
        basin_code = basin_df["code"][idx]
        codedown = basin_df["codedown"][idx]

        prev_q = 0
        if basin_code in agg_q:
            prev_q = agg_q[basin_code]        

        objfunc = HydroSemidistModelGOF(rscript_name_cascade, data_file_cascade, basin_file_cascade, metric, model, basin_code, prev_q)
        c = CRO_SL(objfunc, substrates_real, params)
        
        print(f"Optimizing basin with code {basin_code}.\n\n")

        c.safe_optimize()

        output_name = f"config_cascade_nd_{model}_{metric}_{basin_code}"

        c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")

        basin_params[basin_code] = c.best_solution()[0]

        if codedown not in agg_q:
            agg_q[codedown] = 0

        agg_q[codedown] += get_basin_q(model, basin_params[basin_code], basin_code, prev_q)
    
    
    result = []
    for i in basin_df["code"]:
        result += [basin_params[i]]
    result = np.array(result)

    np.savetxt(f"config_cascade_{basin_n}_{model}_{metric}.csv", result, delimiter=",")

    print(f"\n\n\nFINISHED for {metric} using model {model}")

def execute_hydro_cro_full(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    objfunc = HydroFullModelGOF(rscript_name_cascade, data_file_cascade, basin_file_cascade, metric, model)
    substrates_real = [
        # SubstrateReal("Gauss", {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])}),
        SubstrateReal("Cauchy", {"F": np.tile(np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001]), len(basin_df)).flatten()}),
        SubstrateReal("MutNoise", {"method": "Gauss", "F": 1e-4, "N": 1}),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
        SubstrateReal("BLXalpha", {"F": 0.35}),
        SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
        # SubstrateReal("Perm", {"N": 2}),
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
        "Neval": 4e4,
        "fit_target": 1000,

        "verbose": True,
        "v_timer": 1,

        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 200,
        "prob_amp": 0.015
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_full_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")

def execute_hydro_cro_fullpon(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    objfunc = HydroFullModelGOF(rscript_name_cascade, data_file_cascade, basin_file_cascade, metric, model, weights=weights)
    substrates_real = [
        # SubstrateReal("Gauss", {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])}),
        SubstrateReal("Cauchy", {"F": np.tile(np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001]), len(basin_df)).flatten()}),
        SubstrateReal("MutNoise", {"method": "Gauss", "F": 1e-4, "N": 1}),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr":0.7}),
        SubstrateReal("BLXalpha", {"F": 0.35}),
        SubstrateReal("Firefly", {"a": 0.7, "b": 1, "d": 0.95, "g": 10}),
        # SubstrateReal("Perm", {"N": 2}),
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

        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 200,
        "prob_amp": 0.015
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_fullpon_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name+".eps")
    c.save_solution(output_name+".csv")

def execute_hydro_cro_wrapper(x):
    if x[0] == "single":
        execute_hydro_cro_single(*x[1:])
    elif x[0] == "cascade":
        execute_hydro_cro_cascade(*x[1:])
    elif x[0] == "full":
        execute_hydro_cro_full(*x[1:])
    elif x[0] == "fullpon":
        execute_hydro_cro_fullpon(*x[1:])
    else:
        raise Exception(f"Try using 'single', 'cascade', 'full' or 'fullpon' instead of {x[0]}.")

def main(args):
    pool = multiprocessing.Pool(processes=16)

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    args = [
        ('single', 'NSE', 0),
        ('single', 'NSE', 1),
        ('single', 'KGE', 0),
        ('single', 'KGE', 1),
        ('cascade', 'KGE', 1),
        # ('full', 'NSE', 0),
        # ('full', 'NSE', 1),
        # ('full', 'MSE', 1),
        # ('full', 'KGE', 1),
        # ('fullpon', 'NSE', 0),
        # ('fullpon', 'NSE', 1),
        # ('fullpon', 'MSE', 0),
        # ('fullpon', 'MSE', 1),
        # ('fullpon', 'KGE', 1),
    ]
    
    main(args)
