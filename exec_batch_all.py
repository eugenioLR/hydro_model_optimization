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
import argparse

parser = argparse.ArgumentParser(prog="River basin optimization", description="Optimizes river basins.")

parser.add_argument("-b", "--basin_n")
args = parser.parse_args()

if args.basin_n is None:
    raise Exception("Please indicate the basin with the -b flag. Available basins codes are 3005 and 5043.")

# basin_n = 5043  # 3005 or 5043
basin_n = int(args.basin_n)

rscript_name_single = "exec_optim.R"
rscript_name_cascade = "exec_optim_semidist.R"

if basin_n == 3005:
    data_file = "data/data_CHT_SIMPA_cal.csv"
    basin_file = "data/basins_CHT.csv"
elif basin_n == 5043:
    data_file = "data/data_CHG_SIMPA_cal.csv"
    basin_file = "data/basins_CHG.csv"
else:
    raise Exception("Invalid basin. Try 5043 or 3005")

r = robjects.r
r["source"](rscript_name_single)
r["source"](rscript_name_cascade)
get_basin_q = robjects.globalenv["eval_basin_param"]

weights = []
basin_df = pd.read_csv(basin_file)
basin_df.sort_values(by=["order"])

full_data = pd.read_csv(data_file)
for i, idx in enumerate(basin_df.index):
    basin_code = basin_df["code"][idx]
    total_caudal = full_data[full_data["qhmobs"] != -100][full_data["code"] == basin_code]["qhmobs"].sum()
    weights.append(total_caudal)
print(weights)
weights = np.asarray(weights) / sum(weights)
print(weights)


def execute_hydro_cro_single(metric, model):
    print(f"START for {metric} using model {model}\n\n")
    objfunc = HydroSimpleModelGOF(rscript_name_single, data_file, basin_file, metric, model)
    substrates_real = [
        SubstrateReal(
            "Cauchy",
            {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])},
        ),
        SubstrateReal(
            "MutNoise",
            {
                "method": "Gauss",
                "F": np.array([0.001, 1, 0.001, 0.1, 0.001, 0.01, 0.001]),
                "N": 1,
            },
        ),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr": 0.7}),
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
        "Ngen": 2,
        "Neval": 1e4,
        "fit_target": 1000,
        "verbose": True,
        "v_timer": 1,
        "Njobs": 1,
        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 500,
        "prob_amp": 0.015,
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_single_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name + ".eps")
    c.save_solution(output_name + ".csv")


def execute_hydro_cro_cascade(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    substrates_real = [
        SubstrateReal(
            "Cauchy",
            {"F": np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001])},
        ),
        SubstrateReal(
            "MutNoise",
            {
                "method": "Gauss",
                "F": np.array([0.001, 1, 0.001, 0.1, 0.001, 0.01, 0.001]),
                "N": 1,
            },
        ),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr": 0.7}),
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
        "Ngen": 2,
        "Neval": 1e4,
        "fit_target": 1000,
        "verbose": True,
        "v_timer": 1,
        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 500,
        "prob_amp": 0.015,
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

        objfunc = HydroSemidistModelGOF(
            rscript_name_cascade,
            data_file,
            basin_file,
            metric,
            model,
            basin_code,
            prev_q,
        )
        c = CRO_SL(objfunc, substrates_real, params)

        print(f"Optimizing basin with code {basin_code}.\n\n")

        c.safe_optimize()

        output_name = f"config_cascade_nd_{model}_{metric}_{basin_code}"

        c.display_report(show_plots=False, save_figure=True, figure_name=output_name + ".eps")

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

    objfunc = HydroFullModelGOF(rscript_name_cascade, data_file, basin_file, metric, model)
    substrates_real = [
        SubstrateReal(
            "Cauchy",
            {"F": np.tile(np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001]), len(basin_df)).flatten()},
        ),
        SubstrateReal(
            "MutNoise",
            {
                "method": "Gauss",
                "F": np.tile(np.array([0.001, 1, 0.001, 0.1, 0.001, 0.01, 0.001]), len(basin_df)).flatten(),
                "N": 1,
            },
        ),
        SubstrateReal("MutNoise", {"method": "Gauss", "F": 1e-4, "N": 1}),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr": 0.7}),
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
        "Ngen": 2,
        "Neval": 1e4,
        "fit_target": 1000,
        "verbose": True,
        "v_timer": 1,
        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 200,
        "prob_amp": 0.015,
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.safe_optimize()
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_full_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name + ".eps")
    c.save_solution(output_name + ".csv")


def execute_hydro_cro_fullpon(metric, model):
    print(f"START for {metric} using model {model}\n\n")

    objfunc = HydroFullModelGOF(
        rscript_name_cascade,
        data_file,
        basin_file,
        metric,
        model,
        weights=weights,
    )
    substrates_real = [
        SubstrateReal(
            "Cauchy",
            {
                "F": np.tile(
                    np.array([0.0001, 0.1, 0.0001, 0.01, 0.0001, 0.001, 0.0001]),
                    len(basin_df),
                ).flatten(),
            },
        ),
        SubstrateReal(
            "MutNoise",
            {
                "method": "Gauss",
                "F": np.tile(
                    np.array([0.001, 1, 0.001, 0.1, 0.001, 0.01, 0.001]),
                    len(basin_df),
                ).flatten(),
                "N": 1,
            },
        ),
        SubstrateReal("DE/best/2", {"F": 0.7, "Cr": 0.7}),
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
        "Ngen": 2,
        "Neval": 1e4,
        "fit_target": 1000,
        "verbose": True,
        "v_timer": 1,
        "dynamic": True,
        "dyn_method": "fitness",
        "dyn_metric": "best",
        "dyn_steps": 200,
        "prob_amp": 0.015,
    }
    c = CRO_SL(objfunc, substrates_real, params)

    c.optimize(save_data=False, update_data=False)
    print(f"\n\n\nFINISHED for {metric} using model {model}")

    output_name = f"config_fullpon_{basin_n}_{model}_{metric}"

    c.display_report(show_plots=False, save_figure=True, figure_name=output_name + ".eps")
    c.save_solution(output_name + ".csv")


def execute_hydro_cro_wrapper(x):
    strategy, metric, model_type = x
    if strategy == "single":
        execute_hydro_cro_single(metric, model_type)
    elif strategy == "cascade":
        execute_hydro_cro_cascade(metric, model_type)
    elif strategy == "full":
        execute_hydro_cro_full(metric, model_type)
    elif strategy == "fullpon":
        execute_hydro_cro_fullpon(metric, model_type)
    else:
        raise Exception(f"Try using 'single', 'cascade', 'full' or 'fullpon' instead of {strategy}.")


def main(args):
    pool = multiprocessing.Pool(processes=16)

    pool_results = pool.map_async(execute_hydro_cro_wrapper, args)

    pool_results.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    args_5043 = [
        ("single", "NSE", 0),
        ("cascade", "NSE", 0),
        ("full", "NSE", 0),
        ("fullpon", "NSE", 0),
        ("single", "NSE", 1),
        ("cascade", "NSE", 1),
        ("full", "NSE", 1),
        ("fullpon", "NSE", 1),
        ("single", "MSE", 0),
        ("cascade", "MSE", 0),
        ("full", "MSE", 0),
        ("fullpon", "MSE", 0),
        ("single", "MSE", 1),
        ("cascade", "MSE", 1),
        ("full", "MSE", 1),
        ("fullpon", "MSE", 1),
        ("single", "KGE", 0),
        ("cascade", "KGE", 0),
        ("full", "KGE", 0),
        ("fullpon", "KGE", 0),
        ("single", "KGE", 1),
        ("cascade", "KGE", 1),
        ("full", "KGE", 1),
        ("fullpon", "KGE", 1),
    ]

    args_3005 = [
        ("single", "NSE", 0),
        ("cascade", "NSE", 0),
        ("full", "NSE", 0),
        ("fullpon", "NSE", 0),
        ("single", "NSE", 1),
        ("cascade", "NSE", 1),
        ("full", "NSE", 1),
        ("fullpon", "NSE", 1),
        ("single", "MSE", 0),
        ("cascade", "MSE", 0),
        ("full", "MSE", 0),
        ("fullpon", "MSE", 0),
        ("single", "MSE", 1),
        ("cascade", "MSE", 1),
        ("full", "MSE", 1),
        ("fullpon", "MSE", 1),
        ("single", "KGE", 0),
        ("cascade", "KGE", 0),
        ("full", "KGE", 0),
        ("fullpon", "KGE", 0),
        ("single", "KGE", 1),
        ("cascade", "KGE", 1),
        ("full", "KGE", 1),
        ("fullpon", "KGE", 1),
    ]

    if basin_n == 5043:
        main(args_5043)
    else:
        main(args_3005)
