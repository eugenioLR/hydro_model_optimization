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
    raise Exception("Please indicate the basin with the -b flag")

# basin_n = 5043  # 3005 or 5043
basin_n = int(args.basin_n)

rscript_name_single = "exec_optim.R"
rscript_name_cascade = "exec_optim_semidist.R"

if basin_n == 3005:
    data_file_single = "data/data_CHT_3005.txt"
    data_file_cascade = "data/data_CHT_SIMPA.txt"
    basin_file_single = "data/basins_CHT_3005.txt"
    basin_file_cascade = "data/basins_CHT.txt"
elif basin_n == 5043:
    data_file_single = "data/data_CHG_5043.txt"
    data_file_cascade = "data/data_CHG_SIMPA.txt"
    basin_file_single = "data/basins_CHG_5043.txt"
    basin_file_cascade = "data/basins_CHG.txt"
else:
    raise Exception("Invalid basin. Try 5043 or 3005")

r = robjects.r
r["source"](rscript_name_single)
r["source"](rscript_name_cascade)
get_basin_q = robjects.globalenv["get_basin_q"]

weights = []
basin_df = pd.read_csv(basin_file_cascade)
basin_df.sort_values(by=["order"])

full_data = pd.read_csv(data_file_cascade)
for i, idx in enumerate(basin_df.index):
    basin_code = basin_df["code"][idx]
    total_caudal = full_data[full_data["qhmobs"] != -100][full_data["code"] == basin_code]["qhmobs"].sum()
    weights.append(total_caudal)
print(weights)
weights = np.asarray(weights) / sum(weights)
print(weights)


def execute_hydro_cro_single():    
    r["source"](rscript_name_single)

    # Loading the function we have defined in R.
    robjects.globalenv["init_global"](data_file_single, basin_file_single)
    exec_function_r = robjects.globalenv["hydro_prob"]

    result_path = ""
    problem_type = "single"
    exec_types = [("NSE", 1), ("MSE", 1), ("KGE", 1)]

    result_files = [(f"./{result_path}config_{problem_type}_{basin_n}_{model}_{metric}.csv", model, metric) for metric, model in exec_types]

    eval_df = pd.DataFrame(columns=["model", "target", "MSE", "RMSE", "Pbias", "NSE", "R2", "KGE", "params"])
    params_list = []
    for file_name, model, target in result_files:
        # print(file_name, "Yes" if os.path.exists(file_name) else "No")
        params = np.loadtxt(file_name, delimiter=",", max_rows=1)
        params_list.append(",".join(map(str, params)))

        metrics = exec_function_r(model, params)
        param_str = np.array2string(params, max_line_width=np.inf, separator=";").replace(" ", "")
        eval_df.loc[len(eval_df)] = [model, target] + list(metrics) + [param_str]
    eval_df = eval_df[eval_df["target"] != "R2"]
    eval_df

    # code,order,codedown,supha,MSE,RMSE,PBIAS,NSE,R2,KGE,f0,param
    eval_df["code"] = 5043
    eval_df["order"] = 1
    eval_df["codedown"] = 0
    eval_df["supha"] = 329429.24
    eval_df["f0"] = -100
    eval_df["PBIAS"] = -100
    eval_table_df = eval_df[["code", "order", "codedown", "supha", "MSE", "RMSE", "PBIAS", "NSE", "R2", "KGE", "f0", "params"]]
    eval_table_df.loc[:, "params"] = eval_table_df.loc[:, "params"].str.replace("[", "")
    eval_table_df.loc[:, "params"] = eval_table_df.loc[:, "params"].str.replace("]", "")
    eval_table_df

    eval_table_df.to_csv(f"table_{problem_type}_{basin_n}.txt", index=False, quoting=csv.QUOTE_NONE)


def execute_hydro_cro_cascade():
    print("hello2")

def execute_hydro_cro_full():
    print("hello3")

def execute_hydro_cro_fullpon():
    print("hello4")

def main():
    execute_hydro_cro_single()
    execute_hydro_cro_cascade()
    execute_hydro_cro_full()
    execute_hydro_cro_fullpon()

if __name__ == "__main__":
    main()
