import random
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

from PyCROSL import AbsObjectiveFunc


class HydroSimpleModelGOF(AbsObjectiveFunc):
    def __init__(self, rscript_name='exec_optim.R', data_file="data/CHGdataSIMPA5043AG.txt", basin_file="data/CHGbasins5043AG.txt",
                 metric="MSE", model_used=0, basin_code=None):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)

        robjects.globalenv['init_global'](data_file, basin_file)

        inf_lim = np.array([1e-4,10,0,1e-4,1e-4,2,-6])
        sup_lim = np.array([1,2000,1,100,1,11,1])
        self.input_size = 7

        self.metric = metric
        self.model_used = model_used

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"
        else:
            opt = "min"

        # Loading the function we have defined in R.
        self.exec_function_r = robjects.globalenv['hydro_prob']

        super().__init__(self.input_size, opt, sup_lim, inf_lim)
        self.name = f"Hydro Simple GOF (metric={metric}, model_type={model_used}{', basin_code=' + str(basin_code) if basin_code else ''})"

    def objective(self, solution):
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
        return (self.sup_lim-self.inf_lim) * np.random.random(self.input_size) + self.inf_lim
    
    def repair_solution(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)


class HydroSemidistModelGOF(AbsObjectiveFunc):
# class HydroChainedModelGOF(AbsObjectiveFunc):
    def __init__(self, rscript_name='exec_optim_semidist.R', data_file="data/CHGdataSIMPA.txt", basin_file="data/CHGbasins.txt",
                 metric="MSE", model_used=0, basin_code=-1, prev_q=0):         
        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)
        robjects.globalenv['init_global'](data_file, basin_file)

        self.metric = metric
        self.model_used = model_used
        self.basin_code = basin_code
        self.prev_q = prev_q

        # Loading the function we have defined in R.
        self.exec_function_r = robjects.globalenv['eval_basin_param']
        inf_lim = np.array([1e-4,10,0,1e-4,1e-4,2,-6])
        sup_lim = np.array([1,2000,1,100,1,11,1])

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"
        elif metric in ["MSE", "RMSE"]:
            opt = "min"

        size = 7

        # self.input_size = 0
        super().__init__(size, opt, sup_lim, inf_lim)
        self.name = f"Hydro Semidist GOF (metric={metric}, basin={basin_code}, model_type={model_used})"

    def objective(self, solution):
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
        return (self.sup_lim-self.inf_lim) * np.random.random(self.input_size) + self.inf_lim
    
    def repair_solution(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)


class HydroFullModelGOF(AbsObjectiveFunc):
    def __init__(self, rscript_name='exec_semidist.R', data_file="data/CHGdataSIMPA.txt", basin_file="data/CHGbasins.txt",
                 metric="MSE", model_used=0, weights=None):#, basin_code=3005):

        # Defining the R script and loading the instance in Python
        r = robjects.r
        r['source'](rscript_name)
        robjects.globalenv['init_global'](data_file, basin_file)

        # Loading the function we have defined in R.
        self.get_basin_q = robjects.globalenv['get_basin_q']
        self.exec_function_r = robjects.globalenv['eval_basin']
        
        basin_df = pd.read_csv(basin_file)
        basin_df.sort_values(by=["order"])
        

        self.metric = metric
        self.model_used = model_used
        self.basin_df = basin_df
        # self.basin_code = basin_code

        if weights is None:
            weights = np.full(len(basin_df), 1/len(basin_df))
        self.weights = np.asarray(weights)

        inf_lim = np.tile(np.array([1e-4,10,0,1e-4,1e-4,2,-6]), len(basin_df))
        sup_lim = np.tile(np.array([1,2000,1,100,1,11,1]), len(basin_df))
        self.param_len = 7

        opt = "min"
        if metric in ["KGE", "NSE", "R2"]:
            opt = "max"
        elif metric in ["MSE", "RMSE"]:
            opt = "min"     

        size = 7*len(basin_df)

        super().__init__(size, opt, sup_lim, inf_lim)
        self.name = f"Hydro Full GOF (metric={metric}, basin={basin_df['code'].iloc[-1]}, model_type={model_used})"

    def objective(self, solution):
        agg_q = {}
        result_q = {}
        
        results = []

        for i, idx in enumerate(self.basin_df.index):
            basin_code = self.basin_df["code"][idx]
            codedown = self.basin_df["codedown"][idx]

            prev_q = 0
            if basin_code in agg_q:
                prev_q = agg_q[basin_code]
            solution_basin = solution[i*self.param_len:(i+1)*self.param_len]
            
            result_q[basin_code] = self.get_basin_q(self.model_used, solution_basin, basin_code, prev_q)
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

            # result += result_aux * self.weights[i]
            results.append(result_aux)
        results = np.asarray(results)
        return np.average(results, weights=self.weights)
        

    def random_solution(self):
        return (self.sup_lim-self.inf_lim) * np.random.random(self.input_size) + self.inf_lim
    
    def repair_solution(self, solution):
        return solution.clip(self.inf_lim, self.sup_lim)

if __name__ == "__main__":
    from os import system
    from ascii_magic import AsciiArt

    try:
        print(f"\n\nLaunching simple model optimizaion example\n\n")
        f1 = HydroSimpleModelGOF()
        s1 = f1.random_solution()
        err = f1.objective(s1)
        print(f"\n\nERROR OBTAINED IN F1: {err}\n\n")


        print(f"\n\nLaunching cascade model optimizaion example\n\n")
        f2 = HydroSemidistModelGOF(basin_code=5043)
        s2 = f2.random_solution()
        err = f2.objective(s2)
        print(f"\n\nERROR OBTAINED IN F2: {err}\n\n")


        print(f"\n\nLaunching Full model optimizaion example\n\n")
        f3 = HydroFullModelGOF()
        s3 = f3.random_solution()
        err = f3.objective(s3)
        print(f"\n\nERROR OBTAINED IN F3: {err}\n\n")

        my_art = AsciiArt.from_url("http://i0.kym-cdn.com/photos/images/facebook/001/244/891/d1f.png")
        my_art.to_terminal(columns=70, monochrome=False)
    except Exception as e:
        
        print(f"NOOOOOOOO. \n{e}\n")
        my_art = AsciiArt.from_url("https://cdn.frankerfacez.com/emoticon/620942/4")
        my_art.to_terminal(columns=70, monochrome=False)

    