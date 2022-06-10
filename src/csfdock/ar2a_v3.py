import re
import sys
from collections import Counter

import ipywidgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import py3Dmol
from IPython.display import HTML, display
from ipywidgets import (
    FileUpload,
    IntSlider,
    fixed,
    interactive,
    interactive_output,
    widgets,
    Layout,
)
from ipywidgets.embed import embed_minimal_html
from matplotlib.offsetbox import AnchoredText
from rdkit import Chem
from rdkit.Chem import AllChem
from rich.console import Console
from rich.table import Table

console = Console()
# from rich import print

from sklearn.datasets import make_classification
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import KFold, train_test_split

from csfdock.DockingTools import *
from csfdock.DVisualize import *
from csfdock.KinaseModules import *
from csfdock.MolView import *
from csfdock.Project import *
from csfdock.utils import *
from csfdock.xg_mod import *

#


def view(structure, ligand=None, color="grey", save=False):
    """3d visualization of pdb
    Args:
        structure (TYPE): Description
        ligand (None, optional): small molecule
        color (str, optional): color of wish, default: grey
    Returns:
        TYPE: structure view.
    """
    structure_dir, structure_name, structure_format = give_id(structure)
    v = py3Dmol.view(width=900, height=500)
    if structure_format.lower() == "sdf":
        mol = Chem.MolFromMolFile(structure, removeHs=False)
        mol = Chem.MolToMolBlock(mol)
        v.addModel(mol, f"{structure_format}")
    else:
        v.addModel(open(structure).read())
    v.setStyle({"cartoon": {"color": f"{color}"}})
    if ligand is not None:
        v.setStyle({"resn": f"{ligand}"}, {"stick": {"colorscheme": "greenCarbon"}})
    v.zoomTo()
    v.show()
    if save:
        prefix = "image"
        while os.path.exists(f"./images/{prefix}.png"):
            suffix += 1
            name = f"{prefix}{suffix}.png"
        v.save_fig(f"./Images/{name}", dpi=600)
    return structure


def update_exp_data(new_data):
    """Enter new data to already generated experimental data
    Args:
        new_data (list|dict): New experimental data
    Returns:
        pd.DataFrame: Latest data
    """
    try:
        order_list = ["Elec", "Vdw", "exp"]
        if not new_data:
            return "Enter valid data"
        if isinstance(new_data, list):
            CONFIRMED = input(f"Is the list in order(yes|no)\n{order_list}: ")
            if CONFIRMED.lower() != "yes":
                return "Enter valid order experimental data"
            df = pd.DataFrame(data=new_data)
            df = df.T
            new_columns = {0: "Elec", 1: "Vdw", 2: "exp"}
            df.rename(columns=new_columns, inplace=True)
        else:
            df = pd.DataFrame(data=new_data)
        print(
            "[bold magenta]Staged for updating previous data with[/bold magenta]"
            f" \n{df}\n"
        )
        old_df = pd.read_pickle("./DATA/experimental_data.pickle")
        latest_df = pd.concat([old_df, df], ignore_index=True)
        print(
            "[bold green]Successfully save!! [/bold green]\n\nlatest experimental_data"
            f" :\n {latest_df}"
        )
        latest_df.to_pickle("./DATA/experimental_data.pickle")
        return latest_df
    except Exception as e:
        print(e)
        return


def parse_log(file):
    """Parse LIE log file in return delta Vwd and Elec
    Args:
        file (str): Log file path
    Returns:
        pd.DataFrame: Vdw and Elect DataFrame.
    """
    try:
        with open(file, "r") as file:
            info = []
            lines = file.readlines()
            extract = False
            for index, line in enumerate(lines):
                # print(f" {line.strip()}" )
                if line[:6].strip() == "Delta":
                    extract = True
                if extract and (
                    line[:6].strip() == "Vdw" or line[:7].strip() == "Elec"
                ):
                    info.append(line.split())
    except Exception as e:
        print(e)
    df = pd.DataFrame(info)
    df = df.T.reset_index(drop=True)
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace=True)
    return df


def exp_model_score(file_path, num_features=2, intercept=False, tsize=0.3, plot=False):
    """Generates regression model using sklearn. Will Print out coefficients
    Args:
        file_path (str): csv/excel file path
        num_features (int, optional): Number of features to use. Default: 3
        intercept (bool, optional): Mean Error
        tsize (float, optional): Percentage of datato use for test. Default 0.3(30%)
        plot_auc (bool, optional): Plot ROC AUC curve
        plot_save (bool, optional): Save ROC AUC plot
    """
    data, file_name, file_format = give_id(file_path)
    supported_file_format = ["csv", "excel"]  # for now | smina result file
    assert (
        file_format in supported_file_format
    ), f"Note: FileType Error: Not supported file format. Use {supported_file_format}"
    try:
        if file_format == "excel":
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path)
        print(df)
        X, y = df.iloc[:, :-1], df.iloc[:, -1]
        X = pd.DataFrame(X)
        header = X.iloc[:0, :]
        # print(f"----------\n {X}")
        # X, y = make_classification(n_samples=df.shape()[0],n_features=num_features)
        # TODO //include K-Fold Test
        model = LinearRegression(fit_intercept=intercept)
        model.fit(X, y)  # TODOD accept use model input
        print(f"Alpha: {model.coef_[-1]}, Beta= {model.coef_[0]}")
        # weight = model.coef_
        # weight = [item for i in weight for item in i]
        # for head, coeff in zip(header, weight):
        #    print(coeff, head, end="\n")
    except Exception as er:
        print(er)
    if plot:
        plt.scatter(X, y, color="black")
        plt.plot(X, y, color="blue", linewidth=3)
        plt.xticks(())
        plt.yticks(())
        plt.show()
    return model
