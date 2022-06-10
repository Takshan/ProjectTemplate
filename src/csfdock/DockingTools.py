# Functions invloved in docking using smina

import itertools
import os
import re
import subprocess

import ipywidgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, display
from ipywidgets import (
    FileUpload,
    IntSlider,
    Layout,
    fixed,
    interactive,
    interactive_output,
    widgets,
)
from matplotlib.offsetbox import AnchoredText
from rich.console import Console
from scipy.special import expit
from sklearn.datasets import make_classification
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import KFold, train_test_split

from csfdock.utils import give_id

console = Console()


def add_hydrogen(list_xyz):
    """A function to add hydrogen atoms using
    openbabel to a list of molecules.
    Creates sub folder "protein_id/poses" in
    input file directory and dump there.
    Args:
        list_xyz (list): Molecules to be added Hydrogen
    Returns:
        xyz: saves xyz in protein_id/poses with suffix "_addHs.xyz"
    """
    # FOR SERVER or use own path..
    OBABEL_PATH = "/share/openbabel-3.1.1/bin/obabel"
    for i in list_xyz:
        dir_name, id, file_format = give_id(i)
        if not os.path.exists(f"{dir_name}/addH/{id[:4]}"):
            os.makedirs(f"{dir_name}/addH/{id[:4]}")
        command = f"{OBABEL_PATH} {i} -O {dir_name}/addH/{id[:4]}/{id}_addHs.xyz -h"
        subprocess.run(command, cwd=f"{dir_name}", shell=True)
    return "Successfully Completed."


def rmsd_calculator(reference, poses_list, key=None, nomatch=False, verbose=True):
    """Calculates RMSD between reference and pose using obrms in server.
    Args:
        reference (reference molecule): Reference molecule
        poses_list (test_poses): List of poses to test with the reference.
        key= Any suffix to add to the name of the file.
    Returns:
        Txt file: Writes reference, poses name and RMSD to a log file.
    """
    count = 0
    try:
        for ref in reference:
            _, ref_name, _dir = give_id(ref)
            ref_id = os.path.basename(ref)
            ref_dir = os.path.dirname(ref)
            ref_dir = os.path.dirname(ref_dir)
            for pose in poses_list:
                pose_id = os.path.basename(pose)
                if ref_id[:4].lower() == pose_id[:4].lower() or (nomatch is True):
                    command = f"/share/openbabel-3.1.1/bin/obrms {ref} \t {pose}"
                    pose_id = os.path.basename(pose)
                    rmsd_out = subprocess.run(
                        command,
                        cwd=f"{ref_dir}",
                        capture_output=True,
                        text=True,
                        shell=True,
                    )
                    if not os.path.exists(f"{ref_dir}/result/RMSD"):
                        os.makedirs(f"{ref_dir}/result/RMSD")
                    with open(
                        f"{ref_dir}/result/RMSD/{ref_id[:4]}_RMSD_{key}.txt", "a+"
                    ) as write_out:
                        rmsd_result = f"{pose_id.lower()}\t" + str(rmsd_out.stdout)
                        temp_rmsd = f"{ref_name}\t{pose_id.lower()}\t" + str(
                            float(rmsd_out.stdout.split()[-1])
                        )
                        print(rmsd_result, file=write_out)
                    count += 1
                else:
                    return "No matched header found"
        if verbose:
            print(f"Total of {count} rmsd calculated.")
        return temp_rmsd
    except Exception as er:
        print(er)


def rmsd_matrix_prep(rmsd_results, print_it=True, return_df=True, only_best=False):
    default_scoring_ = {}
    ad4_scoring_ = {}
    dkoes_fast_scoring_ = {}
    dkoes_scoring_old_scoring_ = {}
    vina_scoring_ = {}
    vinardo_scoring_ = {}
    custom_scoring_ = {}
    best = {}
    with open(rmsd_results[0], "r") as rmsd_results_read:
        for count, line in enumerate(rmsd_results_read):
            line_info = line.rsplit("/")[-1]
            try:
                key, value = line_info.split(".pdb")
            except Exception as er:
                # print(er)
                pass
            if "_ad4_scoring_" in line:
                ad4_scoring_[key] = value.strip()
            elif "_default_" in key:
                default_scoring_[key] = value.strip()
            elif "_dkoes_fast_" in key:
                dkoes_fast_scoring_[key] = value.strip()
            elif "_dkoes_scoring_old_" in key:
                dkoes_scoring_old_scoring_[key] = value.strip()
            elif "_vina_" in key:
                vina_scoring_[key] = value.strip()
            elif "_vinardo_" in key:
                vinardo_scoring_[key] = value.strip()
            else:
                custom_scoring_[key] = value.strip()
        # print(count)
    try:
        ad4_df = pd.DataFrame.from_dict(
            ad4_scoring_, orient="index", columns=(["ad4_scoring"])
        )
        default_df = pd.DataFrame.from_dict(
            default_scoring_, orient="index", columns=(["default_scoring"])
        )
        dkoes_fast_df = pd.DataFrame.from_dict(
            dkoes_fast_scoring_, orient="index", columns=(["dkoes_fast_scoring"])
        )
        dkoes_scoring_old_df = pd.DataFrame.from_dict(
            dkoes_fast_scoring_, orient="index", columns=(["dkoes_fast_scoring"])
        )
        vina_df = pd.DataFrame.from_dict(
            vina_scoring_, orient="index", columns=(["vina_scoring"])
        )
        vinardo_df = pd.DataFrame.from_dict(
            vinardo_scoring_, orient="index", columns=(["vinardo_scoring"])
        )
        custom_df = pd.DataFrame.from_dict(
            custom_scoring_, orient="index", columns=(["custom_scoring"])
        )
        return_df = pd.concat(
            [
                ad4_df,
                default_df,
                dkoes_fast_df,
                dkoes_scoring_old_df,
                vina_df,
                vinardo_df,
                custom_df,
            ]
        )
    except Exception as er:
        print(f"{er}\n error in data frame")
    if print_it:
        try:
            ad4_best = min(ad4_scoring_.items(), key=lambda x: x[1])
            print("===========BEST RMSD===================")
            best[ad4_best[0]] = ad4_best[1]
            print(f"{ad4_best[0]} : {ad4_best[1]}")
            default_best = min(default_scoring_.items(), key=lambda x: x[1])
            best[default_best[0]] = default_best[1]
            print(f"{default_best[0]} : {default_best[1]}")
            dkoes_fast_best = min(dkoes_fast_scoring_.items(), key=lambda x: x[1])
            best[dkoes_fast_best[0]] = dkoes_fast_best[1]
            print(f"{dkoes_fast_best[0]} : {dkoes_fast_best[1]}")
            dkoes_scoring_old_best = min(
                dkoes_scoring_old_scoring_.items(), key=lambda x: x[1]
            )
            best[dkoes_scoring_old_best[0]] = dkoes_scoring_old_best[1]
            print(f"{dkoes_scoring_old_best[0]} :{dkoes_scoring_old_best[1]}")
            vina_best = min(vina_scoring_.items(), key=lambda x: x[1])
            best[vina_best[0]] = vina_best[1]
            print(f"{vina_best[0]} : {vina_best[1]}")
            vinardo_best = min(vinardo_scoring_.items(), key=lambda x: x[1])
            best[vinardo_best[0]] = vinardo_best[1]
            print(f"{vinardo_best[0]} : {vinardo_best[-1]}")
            custom_best = min(custom_scoring_.items(), key=lambda x: x[1])
            best[custom_best[0]] = custom_best[1]
            print(f"{custom_best[0]} : {custom_best[1]}")
        except Exception as er:
            print(er)
    if only_best:
        return best
    if return_df is True:
        return return_df
    else:
        return (
            default_scoring_,
            ad4_scoring_,
            dkoes_fast_scoring_,
            dkoes_scoring_old_scoring_,
            vina_scoring_,
            vinardo_scoring_,
            custom_scoring_,
        )


def conformer_split(filenames, target):
    for i in filenames:
        file_dir, file_id, file_format = give_id(i)
        if not os.path.exists(f"{file_dir}/poses"):
            os.makedirs(f"{file_dir}/poses")
        command = (
            f"/share/openbabel-3.1.1/bin/obabel {i} -o{target} -O"
            f" ./poses/{file_id}_.{target} -m"
        )
        subprocess.run(command, cwd=f"{file_dir}", shell=True)
    return "Successfully completed."


def smina_histogram(sorted_RES, save=False):
    """Creates histogram of docking from smina output
    Args:
        sorted_RES (list): Sorted list of affinity values from smina output.
        save (bool, optional): Save plot in ./image/smina_histogram.png
    """
    name, benergy = zip(*sorted_RES.items())
    benergy = np.array((benergy), dtype=np.float32)
    mean = benergy.mean()
    best = benergy.min()
    worst = benergy.max()
    fig, axs = plt.subplots(1, sharey=False, sharex=False, tight_layout=True)
    axs.add_artist(
        AnchoredText(
            f"Total: {len(name)}\nMean: {mean:.2f}\nBest: {best:.2f}\nWorst:"
            f" {worst:.2f}",
            loc=1,
        )
    )
    axs.hist(benergy)
    axs.yaxis.set_label_text("Number of Datasets")
    axs.xaxis.set_label_text("Binding Energy Range")
    axs.set_title("Distribution of binding energy")
    if save:
        if not os.path.exists("./images"):
            os.makedirs("./images")
        plt.savefig("./images/smina_histogram.png", dpi=600)
    plt.show()


def smina_monitor(smina_output_monitor, plot=False, save=False):
    """Display smina process while enter smina stdout
    Args:
        smina_output_monitor (stdout): stdout from qstat/qsub
    Returns:
        dict: Displays result in jupyter
    """
    RES = {}
    count = 0
    for i in smina_output_monitor:
        algo_dir, algo_name, algo_format = give_id(i)
        if algo_format.lower() == "pdb":
            with open(i, "r") as read_smina:
                for line in read_smina:
                    if line[:5] == "MODEL":
                        number = line[5:].strip()
                    elif "REMARK" in line:
                        energy = float(line.rsplit(" ")[-1].strip())
                        RES[f"{algo_name}_{number}"] = f"{energy}"
                        count += 1
        elif algo_format.lower() == "sdf":
            pattern_id = r"^[a-zA-Z]\S"
            pattern_affinity = r"^>\s<[a-zA-Z]+>"
            with open(i, "r") as read_smina:
                write_affinity = False
                for line in read_smina:
                    if write_affinity:
                        energy = float(line.strip())
                        # print(energy)
                        RES[f"{algo_name}_{number}"] = f"{energy}"
                        count += 1
                        write_affinity = False
                    if re.match(pattern_id, line):
                        number = line.strip()
                    elif re.match(pattern_affinity, line):
                        write_affinity = True
        else:
            return "File format not supported yet. ['pdb', 'sdf']"

    print(f"Total number of poses generated: {count}")
    sorted_RES = dict(sorted(RES.items(), key=lambda x: x[1]))
    print("_____________Detail list________________\n")
    for key, value in sorted_RES.items():
        print(key, ":", value)
    if plot:
        smina_histogram(sorted_RES, save=save)

    return sorted_RES


def auc_plot(model, X_test, y_test, save=False):
    """Plot AUC plot from model, and X_text(data point) y_test(label).
    Args:
        model (sklearn obj)_: Model object from sklearn trained model.
        X_test (pd.DataFrame): Data point for test.
        y_test (pd.DataFrame): Label for the data point.
        save (bool, optional): Save AUC plot.
    """
    # assert isinstance(model, LinearRegression), f"{model} is not a valid model"
    # assert isinstance(X_test, pd.DataFrame), f"{X_test} is not a DataFrame"
    # assert isinstance(y_test, pd.DataFrame), f"{y_test} is not a DataFrame"
    model_regression_probability = model.predict_proba(X_test)
    model_regression_probability = model_regression_probability[:, 1]
    random_probability = [0 for _ in range(len(y_test))]
    random_auc = roc_auc_score(y_test, random_probability)
    model_auc = roc_auc_score(y_test, model_regression_probability)
    print(f"Random: ROC AUC={random_auc}")
    print(f"Model: ROC AUC={model_auc}")
    random_false_positive_rate, random_true_positive_rate, _ = roc_curve(
        y_test, random_probability
    )
    model_false_positive_rate, model_true_positive_rate, _ = roc_curve(
        y_test, model_regression_probability
    )
    plt.plot(
        random_false_positive_rate,
        random_true_positive_rate,
        linestyle="--",
        label="Random",
    )
    plt.plot(
        model_false_positive_rate,
        model_true_positive_rate,
        marker=".",
        label=f"Model (AUC:{model_auc:.2f})",
    )
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)
    plt.title("ROC")
    plt.legend()
    if save:
        prefix = "image"
        while os.path.exists(f"./Generated/images/{prefix}.png"):
            suffix += 1
            name = f"{prefix}{suffix}.png"
        plt.savefig("./Generated/images/{name}", dpi=600)
    plt.show()


def smina_model_score(
    file_path,
    num_features=3,
    intercept=False,
    tsize=0.3,
    plot_auc=False,
    plot_save=False,
):
    """Generates regression model using sklearn. Will Print out coefficients
    Args:
        file_path (str): csv/excel file path
        num_features (int, optional): Number of features to use. Default: 3
        intercept (bool, optional): Mean Error
        tsize (float, optional): Percentage of datato use for test. Default 0.3(30%)
        plot_auc (bool, optional): Plot ROC AUC curve
        plot_save (bool, optional): Save ROC AUC plot
    """
    try:
        if isinstance(file_path, str):
            file_dir, file_name, file_format = give_id(file_path)
            # for now | smina result file
            supported_file_format = ["csv", "excel"]
            assert file_format in supported_file_format, (
                "Note: FileType Error: Not supported file format. Use"
                f" {supported_file_format}"
            )
            if file_format == "excel":
                df = pd.read_excel(file_path)
            else:
                df = pd.read_csv(file_path)

    except Exception as e:
        print(e)

    if isinstance(file_path, pd.DataFrame):
        df = file_path
    try:
        X, y = df.iloc[1:, 1:-2], df.iloc[1:, -1]
        X = pd.DataFrame(X)
        header = X.iloc[:0, :]
        X, y = make_classification(n_features=num_features)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=tsize
        )  # TODO //include K-Fold Test
        model = LogisticRegression(fit_intercept=intercept)
        model.fit(X_train, y_train)  # TODOD accept use model input
        weight = model.coef_
        weight = [item for i in weight for item in i]
        console.print("[bold cyan]Model weights are :~[/bold cyan]\n")
        for head, coeff in zip(header, weight):
            print(coeff, head, end="\n")
        model.predict(X_test)
        model.predict_proba(X_test)
        score = model.score(X_test, y_test)
        print(f"\nModel score: {score}")
    except Exception as er:
        print(er)
    if plot_auc:
        auc_plot(model, X_test, y_test, save=plot_save)
        # train_plot(model, X, y, X_test, y_test)
    return model


def train_plot(model, X, y, X_test, y_test):
    plt.figure(1, figsize=(4, 3))
    plt.clf()
    print(len(X))
    print(len(y))
    plt.scatter(X.ravel(), y, color="black", zorder=20)
    # plt.scatter(y_test, X_test.iloc[:,0].values)
    loss = expit(X_test * model.coef_ + model.intercept_).ravel()
    plt.plot(X_test, loss, color="red", linewidth=3)

    ols = LinearRegression()
    ols.fit(X, y)
    plt.plot(X_test, ols.coef_ * X_test + ols.intercept_, linewidth=1)
    plt.axhline(0.5, color=".5")

    plt.ylabel("y")
    plt.xlabel("X")
    plt.xticks(range(-5, 10))
    plt.yticks([0, 0.5, 1])
    plt.ylim(-0.25, 1.25)
    plt.xlim(-4, 10)
    plt.legend(
        ("Logistic Regression Model", "Linear Regression Model"),
        loc="lower right",
        fontsize="small",
    )
    plt.tight_layout()
    plt.show()


def input_custom_scoring():
    """GUI window to enter custom scoring function"""
    # initialize some msg and output env
    output_csf = widgets.Output()
    msg_empty_name = "Enter any name for the file."
    warn_empty_name = widgets.HTML(value=f"<b><font color='red'>{msg_empty_name}</b>")
    # save the input custom scoring value

    def save_scoring(data):
        output_csf.clear_output()
        msg_confirm_warn = "Please confirm if the values are right."
        information = ipywidgets.widgets.HTML(
            value=f"<b><font color='red'>{msg_confirm_warn}</b>"
        )
        global scoring_data
        global temp_name
        temp_name = file_name.value
        if not temp_name:
            with output_csf:
                display(warn_empty_name)
        else:
            splitted = custom_scoring_area.value.split("\n")
            scoring_data = []
            for split in splitted:
                split = split.strip()
                scoring_data.append(split)
            with output_csf:
                print(f"Entered file name : {temp_name}")
                for line in scoring_data:
                    if select == "custom_scoring":
                        if len(line.rstrip()) != 0:
                            value, item = line.split()
                            print(f"{value}\t{item}")
                    else:
                        print(line)
                display(information)

    # writes the save scoring data to a file
    def confirm_scoring(data):
        output_csf.clear_output()
        if not temp_name:
            with output_csf:
                display(warn_empty_name)
        else:
            msg_success = "Confirmed and Saved!"
            information = ipywidgets.widgets.HTML(
                value=f"<b><font color='green'>{msg_success}</b>"
            )
            file_name_save = temp_name
            file_content = scoring_data
            if select == "custom_scoring":
                if not os.path.exists("./Generated/custom_function"):
                    os.makedirs("./Generated/custom_function")
                with open(
                    f"./Generated/custom_function/{file_name_save}_csf.txt", "w+"
                ) as write_scoring_function:
                    for line in file_content:
                        if len(line.rstrip()) != 0:
                            value, item = line.split()
                            print(f"{value}\t{item}", file=write_scoring_function)

                info = (
                    f"file save at ./Generated/custom_function/{file_name_save}_csf.txt"
                )
            else:
                if not os.path.exists("./Generated/smina_input"):
                    os.makedirs("./Generated/smina_input")
                with open(
                    f"./Generated/smina_input/{file_name_save}_mconfig.txt", "w+"
                ) as write_config:
                    for line in file_content:
                        print(f"{line}", file=write_config)

                info = (
                    f"file save at ./Generated/smina_input/{file_name_save}_mconfig.txt"
                )
            with output_csf:
                print(info)
                display(information)

    # def all_clear(data):
    #    with output_csf:
    #        output_csf.clear_output()
    # #text area to observe all input text
    config_placeholder = (
        "Paste here\n        \n        Sample config.txt Docking parameters file\n     "
        "   -------------------------------------\n        #Inputs\n        receptor ="
        " ./3L6B_prot.pdbqt\n        ligand = ./3L6B_lig.pdbqt\n        #Outputs\n     "
        "   out = 3L6B-nowat-Vina.pdbqt\n        log = 3L6B-nowat-Vina.log\n       "
        " #Box center\n        center_x = 4.500\n        center_y = -2.944\n       "
        " center_z = -5.250\n        #Box size\n        size_x = 50\n        size_y ="
        " 50\n        size_z = 50\n        #Parameters\n        exhaustiveness = 8\n   "
        "     seed = 123456\n"
    )
    csf_placeholder = (
        " Paste here\n        \n            Sample format of custom scoring\n       "
        " -------------------------------------\n        -0.035579   "
        " gauss(o=0,_w=0.5,_c=8)\n        -0.005156    gauss(o=3,_w=2,_c=8\n       "
        " 0.840245     repulsion(o=0,_c=8)\n        -0.035069   "
        " hydrophobic(g=0.5,_b=1.5,_c=8)\n        -0.587439   "
        " non_dir_h_bond(g=-0.7,_b=0,_c=8)\n        1.923        num_tors_div\n       "
        " -100.0       atom_type_gaussian(t1=Chlorine,t2=Sulfur,o=0,_w=3,_c=8)\n"
    )

    def evaluate(selected):
        output_csf.clear_output()
        area_layout = Layout(width="100%", height="400px", flex="row")
        global select
        select = selected
        if selected == "custom_scoring":
            global custom_scoring_area
            custom_scoring_area = widgets.Textarea(
                placeholder=csf_placeholder,
                description="Enter:",
                disabled=False,
                justify_content="space_between",
                continuous_update=True,
                layout=area_layout,
            )
        else:
            custom_scoring_area = widgets.Textarea(
                placeholder=config_placeholder,
                description="Enter:",
                disabled=False,
                justify_content="space_between",
                continuous_update=True,
                layout=area_layout,
            )
        display(custom_scoring_area)
        output_csf.clear_output()

    select_option = widgets.RadioButtons(
        options=["custom_scoring", "manual_config"],
        value="custom_scoring",
        description="What:",
        disabled=False,
    )
    ui = widgets.HBox([select_option])
    options = widgets.interactive_output(evaluate, {"selected": select_option})
    instruction = ipywidgets.widgets.HTML(
        "<font size = 4><b><font family: Times New Roman>Copy and Paste the scoring"
        " function below and enter</b></font-family></font size>"
    )
    display(instruction)
    # buttons widgets
    file_name = widgets.Text(description="Filename:", placeholder="file name ")
    save_button = widgets.Button(description="Save")
    save_button.style.button_color = "lightgreen"
    confirm_button = widgets.Button(description="Confirm")
    confirm_button.style.button_color = "salmon"
    # clear_button = widgets.Button(description="Clear")
    # clear_button.style.button_color = "lightgreen"
    display((widgets.VBox([file_name, ui, options])), output_csf)
    display((widgets.HBox([save_button, confirm_button])))
    save_button.on_click(save_scoring)
    confirm_button.on_click(confirm_scoring)
    # clear_button.on_click(all_clear)


def xg_model(X, y):
    from sklearn.datasets import make_classification

    num_classes = 3
    X, y = make_classification(n_samples=1000, n_informative=5, n_classes=num_classes)
    dtrain = xgb.DMatrix(data=X, label=y)
    num_parallel_tree = 4
    num_boost_round = 16
    # total number of built trees is num_parallel_tree * num_classes * num_boost_round

    # We build a boosted random forest for classification here.
    booster = xgb.train(
        {"num_parallel_tree": 4, "subsample": 0.5, "num_class": 3},
        num_boost_round=num_boost_round,
        dtrain=dtrain,
    )

    # This is the sliced model, containing [3, 7) forests
    # step is also supported with some limitations like negative step is invalid.
    sliced: xgb.Booster = booster[3:7]

    return [_ for _ in booster]


def xgb_boost(file_train, file_test):
    import xgboost as xgb

    CURRENT_DIR = os.path.dirname(__file__)
    dtrain = xgb.DMatrix(os.path.join(CURRENT_DIR, "file_train"))
    dtest = xgb.DMatrix(os.path.join(CURRENT_DIR, "file_test"))
    param = {
        "objective": "binary:logistic",
        "booster": "gblinear",
        "alpha": 0.0001,
        "lambda": 1,
    }
    watchlist = [(dtest, "eval"), (dtrain, "train")]
    num_round = 4
    bst = xgb.train(param, dtrain, num_round, watchlist)
    preds = bst.predict(dtest)
    labels = dtest.get_label()
    print(
        "error=%f"
        % (
            sum(int(preds[i] > 0.5) != labels[i] for i in range(len(preds)))
            / float(len(preds))
        )
    )


def run_smina(dir_name_, config_file_name, **kwargs):
    """ Creates folder within the cwd with the name id and sub\
    folder run where it will write sh.Also creates a dump folder\
    where error and output log will be dumped.

    """
    mode = kwargs.get("mode", False)
    log = kwargs.get("log", "log.txt")
    output = kwargs.get("output", "output.sdf")
    local = kwargs.get("local", False)
    cpu_num = kwargs.get("cpu", 2)
    job_name = kwargs.get("job_name", None)
    scoring = kwargs.get("scoring")
    custom = kwargs.get("custom", False)
    enter_output = kwargs.get("enter_output", True)
    enter_log = kwargs.get("enter_log", True)
    cluster = kwargs.get("cluster", None)
    cluster_grp = ["all.q", "gp1", "gp2"]
    if (cluster is not None) and (cluster not in cluster_grp):
        return f"Invalid cluster name. Available cluster names: {cluster_grp}"
    name_id = config_file_name[:4].lower()  # FIX
    dir_name_cwd = os.getcwd()
    dir_name = os.path.dirname(dir_name_)
    PATH = kwargs.get("PATH", False)
    dir_name = f"{dir_name}" if PATH else f"{dir_name_cwd}"

    if not os.path.exists(f"{dir_name}/Generated/jobs/{name_id}/run"):
        os.makedirs(f"{dir_name}/Generated/jobs/{name_id}/run")
    if local is False:
        SMINA_PATH = "/share/vina/smina"
        with open(
            f"{dir_name}/Generated/jobs/{name_id}/run/{name_id}_SMina.sh", "w"
        ) as out:
            if job_name is None:
                job_name = name_id
            if job_name[0].isdigit():
                job_name = "S" + job_name
            print(f"#$ -N {job_name}", file=out)
            print("#$ -V", file=out)
            print("#$ -S /bin/bash", file=out)
            if cluster is not None:
                print(f"#$ -q {cluster}", file=out)
            print(f"#$ -pe {cpu_num}cpu {cpu_num}", file=out)
            if not os.path.exists(f"{dir_name}/Generated/jobs/{name_id}/dump/"):
                os.makedirs(f"{dir_name}/Generated/jobs/{name_id}/dump/")
            print(f"#$ -o {dir_name}/Generated/jobs/{name_id}/dump/", file=out)
            print(f"#$ -e  {dir_name}/Generated/jobs/{name_id}/dump/", file=out)
            print("#$ -cwd", file=out)

            # Conditional to write log and output
            enter_log = f"--log {log}" if enter_log else ""
            enter_output = f"--out {output}" if enter_output else ""
            if (mode is True) and (custom is False):
                print(
                    f"{SMINA_PATH} --config"
                    f" {dir_name}/Generated/smina_input/{config_file_name}"
                    f" --scoring {scoring} --score_only {enter_log} {enter_output}",
                    file=out,
                )

            elif (mode is True) and (custom is True):
                print(
                    f"{SMINA_PATH} --config"
                    f" {dir_name}/Generated/smina_input/{config_file_name}"
                    f" --custom_scoring {scoring} --score_only  {enter_output}"
                    f" {enter_log}",
                    file=out,
                )

            elif (mode is False) and (custom is False):
                print(
                    f"{SMINA_PATH} --config"
                    f" {dir_name}/Generated/smina_input/{config_file_name} --scoring"
                    f" {scoring}  {enter_log} {enter_output} ",
                    file=out,
                )

            elif (mode is False) and (custom is True):
                print(
                    f"{SMINA_PATH} --config"
                    f" {dir_name}/Generated/smina_input/{config_file_name}"
                    f" --custom_scoring {scoring}  {enter_log} {enter_output} ",
                    file=out,
                )

        command = f"qsub {dir_name}/Generated/jobs/{name_id}/run/{name_id}_SMina.sh"
    else:
        if custom is False:
            command = (
                f"smina --config {dir_name}/Generated/smina_input/{config_file_name}"
                f" --scoring {scoring} {enter_log} {enter_output}"
            )
        else:
            command = (
                f"smina --config {dir_name}/Generated/smina_input/{config_file_name}"
                f" --custom_scoring {scoring}  {enter_log} {enter_output}"
            )

    subprocess.run(command, cwd=f"{dir_name}/Generated/jobs/{name_id}", shell=True)

    return "Succesfully completed."


# RECORDS OF ALL SCORING FUNCTION
# SF = [  "ad4_scoring",
#        "default",
#        "dkoes_fast",
#        "dkoes_scoring",
#        "dkoes_scoring_old",
#        "vina",
#        "vinardo",
#    ]

# CSF = ["custom"] # for custom scoring  and pass custom_scoring_file=PATH to the function
# Otherwise all the other SF will be run but will be calculate with csf
#  for numerous time


def smina_run(protein_list, ligand_list, **kwargs):
    """Prepares config file for smina when enter proterin and ligand
    Above sh_run function must be initialezed before in notebook"""

    SF = kwargs.get("SF")
    if not isinstance(SF, list):
        return "Supplied SF is not a list"
    cluster = kwargs.get("cluster", None)
    cluster_grp = ["all.q", "gp1", "gp2"]
    if (cluster is not None) and (cluster not in cluster_grp):
        return f"Invalid cluster name. Available cluster names: {cluster_grp}"
    autobox = kwargs.get("autobox", None)
    if isinstance(autobox, list):
        autobox = autobox[0]
    manual_config = kwargs.get("manual_config", None)
    if isinstance(manual_config, list):
        manual_config = manual_config[0]
    run = kwargs.get("run", False)
    mode = kwargs.get("mode", False)
    local = kwargs.get("local", False)
    job_name = kwargs.get("job_name", None)
    NUM_MODES = kwargs.get("num_modes", 10)
    EXHAUSTIVE = kwargs.get("exhaustive", 50)
    ENERGY_RANGE = kwargs.get("energy_range", 10)
    SEED = kwargs.get("seed", None)
    AUTOBOX_PAD = kwargs.get("pad", 4)
    CPU_NUM = kwargs.get("cpu", 8)
    nomatch = kwargs.get("match", False)
    OUT_FORMAT = kwargs.get("out_format", "sdf")
    custom = kwargs.get("custom", False)

    # if CUSTOM_SCORE is not None:
    #    CSF_FLAG = True
    # else:
    #    CSF_FLAG = False

    def write_config(
        receptor,
        ligand,
        config_file_name,
        output_file_name,
        log_file_name,
        # scoring=None, # Moved to CLI
    ):
        """Writes config into new files, if already \
        exist append to it."""

        # dir_name = os.path.dirname(receptor)
        dir_name = os.getcwd()

        if not os.path.exists(f"{dir_name}/Generated/smina_input"):
            os.makedirs(f"{dir_name}/Generated/smina_input/")
        with open(
            f"{dir_name}/Generated/smina_input/{config_file_name}", "w+"
        ) as config_file:

            # required config arguments
            # ------------------------
            print(f"receptor = {receptor} ", file=config_file)
            print(f"ligand = {ligand}", file=config_file)
            if autobox is not None:
                print(f"autobox_ligand = {autobox}", file=config_file)
                print(f"autobox_add = {AUTOBOX_PAD}", file=config_file)

            # Optionals con  # MOVED TO CLI
            # ------------------------------
            print(f"out = {output_file_name}", file=config_file)
            print(f"log = {log_file_name}", file=config_file)
            # print(f"scoring = {scoring}", file=config_file)  ## change to run in CLI

            # Misc(optional) configs
            # -----------------------------
            # if CUSTOM_SCORE is not None:  ## MOVED TO CLI
            #    print(f"custom_scoring = {CUSTOM_SCORE}", file=config_file)
            # if SMINA_MODE is not None:
            #    print(f"{SMINA_MODE}", file=config_file)

            print(f"cpu = {CPU_NUM}", file=config_file)
            if SEED is not None:
                print(f"\n\nseed = {SEED}", file=config_file)
            print(f"exhaustiveness = {EXHAUSTIVE}", file=config_file)
            # if CUSTOM_SCORE is None: ## change to run in CLI
            print(f"num_modes = {NUM_MODES}", file=config_file)
            print(f"energy_range = {ENERGY_RANGE }", file=config_file)

    for scoring in SF:
        for protein in protein_list:
            for ligand in ligand_list:
                protein_dir, protein_id, prot_format = give_id(protein)
                ligand_dir, ligand_id, lig_format = give_id(ligand)
                if custom:
                    _dir, _name, _format = give_id(scoring)
                    scoring = _name
                if protein_id[:4].lower() == ligand_id[:4].lower() or (nomatch == True):
                    output_file_name = ligand_id + f"_output_{scoring}.{OUT_FORMAT}"
                    log_file_name = ligand_id + f"_log_{scoring}.txt"
                    if manual_config is None:
                        config_file_name = ligand_id.lower() + f"_config_{scoring}.txt"
                        enter_output = True
                        enter_log = True
                        write_config(
                            protein_list[0],
                            ligand,
                            config_file_name,
                            output_file_name,
                            log_file_name,
                        )
                    else:
                        (
                            manual_config_dir,
                            manual_config_name,
                            manual_config_format,
                        ) = give_id(manual_config)
                        config_file_name = (
                            f"{manual_config_name}.{manual_config_format}"
                        )
                        if "output" in open(manual_config).read():
                            enter_output = False
                        else:
                            enter_output = True
                        if "log" in open(manual_config).read():
                            enter_log = False
                        else:
                            enter_log = True
                else:
                    print("Protein and ligand prefix [4 letter] didnt match")
                if custom:
                    scoring = f"{_dir}/{_name}.{_format}"
                if run is True and mode is True:
                    run_smina(
                        protein_dir,
                        config_file_name,
                        scoring=scoring,
                        mode=mode,
                        local=local,
                        # cpu=CPU_NUM,
                        # job_name=job_name,
                        custom=custom,
                        log=log_file_name,
                        output=output_file_name,
                        enter_output=enter_output,
                        enter_log=enter_log,
                    )
                elif run is True and mode is False:
                    run_smina(
                        protein_dir,
                        config_file_name,
                        scoring=scoring,
                        local=local,
                        cpu=CPU_NUM,
                        job_name=job_name,
                        custom=custom,
                        output=output_file_name,
                        log=log_file_name,
                        enter_output=enter_output,
                        enter_log=enter_log,
                        PATH=True,
                        cluster=cluster,
                    )
                else:
                    print(
                        "Run command was not passed so only created the"
                        " config and sh file but not executed"
                    )
            # print(base_name)
    return "Succesfully completed."


def view_affinity(sorted_RES, keyword):
    return pd.DataFrame(
        [
            (key, value)
            for key, value in sorted_RES.items()
            if f"{keyword}" in key.lower()
        ],
        columns=["Pose", "Affinity"],
    )


def smina_output_df(sorted_RES):
    return pd.DataFrame(
        [(key, value) for key, value in sorted_RES.items()],
        columns=["Pose", "Affinity"],
    )


def rmsd_matrix(
    ref, length=2, key="MATRIX", verbose=False, plot=True, save=False, annot=False
):
    _rmsd_list = itertools.combinations(ref, length)
    cols = ["Reference", "Pose", "RMSD"]
    df = pd.DataFrame(columns=cols)
    for count, i in enumerate(_rmsd_list):
        ref, conf = list([i[0]]), list([i[1]])
        x = rmsd_calculator(ref, conf, nomatch=True, key=key, verbose=False)
        if verbose:
            print(f"{count}. {x}", sep=" ", flush=True)
        x1, x2, x3 = x.split("\t")
        df.loc[count] = [x1, x2, x3]

    mdf = df.pivot(index="Reference", columns="Pose", values="RMSD")
    mdf.fillna(0)
    mdf = mdf.astype(float)
    if plot:
        # plt.figure(figsize=[15, 8])
        hmap = sns.heatmap(mdf, annot=annot)
        hmap.set_title("RMSD MATRIX")
        if save:
            if not os.path.exists("./Generated/images/"):
                os.makedirs("./Generated/images/")
            fig = hmap.figure
            image = f"./Generated/images/{key}.jpg"
            fig.savefig(image, dpi=600)
            print(f"Saved ! {image}")
    return mdf
