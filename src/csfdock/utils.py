# Some often used utilities..
from glob import glob
import os
import subprocess


def file_search(type=None, target="*", specific=None):
    """searches files in sub dir
    Args:
        type (str, optional): Search file format
        target (str, optional): Identifier to search
        specific (str, optional): Specific folder to search
    Returns:
        list: Search result
    """
    BASE_DIR = os.getcwd()
    try:
        if specific is None:
            return sorted(glob(f"{BASE_DIR}/**/{target}.{type}", recursive=True))
        else:
            return sorted(
                glob(f"{BASE_DIR}/**/{specific}/{target}.{type}", recursive=True)
            )
    except Exception as error:
        print(f"{error} \n File not found anywhere.")


def give_id(input_file):
    """Function to return the main file name excluding "." extension.
    Args:
        file (list): Name of file with "." extension.
    Returns:
        Name: Name without extension.
    """
    file_name = os.path.basename(input_file)
    file_name, file_format = file_name.rsplit(".")
    file_dir = os.path.dirname(input_file)
    return file_dir, file_name, file_format


def get(id, molecule="protein", prot_id="", type_="pdb"):
    """Downloads structure from RCSB and save in Generated sub folder.
    Args:
        id (TYPE): PDB ID
        molecule (str, optional): default:Protein or Small Molecule
        prot_id (str, optional): Description
        type_ (str, optional): Structure type to download.
    Returns:
        pdb/sdf/**: 3D coordinate file
    """
    try:
        assert molecule in [
            "protein",
            "ligand",
        ], 'Note: molecule parameter must be either "protein" or "ligand" only'
        if not os.path.exists("./Generated/"):
            os.makedirs("./Generated/")
        if molecule.lower() == "protein":
            assert type_ in [
                "pdb"
            ], "Note: \n Only PDB format supported for protein for now."
            command = (
                f"wget https://files.rcsb.org/download/{id}.{type_} -q -P ./Generated/"
            )
            msg = f"downloading of {id}.{type_}"
        elif molecule.lower() == "ligand":
            command = (
                f"wget -c -O https://files.rcsb.org/ligands/download/{id}_ideal.{type_}"
                f" > ./Generated/{prot_id}_{id}_ligand.{type_}"
            )
            msg = f"downloading of {prot_id}_{id}_ligand.{type_}"
        if os.path.exists(f"./Generated/{id}.{type_}"):
            msg = "but not downloaded as it already exists"
        else:
            subprocess.run(command, shell=True)
        return f"Succcesfully executed {msg} in ./Generated folder."
    except Exception as er:
        print(er)


def search_gui():
    output = widgets.Output()

    def f(File_type):
        global file_type
        file_type = File_type

    def search(file):
        output.clear_output()
        try:
            target = target_in.value
            if len(target) == 0:
                with output:
                    print("** Cannot be empty target!")
            else:
                global search_result
                specific = folder_in.value
                search_result = folder_search(
                    type=file_type, target=target, specific=specific
                )
                with output:
                    print(f"Total files found: {len(search_result)}")
                    print(search_result)
                    return search_result
        except:
            with output:
                print("something wrong")

    usage_information = widgets.HTML(
        "<b><font size = 3>Enter target file type and target name and specify"
        " the folder.</font size></b>"
    )
    display(usage_information)
    interact(f, File_type=["pdb", "sdf", "xyz", "txt"])
    target_in = widgets.Text(
        placeholder="Enter name", description="Target: ", disable=False
    )
    folder_in = widgets.Text(
        placeholder="specific folder name(optional)",
        description="Folder: ",
        disable=False,
    )
    search_button = widgets.Button(description="Search")
    search_button.style.button_color = "lightgreen"
    display(widgets.HBox([target_in, folder_in, search_button]), output)
    search_button.on_click(search)


def PDBParse(target):
    protein = []
    membrane = []
    tip3 = []
    ligand = []
    with open(target, "r") as target:
        temp = target.readlines()
        for line in temp:
            if line.startswith("ATOM") and line[21:22].strip() == "P":
                protein.append(line)
            elif line.startswith("ATOM") and line[21:22].strip() == "M":
                membrane.append(line)
            elif line.startswith("ATOM") and line[21:22].strip() == "T":
                tip3.append(line)
            elif line.startswith("HETATM"):
                ligand.append(line)

    return protein, membrane, tip3, ligand
