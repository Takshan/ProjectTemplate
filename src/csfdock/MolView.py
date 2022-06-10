import ipywidgets
import py3Dmol
from IPython.display import HTML, display
from ipywidgets import (
    FileUpload,
    IntSlider,
    fixed,
    interact,
    interactive,
    interactive_output,
    widgets,
    Layout,
)
from rich.console import Console

console = Console()
from rdkit import Chem
from rdkit.Chem import AllChem
from csfdock.utils import *


class MolView:
    """3D molecular view
    Args:
        mol (AllChem Obj): rdkit return object
        size (tuple, optional): window size to display 3d view
        style (str, optional): "strick | line | ribbon"
        surface (bool, optional): surface view
        opacity (float, optional): opacity of view
    Returns:
        py3dmol: 3d visual
    """

    def __init__(self, *args, **kwargs):
        self.molecule = kwargs.get("molecule")
        self.size = kwargs.get("size", (800, 600))
        self.style = kwargs.get("style", "stick")
        self.surface = kwargs.get("surface", False)
        self.opacity = kwargs.get("opacity", 0.5)
        self.conformers_list = []
        self.all_select = False
        self.default_name = False
        self.conformer = None
        # header information and display section
        self.msg_header = "Molecule Visualizer"
        self.msg_header_note = "Use upload option if want to import smiles from a file."
        self.msg_upload_info = (
            "Note: Upload option will copy your data on the server @"
            " ./Generated/smiles folder. "
        )
        self.msg_note = "Note: Use 4letter prefix similar to protein for ligand name"
        self.header = widgets.HTML(
            value=(
                "<b> <p style='text-align:center'><font size = 6"
                " > <text align = 'center'><font color"
                f" ='blue'>{self.msg_header}</text align></font"
                " color></font size></p></b>"
            )
        )
        self.header.add_class("header_bg")
        self.upload_info = widgets.HTML(
            value=f"{self.msg_upload_info}<br>{self.msg_note}"
        )
        self.output = widgets.Output()
        self.style_available = ["line", "stick", "sphere", "carton"]
        assert self.style in self.style_available, "Style Not Supported yet"
        # button sections
        self.add_button = widgets.Button(description="Add")
        self.add_button.style.button_color = "lightgreen"
        self.add_button.on_click(self.add_mol)
        self.remove_button = widgets.Button(description="Remove")
        self.remove_button.style.button_color = "salmon"
        self.remove_button.on_click(self.remove_mol)
        self.save_structure = widgets.Button(description="Save")
        self.save_structure.style.button_color = "lightblue"
        self.save_structure.on_click(self.write_mol)
        self.delete_button = widgets.Button(description="Delete")
        self.delete_button.style.button_color = "brown"
        self.delete_button.on_click(self.delete_mol)
        self.caption = widgets.Label(value="File Browser")
        self.upload_button = widgets.FileUpload(
            accept="", multiple=True, continuous_update=True
        )
        self.upload_button.observe(
            self.upload, names=["value", "content", "type", "name", "size"]
        )
        self.default_checkbox = widgets.Checkbox(
            value=False, description="Default names", disabled=False, indent=False
        )
        self.default_checkbox.observe(self.default, names="value")
        self.all_checkbox = widgets.Checkbox(
            value=False, description="Select All", disabled=False, indent=False
        )
        self.all_checkbox.observe(self.all_check, names="value")
        self.prefix_in = widgets.Text(
            placeholder="Enter atleast 4 letter name",
            description="Name: ",
            disable=False,
        )
        self.server_file_selected = widgets.Text(
            description="Selected File: ", disable=False
        )
        self.server_file_selected.on_submit(self.upload)
        self.prefix_in.on_submit(self.prefix_input)
        self.smile_in = widgets.Text(
            placeholder="Enter smile", description="Smile Code: ", disable=False
        )
        self.smile_in.on_submit(self.smile_input)
        self.INPUT_FLAG = False
        # Link accordian and file upload

    def __str__(self):
        return f"Total Molecules = {len(self.conformer)} {self.msg_header}"

    # widgets function section
    def smile_input(self, smi):
        self.output.clear_output()
        self.prefix_in.layout.visibility = None
        self.default_checkbox.layout.visibility = None
        self.smile_in.layout.visibility = "hidden"
        smile_name = self.smile_in.value if not isinstance(smi, list) else smi
        max = 0 if self.conformer is None else len(smi) - 1
        self.index_slider = IntSlider(
            value=0,
            min=0,
            max=max,
            step=1,
            disable=False,
            continuous_update=True,
            orientation="horizontal",
            layout=Layout(width="100%"),
        )
        smile_obj = interactive(
            self.smi_viewer,
            smile=smile_name,
            style=self.style_available,
            index=self.index_slider,
        )
        # self.smile_in.layout.visibility = "hidden"
        return display(smile_obj)

    def prefix_input(self, a):
        self.output.clear_output()
        prefix = self.prefix_in.value
        if len(prefix) < 4:
            with self.output:
                warnings.warn("Should be atleast 4 letter!.")
                exit()
        with self.output:
            print(
                f"Entered name: {prefix} .\nIf structure looks"
                " OK!\nYou can Save Now.\nElse delete it using Remove"
            )

    def add_mol(self, x):
        # self.output.clear_output()
        if self.INPUT_FLAG:
            self.molecule = self.conformer
        # TODO:  None when called through infunction call.
        try:
            if self.molecule in self.conformers_list:
                with self.output:
                    print("Already exist in the list.")
            else:
                self.conformers_list.append(self.molecule)
                with self.output:
                    # print(self.conformers_list)
                    console.print("Successfully Added!.")
        except:
            with self.output:
                print("Not able to Add")

    def remove_mol(self, y):
        self.output.clear_output()
        prefix = self.prefix_in.value
        try:
            self.conformers_list.remove(self.molecule)
            with self.output:
                console.print(" Successfully Remove from the temp list to save.!")
        except:
            with self.output:
                print("Molecule was not found!.")

    def delete_mol(self, m):
        # self.output.clear_output()
        prefix = self.prefix_in.value
        try:
            os.remove(f"./Generated/data/{prefix}.sdf")
            with self.output:
                print(f"./Generated/data/{prefix}.sdf Successfully deleted!.")
        except:
            with self.output:
                print(f"./Generated/data/{prefix}.sdf file not found!.")

    def write_all_mol(self, suffix, prefix):
        try:
            if all(isinstance(i, str) for i in self.conformer):
                conformers = [self.smile2conf(x) for x in self.conformer]
            elif isinstance(self.conformer, list):
                conformers = self.conformer
            conformers = list(filter(None, conformers))
            for conf in conformers:
                print(
                    Chem.MolToMolBlock(conf),
                    file=open(f"./Generated/data/{prefix}{suffix}.sdf", "w+"),
                )
                suffix += 1
            with self.output:
                print(f"{len(conformers)} molecules save in ./Generated/data folder. ")
            self.write = False
        except Exception as e:
            with self.output:
                print(
                    f"{e}\nCannot write all molecules\n"
                    "Presiding str maybe in the smiles code."
                )

    def default(self, m):
        # self.output.clear_output()
        self.default_name = m["new"]
        # print(f"Use default name: {default_name}")

    def all_check(self, n):
        # self.output.clear_output()
        self.all_select = n["new"]

    # Issue : invalid smiles upload saves the default "C"
    # TODO : smiles validity check and warn
    def upload(self, z):
        self.INPUT_FLAG = False
        if isinstance(z, str):
            with open(z, "r") as input_file_path:
                file_content = input_file_path.readlines()
                file_content = [x.strip() for x in file_content]
                input_file_dir, input_file_name, file_format = give_id(z)
                file_detail = f"{input_file_name}.{file_format}"
                self.INPUT_FLAG = True
                self.output = widgets.Output()
        else:
            try:
                file_detail = next(iter(self.upload_button.value))
                file_name, file_format = file_detail.rsplit(".", 1)
                file_content = self.upload_button.data
                file_content = [i.decode("utf-8") for i in file_content]
                file_content = "".join(str(i) for i in file_content)
                file_content = file_content.split()
            except StopIteration as er:
                input_file_dir, input_file_name, file_format = give_id(z.value)
                file_content = open(z.value).readlines()
                file_content = [x.strip() for x in file_content]
                file_detail = f"{input_file_name}.{file_format}"
                # print(file_content)
        # print(file_content)
        SDF = False
        if file_format.lower() == "sdf":
            # print(file_content)
            m = Chem.MolFromMolBlock(file_content)
            smiles = {}
            self.smile_in.layout.visibility = "hidden"
            self.view(m)
            SDF = True
        if self.INPUT_FLAG:
            smiles = file_content
        else:
            try:
                temp = [file_content.split("\r\n")]
                smiles = [
                    item
                    for subitem in temp
                    for item in subitem
                    if len(item.rstrip()) is not None
                ]
            except AttributeError:
                smiles = [item for item in file_content]
            # print(temp)

            smiles = list(filter(None, smiles))
        if not os.path.exists(f"./Generated/upload/"):
            os.makedirs(f"./Generated/upload/")
        self.all_checkbox.layout.visibility = None
        self.prefix_in.layout.visibility = None
        self.default_checkbox.layout.visibility = None
        try:
            with self.output:
                with open(
                    f"./Generated/upload/{file_detail}", "w+"
                ) as server_upload_file:
                    if SDF:
                        print(file_content, file=server_upload_file)
                        return
                    if isinstance(smiles, list):
                        for smi in smiles:
                            print(smi, file=server_upload_file)
            self.conformer = [s for s in smiles]
            # self.conformer = smiles
            self.smile_input(self.conformer)
            print(f"{file_detail} successfully uploaded.")
        except Exception as e:
            print(e)

    # todo // override default or use both custom and default..
    def write_mol(self, z):
        try:
            # self.output.clear_output()
            prefix = self.prefix_in.value
            self.write = True
            if len(prefix) < 4 and self.default_name is False:
                with self.output:
                    # warnings.warn("Check if name is entered or not and should be atleas 4 letter!.")
                    print(
                        "You need to add and enter either name or select default name."
                    )
            else:
                if not os.path.exists("./Generated/data"):
                    os.makedirs("./Generated/data")
                if len(prefix) < 4 and self.default_name is True:
                    suffix = 0
                    prefix = "small_molecule"
                    while os.path.exists(f"./Generated/data/{prefix}.sdf"):
                        suffix += 1
                        prefix = f"small_molecule{suffix}"
                    if self.all_select is False:
                        print(
                            Chem.MolToMolBlock(self.conformers_list[0]),
                            file=open(f"./Generated/data/{prefix}.sdf", "w+"),
                        )
                    else:
                        self.write_all_mol(suffix, prefix)
                elif len(prefix) >= 4 and self.default_name is True:
                    suffix = 0
                    prefix = f"{prefix}_small_molecule"
                    while os.path.exists(f"./Generated/data/{prefix}.sdf"):
                        suffix += 1
                        prefix = f"{prefix}_small_molecule{suffix}"
                    if self.all_select is False:
                        suffix += 1
                        prefix = f"{prefix}_small_molecule{suffix}"
                        print(
                            Chem.MolToMolBlock(self.conformers_list[0]),
                            file=open(f"./Generated/data/{prefix}.sdf", "w+"),
                        )
                        self.write = False
                    else:
                        self.write_all_mol(suffix, prefix)
                if self.write:
                    print(
                        Chem.MolToMolBlock(self.conformers_list[0]),
                        file=open(f"./Generated/data/{prefix}.sdf", "w+"),
                    )
                with self.output:
                    print(f"{prefix} saved in ./Generated/data/{prefix}.sdf.")
        except Exception as er:
            with self.output:
                print(
                    f"{er}\nSorry, check input!Name \nTip: Need to Add/Select All"
                    " first."
                )

    def display(self):
        """Displays widgets in jupyter notebook"""
        self.output.clear_output()
        display(self.header, self.upload_info)
        file_browser = ServerPath()
        display(widgets.VBox([self.caption, file_browser.accord]))
        display(self.server_file_selected)
        # Link upload and file browser
        display(
            (
                widgets.HBox(
                    (
                        self.add_button,
                        self.remove_button,
                        self.save_structure,
                        self.delete_button,
                        self.upload_button,
                        self.all_checkbox,
                    )
                )
            )
        )
        # elf.all_checkbox.layout.visibility = "hidden"
        # self.prefix_in.layout.visibility = "hidden"
        # self.default_checkbox.layout.visibility = "hidden"
        display(widgets.HBox([self.prefix_in, self.default_checkbox, self.smile_in]))
        display(self.output)
        # display(widgets.HBox([self.upload_button, self.all_checkbox]))
        # display(self.smile_in, self.output)
        # return self.view(self.molecule)

    def view(self, mol):
        """Creates py3dmol view
        Args:
            mol (rdkit obj): mol object from rdkit
        """
        try:
            molecular_block = Chem.MolToMolBlock(mol)
            viewer = py3Dmol.view(width=self.size[0], height=self.size[1])
            viewer.addModel(molecular_block, "mol")
            viewer.setStyle({self.style: {}})
            if self.surface:
                viewer.addSurface(py3Dmol.SAS, {"opacity": self.opacity})
            viewer.zoomTo()
            viewer.show()
        except Exception as er:
            with self.output:
                print(er, "in view section")

    def smi_viewer(self, smile, *args, **kwargs):
        """Converts smile to py3dmol view
        Args:
            smile (str): Valid smiles codes
            *args: Description
            **kwargs: style
        Returns:
            TYPE: Description
        """
        self.style = kwargs.get("style", "stick")
        index = kwargs.get("index", 0)
        # self.entered_smiles = kwargs.get("smiles")
        # self.entered_smiles = self.smile_in.value
        try:
            self.molecule = self.smile2conf(smile)
            print("+++++++++++++++++View+++++++++++++++++++++")
            print("Note: Hydrogens are added and MMFF Optimized.")
            # print(f"{Chem.MolToMolBlock(conf)}")
            print("++++++++++++++++++++++++++++++++++++++++++")
            # print(AllChem.EmbedMolecule(conf,randomSeed=0xf00d))
            return self.view(self.molecule)
        except Exception as er:
            with self.output:
                print(er, "pp")

    def smile2conf(self, smiles):
        """Convert SMILES to rdkit.Mol with 3D coordinates
        Args:
            smiles (str): smiles code
        Returns:
            AllChem.Mol: 3d mol object for visualization
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", f"{smiles}")
            if mol is None:
                return
            mol = Chem.AddHs(mol)
            # print(Chem.MolToMolBlock(mol))
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
            return mol
        except Exception as err:
            with self.output:
                print(err)


class ServerPath(MolView):
    def __init__(self, start_dir=".", select_file=True):
        super().__init__()
        self.file = None
        self.select_file = select_file
        self.cwd = start_dir
        self.select = ipywidgets.SelectMultiple(
            value=(), rows=10, description="", disabled=False
        )
        self.accord = ipywidgets.Accordion(children=[self.select])
        self.accord.selected_index = None  # Start closed (showing path only)
        self.refresh(".")
        self.select.observe(self.on_update, "value")
        # widget 1

    def on_update(self, change):
        if len(change["new"]) > 0:
            self.refresh(change["new"][0])

    def refresh(self, item):
        path = os.path.abspath(os.path.join(self.cwd, item))
        if os.path.isfile(path):
            if self.select_file:
                self.accord.set_title(0, path)
                self.file = path
                self.accord.selected_index = None
            else:
                self.select.value = ()
        else:  # os.path.isdir(path)
            self.file = None
            self.cwd = path
            # ipywidgets list of files and dirs
            keys = ["📁.."]
            for item in os.listdir(path):
                if item[0] == ".":
                    continue
                elif os.path.isdir(os.path.join(path, item)):
                    keys.append("📁" + item)
                else:
                    keys.append(item)
            # Sort and create list of output values
            keys.sort(key=str.lower, reverse=True)
            value = []
            for k in keys:
                if k[0] == "📁":
                    value.append(k[1:])  # strip off brackets
                else:
                    value.append(k)
            # Update widget
            self.accord.set_title(0, path)
            self.select.options = list(zip(keys, value))
            with self.select.hold_trait_notifications():
                self.select.value = ()
        if self.file is not None:
            # print(self.file)
            self.upload(self.file)
        # self.smile_in.layout.visibility = "hidden"
        # self.prefix_in.layout.visibility = None
        # self.default_checkbox.layout.visibility = None
        # return self.file
