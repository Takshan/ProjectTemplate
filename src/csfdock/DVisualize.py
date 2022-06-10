import os

import py3Dmol
from ipywidgets import Layout, interactive

from rdkit import Chem
from csfdock.utils import give_id, PDBParse
from rich.console import Console
console = Console()
blue_console = Console(style="white on blue")


class DVisualize:
    """Grid bix view of the docking pocket
    Attributes:
        box_center_x (int): Coordinates of box center x-axis.
        box_center_y (int): Coordinates of box center y-axis.
        box_center_z (int): Coordinates of box center z-axis.
        box_size_x (int): Size of x-axis of grid box
        box_size_y (int): Size of x-axis of grid box
        box_size_z (int): Size of x-axis of grid box
        ligand (str): Path of ligand
        protein (str): Path of Receptor
    """

    def __init__(self, *args, **kwargs):
        self.grid_box_centers = None
        self.receptor = kwargs.get("protein", None)
        self.ligand = kwargs.get("ligand", None)
        self.box_center_x = kwargs.get("box_center_x")
        self.box_center_y = kwargs.get("box_center_y")
        self.box_center_z = kwargs.get("box_center_z")
        self.box_size_x = kwargs.get("box_size_x", 20)
        self.box_size_y = kwargs.get("box_size_y", 20)
        self.box_size_z = kwargs.get("box_size_z", 20)
        self.prot_color = kwargs.get("prot_color", "spectrum")
        self.lig_color = kwargs.get("lig_color", "red")
        self.membrane = kwargs.get("membrane", None)
        self.save = kwargs.get("save", "False")
        self.bg_color = kwargs.get("bg_color", "white")
        self.mem_color = kwargs.get("mem_color", "blue")
        for arg in args:
            if isinstance(arg, str):
                if self.receptor is None:
                    self.receptor = arg
            elif isinstance(arg, list):
                if self.grid_box_centers is None:
                    self.grid_box_centers = arg
                    self.box_center_x = arg[0]
                    self.box_center_y = arg[1]
                    self.box_center_z = arg[2]
                else:
                    self.grid_box_sizes = arg
                    self.box_size_x = arg[0]
                    self.box_size_y = arg[1]
                    self.box_size_z = arg[2]

    def LoadBox(self, *args, **kwargs):
        try:
            self.grid_box_centers
        except AttributeError:
            self.grid_box_centers = None
        for arg in args:
            if isinstance(arg, str):
                if self.receptor is None:
                    self.receptor = arg
            elif isinstance(arg, list):
                if self.grid_box_centers is None:
                    self.grid_box_centers = arg
                    self.box_center_x = arg[0]
                    self.box_center_y = arg[1]
                    self.box_center_z = arg[2]
                else:
                    self.grid_box_sizes = arg
                    self.box_size_x = arg[0]
                    self.box_size_y = arg[1]
                    self.box_size_z = arg[2]

    def __rep__(self):
        return f"Complex_Grid: {self.receptor} and {self.ligand}"

    def __str__(self):
        return f"Protein: {self.receptor} and \nligand :{self.ligand}"

    def __grid_box(self):
        try:
            self.vobj.addBox(
                {
                    "center": {
                        "x": self.box_center_x,
                        "y": self.box_center_y,
                        "z": self.box_center_z,
                    },
                    "dimensions": {
                        "w": self.box_size_x,
                        "h": self.box_size_y,
                        "d": self.box_size_z,
                    },
                    "color": "blue",
                    "opacity": 0.5,
                }
            )
        except Exception as e:
            print("Failed to add Grid")

    def LoadLipid(self, *args, verbose=True,native=False, **kwargs):
        lipid = kwargs.get("lipid")

        for arg in args:
            lipid = arg
        lipid_path = self.LoadReceptor(
            lipid, key="Lipid", verbose=verbose, native=native
        )
        _, lipid, water, lig = PDBParse(lipid_path)
        self.lipid = lipid_path
        with open("./temp.pdb", "w+") as f:
            for i in lipid:
                print(i, end="", file=f)
        # m = Chem.MolFromPDBFile("./temp.pdb", sanitize=False)
        # print(m)
        try:
            if self.vobj:
                pass
        except AttributeError:
            self.vobj = py3Dmol.view(width=800, height=600)
        lipid_mol = open("./temp.pdb").read()
        self.vobj.addModel(lipid_mol, "pdb")
        self.vobj.setStyle({"lipid_mol": 2}, {"cartoon": {}})
        try:
            os.remove("./temp.pdb")
        except Exception:
            pass

        # self.vobj.addModel(lipid, "pdb")
        # self.vobj.setStyle({"model": 3}, {"cartoon": {}})
        # self.vobj.setStyle({"cartoon": {"color": "spectrum"}})

    def __complex_view(self):
        mol1 = open(self.receptor, "r").read()
        file_format = "pdb"
        try:
            mol2 = open(self.ligand, "r").read()
            lig_dir, lig_name, lig_file_format = give_id(self.ligand)
            if lig_file_format == "sdf":
                file_format = "sdf"
            self.vobj.addModel(mol2, f"{file_format}")
            self.vobj.setStyle({"model": 1}, {"stick": {}})
        except TypeError as er:
            self.mol_view.setStyle(
                {"resn": f"{self.resn}"}, {"stick": {"colorscheme": self.lig_color}}
            )
        self.vobj.addModel(mol1, "pdb")
        self.vobj.setStyle({"cartoon": {"color": self.prot_color}})

    def __visualize_mol(self):
        self.vobj = py3Dmol.view(width=800, height=600)
        self.__grid_box()
        self.__box_view()
        try:                                                                
            self.LoadLipid(self.lipid, verbose=False, native=False)
        except AttributeError as er:
            blue_console.print("Lipid maynot be loaded yet")
        try:
             _ = self.bg_color
        except AttributeError:
            self.bg_color= "white" 
        self.vobj.setBackgroundColor(self.bg_color)     
        self.vobj.rotate(90, {"x": 0, "y": 1, "z": 0}, viewer=(0, 1))
        self.vobj.zoomTo()
        return self.vobj.show()

    def ShowMolecules(self, **kwargs):
        """Visualize grid box with protein complex
        Returns:
            py3dmol : 3D Viewer
        """
        self.resn = kwargs.get("resn", "LIG")

        grid_obj = interactive(self.__visualize_mol)
        return display(grid_obj)

    def __show_ligand(
        self, mol_view_object, mol, resn=None, mol_color="blue", style="stick"
    ):
        _, mol_name, mol_format = give_id(mol)

        try:
            mol2 = open(mol, "r").read()
            *_, mol_file_format = give_id(self.ligand)
            mol_view_object.addModel(mol2, f"{mol_file_format}")
            mol_view_object.setStyle({"model": 1}, {"stick": {}})
        except (TypeError, AttributeError) as er:
            print(
                # "Cannot.."
                "Searching name space..."
            )
            mol_view_object.setStyle(
                {"resn": f"{resn}"}, {f"{style}": {"colorscheme": mol_color}}
            )
            # print(er)
        return mol_view_object

    def SimpleView(self, **kwargs):
        """3d visualization of pdb
        Args:
            protein (TYPE): protein
            ligand (None, optional): small molecule
            color (str, optional): color of wish, default: grey
            resn (str): Ligand from pdb file.
        Returns:
            TYPE: structure view.
        """
        resn = kwargs.get("resn", "LIG")
        self.bg_color = kwargs.get("bg_color", "white")
        self.prot_color = kwargs.get("prot_color", "spectrum")
        self.lig_color = kwargs.get("lig_color", "red")
        self.save = kwargs.get("save", False)
        self.show_ligand = kwargs.get("show_ligand", True)
        self.show_receptor = kwargs.get("show_receptor", True)
        vobj = py3Dmol.view(width=900, height=500)
        vobj.setBackgroundColor(self.bg_color)
        if self.show_receptor:
            structure_dir, structure_name, structure_format = give_id(self.receptor)
            if structure_format.lower() == "sdf":
                mol = Chem.MolFromMolFile(self.receptor, removeHs=False)
                mol = Chem.MolToMolBlock(mol)
                vobj.addModel(mol, f"{structure_format}")
                self.clean.addModel(mol, f"{structure_format}")
            else:
                vobj.addModel(open(self.receptor).read())
                self.clean = vobj
        vobj.setStyle({"cartoon": {"color": f"{self.prot_color}"}})

        if self.show_ligand:
            try:
                self.__show_ligand(vobj, self.ligand, mol_color=self.lig_color)
            except AttributeError as er:
                print("Ligand not yet added to the project...")
        vobj.zoomTo()
        if self.save == True:
            prefix = "image"
            while os.path.exists(f"./images/{prefix}.png"):
                suffix += 1
                name = f"{prefix}{suffix}.png"
            vobj.save_fig(f"./Images/{name}", dpi=600)
            print(f"Successfully saved ./Images/{name} ")
        if self.show_ligand is False and self.show_receptor is False:
            return "Nothing to visualize.."

        return vobj.show()

    def __box_view(self, **kwargs):
        """3d visualization of pdb
        Args:
            protein (TYPE): protein
            ligand (None, optional): small molecule
            color (str, optional): color of wish, default: grey
            resn (str): Ligand from pdb file.
        Returns:
            TYPE: structure view.
        """
        self.resn = kwargs.get("resn", "LIG")
        self.membrane = kwargs.get("membrane", None)
        self.lig_color = kwargs.get("resn_color", "yellow")
        self.element = kwargs.get("element", None)
        self.save = kwargs.get("save", False)
        self.mem_color = kwargs.get("mem_color", "blue")
        file_format = "pdb"
        try:
            structure_dir, structure_name, structure_format = give_id(self.receptor)
            if structure_format.lower() == "sdf":
                mol = Chem.MolFromMolFile(self.receptor, removeHs=False)
                mol = Chem.MolToMolBlock(mol)
                self.vobj.addModel(mol, f"{structure_format}")
            else:
                self.vobj.addModel(open(self.receptor).read())
                                                           
            self.vobj.setStyle({"cartoon": {"color": "spectrum"}})


        except Exception as er:
            print("Failed to open protein")
        # try:
        #    if self.ligand and self.resn is not None:
        #        self.vobj.setStyle(
        #            {"resn": f"{self.ligand}"},
        #            {"stick": {"colorscheme": self.lig_color}},
        #        )
        #
        #        #    elif self.ligand is not None and self.resn is None:
        #        #        mol2 = open(self.ligand, 'r').read()
        #        #        lig_dir, lig_name, lig_file_format = give_id(self.ligand)
        #        #        if lig_file_format == "sdf":
        #        #            file_format = "sdf"
        #        #        self.vobj.addModel(mol2, f"{file_format}")
        #        #        self.vobj.setStyle(
        #        #            {'model': 1}, {'stick': {"colorscheme": self.lig_color}}
        #        #        )
        #
        #        #    elif self.ligand is None and self.resn is not None:
        #        #        self.vobj.setStyle(
        #        #            {"resn": f"{self.resn}"},
        #        #            {"sphere": {"colorscheme": self.lig_color}},
        #        #        )
        #        #except TypeError as er:
        #        #    print("failed to open ligand")
        #    pass
        try:
            if self.ligand is not None:
                mol2 = open(self.ligand, "r").read()
                lig_dir, lig_name, lig_file_format = give_id(self.ligand)
                if lig_file_format == "sdf":
                    file_format = "sdf"
                self.vobj.addModel(mol2, f"{file_format}")
                self.vobj.setStyle({"model": 1}, {"stick": {}})
            else:
                self.vobj.setStyle(
                    {"resn": f"{self.resn}", "clickable": True},
                    {"stick": {"colorscheme": self.lig_color}},
                )
        except Exception as er:
            pass

        self.vobj.zoomTo()
        # self.vobj.setStyle({"clickable": True})
        if self.save is True:
            prefix = "image"
            while os.path.exists(f"./images/{prefix}.png"):
                suffix += 1
                name = f"{prefix}{suffix}.png"
            self.vobj.save_fig(f"./Images/{name}", dpi=600)
            print(f"Successfully saved ./Images/{name} ")
        # return self.vobj.show()

