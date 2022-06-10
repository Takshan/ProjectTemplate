import re
from collections import Counter
from os.path import join, splitext
from rdkit import Chem
from rich.console import Console
from rich.table import Table

from csfdock.DVisualize import *
from csfdock.MolView import *
from csfdock.utils import get, PDBParse

console = Console()
blue_console = Console(style="white on blue")


class ProjectStart(MolView, DVisualize):
    """Creates a Project in Optimizing Scores and Docking
    Args:
        *args: Description
        **kwargs: Description
    """

    def __init__(self, *args, **kwargs):
        self.PROJECT_DIR = kwargs.get("path", os.getcwd())
        super().__init__(*args, **kwargs)
        self.AA = [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ]

    def SetFolders(self, *args, **kwargs):
        actual_cwd = os.getcwd()
        self.PROJECT_DIR = kwargs.get("path", os.getcwd())
        for i in args:
            self.PROJECT_DIR = i
        if self.PROJECT_DIR == ".":
            return console.print(f"Project Base Directory: {actual_cwd}")
        if actual_cwd != self.PROJECT_DIR:
            try:
                os.chdir(self.PROJECT_DIR)
                working_dir = self.PROJECT_DIR
            except Exception as err:
                working_dir = f"{actual_cwd}/{self.PROJECT_DIR}"
                os.chdir(working_dir)
        console.print(f"Project Base Directory: {working_dir}")

    def ProjectTree(self, *args, **kwargs):
        path = os.getwd() if self.PROJECT_DIR is None else self.PROJECT_DIR
        verbose = kwargs.get("verbose",)
        self.directory_tree(path, verbose=verbose)

    def __actual_dir_name(self, path, root=None):
        """helper function for directory tree generation"""
        if root is not None:
            path = os.path.join(root, path)
        result = os.path.basename(path)
        if os.path.islink(path):
            realpath = os.readlink(path)
            result = f"{os.path.basename(path)} -> {realpath}"
        return result

    def directory_tree(self, startpath, verbose=True, depth=-1):
        """Tree view of the project directory tree"
        directory_tree(path)
        """
        supported_file_format = {"txt", "pdb",
                                 "pdbqt", "sdf", "csv", "excel", "pickle"}
        console.print(
            f"Supported File Format :{supported_file_format}", style="bold green")
        table = self.__create_table("bold magenta", "File Type", "Total Files")
        c = Counter(
            [splitext(i)[1][1:] for i in glob(join(startpath, "**"),
             recursive=True) if splitext(i)[1][1:] in supported_file_format]
        )
        console.print("============Details of files====================")
        for ext, count in c.most_common():
            table.add_row(
                f"[bold green]{str(ext)}[/bold green]", f"[red]{str(count)}[/red]"
            )
        console.print(table)
        if verbose:
            console.print("============Directory Tree====================")
            prefix = 0
            if startpath != "/":
                if startpath.endswith("/"):
                    startpath = startpath[:-1]
                prefix = len(startpath)
            for root, dirs, files in os.walk(startpath):
                level = root[prefix:].count(os.sep)
                if depth > -1 and level > depth:
                    continue
                indent = subindent = ""
                if level > 0:
                    indent = "|   " * (level - 1) + "|-- "
                subindent = "|   " * (level) + "|-- "
                print(
                    f"{indent}📂{self.__actual_dir_name(root)}/"
                )  # print dir only if symbolic link; otherwise, will be printed as root
                for d in dirs:
                    if not d.startswith("."):
                        if os.path.islink(os.path.join(root, d)):
                            print(
                                f"{subindent}📃{self.__actual_dir_name(d, root=root)}")
                for f in files:
                    _format = f.rsplit(".")[-1]
                    if _format in supported_file_format:
                        print(f"{subindent}📃{self.__actual_dir_name(f, root=root)}")
                    else:
                        pass

    def __receptor_contents_print(self, receptor, receptor_content):
        number_of_residues = []
        number_of_membrane_molecules = []
        number_of_water_molecule = 0
        present_ions = []
        number_of_chains = []
        number_of_ligands = []
        number_of_ligands_atoms = 0
        for index, line in enumerate(receptor_content):
            if line.startswith("ATOM"):
                if line[17:20] in self.AA or line.split()[-1] == "PROA":
                    number_of_chains.append(line[21])
                    number_of_residues.append(line[22:26])
                elif line.split()[-1] == "MEMB":
                    number_of_membrane_molecules.append(line[22:26])
                elif line.split()[-1] == "TIP3" or line[17:20] == "HOH":
                    number_of_water_molecule += 1
                elif line.split()[-1] == "HETA":
                    number_of_ligands.append(line.split()[3])
                else:
                    present_ions.append(line[17:20])
            elif line.startswith("HETATM"):
                if line[17:20] == "HOH":
                    number_of_water_molecule += 1
                elif len(line[17:20].strip()) < 3:
                    present_ions.append(line[17:20])
                elif len(line[17:20].strip()) == 3:
                    # number_of_ligands.append(line.split()[3])
                    number_of_ligands.append(line[21])
                    number_of_ligands_atoms += 1
        if not present_ions:
            max_number_of_single_ions = 0
        else:
            max_number_of_single_ions = max(
                present_ions, key=present_ions.count)
        types_of_ions = set(present_ions) if present_ions else 0
        number_of_membrane_molecules = (
            number_of_membrane_molecules[-1]
            if len(number_of_membrane_molecules) > 1 and not None
            else 0
        )
        table = self.__create_table("bold blue", "Record", "Counts")
        table.add_row("[bold green]Chains:[/]",
                      f" {len(set(number_of_chains))}")
        table.add_row("[bold green]Ligands:[/]",
                      f"{len(set(number_of_ligands))}")
        try:
            table.add_row(
                "[bold green]Number of ligand atoms :[/]",
                f"{number_of_ligands.count(max(number_of_ligands, key=number_of_ligands.count))}",
            )
        except ValueError:
            table.add_row("[bold green]Number of ligand atoms :[/]", "0")

        table.add_row("[bold green]Protein residues:[/]",
                      f"{number_of_residues[-1]}")

        try:
            table.add_row(
                "[bold green]Lipids molecules :[/]", f"{number_of_membrane_molecules}"
            )
        except ValueError:
            table.add_row("[bold green]Lipids molecules :[/]", "0")

        try:
            table.add_row(
                "[bold green]Water molecules :[/]", f" {number_of_water_molecule}"
            )
        except ValueError:
            table.add_row("[bold green]Water molecules :[/]", "0")
        try:
            table.add_row("[bold green]Ions:[/]", f"{len(present_ions)}")
            table.add_row("[bold green]Ion types :[/]", f"{types_of_ions}")
        except ValueError:
            table.add_row("[bold green]Ions:[/]", "0")
            table.add_row("[bold green]Ion types :[/]", "None")

        console.print(f"\nFor[bold red] {receptor}[/]:")
        console.print(table)

    def LoadReceptor(self, *args, native=True, verbose=True, **kwargs):
        found_what = kwargs.get("key", "Receptor")

        for i in args:
            receptor = i

        info = True
        if not os.path.exists(receptor):
            try:
                _receptor = file_search(type="pdb", target=receptor)
                if len(_receptor) == 1:
                    console.print(f"{found_what}: [bold]{_receptor[0]}[/bold]")
                    info = False
                    receptor = _receptor[0]
                elif len(_receptor) > 1:
                    print(f"{found_what} {_receptor}")
                    _receptor_number = int(
                        input(f"Which {found_what} do you like to select: ")
                    )
                    receptor = _receptor[_receptor_number]
                    console.print(f"Select {found_what} : {receptor}")
                else:
                    print(f"No {found_what} found in Local directory")
                    _download = input(
                        f"Would you like to download the {receptor} from RCSB (Y/n)? "
                    )
                    confirm = ["yes", "y", "YES", "Y"]
                    if _download in confirm:
                        download_protein_msg = get(receptor)
                        # TODO path pass forward
                        console.print(download_protein_msg)
                        receptor = f"./Generated/{receptor}.pdb"
                    else:
                        console.print(
                            f"😞 {found_what}: [bold red]{receptor} [/]to process"
                            " further.."
                        )
            except Exception as er:
                print(er)

        if info and (native == False):
            console.print(f"{found_what}: [bold]{receptor}[/bold]")
            info = False
        receptor_content = open(receptor, "r")
        receptor_content = receptor_content.readlines()
        if verbose:
            self.__receptor_contents_print(receptor, receptor_content)
        # console.print(self.receptor)
        if native:
            self.receptor = receptor
        return receptor

    def __create_table(self, header_style, arg1, arg2):
        result = Table(show_header=True, header_style=header_style)
        result.add_column(arg1, style="dim", width=40)
        result.add_column(arg2, justify="right")
        return result

    def LoadLigand(self, *args, **kwargs):
        for arg in args:
            self.ligand = arg
        *_, _file_format= give_id(self.ligand)
        if _file_format.lower() == "sdf":
            inf = open(f"{self.ligand}", "rb")
            with Chem.ForwardSDMolSupplier(inf) as fsuppl:
                for mol in fsuppl:
                    if mol is None:
                        continue
                    console.print(f"{self.ligand} has {mol.GetNumAtoms()} atoms")
        elif _file_format.lower() == "pdb":
            mol = Chem.MolFromPDBFile(self.ligand, sanitize=False)
            console.print(f"{self.ligand} has {mol.GetNumAtoms()} atoms")
        else:
            return "Unknow ligand file format"

        self.ligand_export = mol
        return self.ligand_export

    def SaveComplex(self, **kwargs):
        ligand = kwargs.get("lig", None)
        receptor = kwargs.get("pro", None)
        lipid = kwargs.get("lipid", None)
        out_file = kwargs.get("out", "complex_out.pdb")
        ligand_mol = Chem.MolToPDBBlock(self.ligand_export, flavor=32)
        # print(out_file)
        *_, structure_format = give_id(self.receptor)
        lipid = lipid if lipid is not None else self.lipid
        receptor = receptor if receptor is not None else self.receptor
        ligand = ligand if ligand is not None else self.ligand

        try:
            if structure_format.lower() == "sdf":
                receptor_mol = Chem.MolFromMolFile(receptor, removeHs=False)
            elif structure_format.lower() == "pdb":
                receptor_mol, *_ = PDBParse(receptor)

        except ValueError as er:
            console.print(f"{self.receptor} Error : {er}")
            blue_console.print("Not able to parse..")

        write_list = [receptor_mol, ligand_mol]
        try:
            _prot, lipid_mol, _tip, _lig = PDBParse(lipid)
            write_list.append(lipid_mol)

        except Exception as er:
            blue_console.print(f"Cannot able to parse {lipid}.\n Error: {er}")
            blue_console.print(f"Not writing Lipid in the {out_file}")

        self.__write_to_file(write_list, out_file, check=True)

        # self.__write_to_file(lipid_mol, out_file)
        # self.__write_to_file(Chem.MolToPDBBlock(self.ligand_export, flavor=32), out_file)

    def __write_to_file(self, content, filename, check=False):
        if check:
            count = 0
            _file, _format = filename.rsplit(".")
            while os.path.exists(filename):
                count += 1
                filename = f"{_file}_{count}.{_format}"
        with open(filename, "a+") as f:
            if isinstance(content, list):
                for _ in content:
                    for line in _:
                        print(line, end="", file=f)
            else:
                for line in content:
                    print(line, end="", file=f)
        console.print(f"{filename} Saved Successfully!", style="bold green")

    def __call__(self):
        raise TypeError(
            'Project must be accessed through "instance=ProjectStart()".')

    def __str__(self):
        ligand = self.ligand if self.ligand != None else "Not implemented"
        receptor = self.receptor if self.receptor != None else "Not implemented"
        lipid = self.lipid if self.lipid != None else "Not implemented"
        return (
            f"Project : \t\nProtein: {receptor}\t\nligand :{ligand}\t\nLipid : {lipid} "
        )
