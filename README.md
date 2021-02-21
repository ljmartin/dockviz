# Docking curation with streamlit and py3dmol

see [the app](https://share.streamlit.io/ljmartin/dockviz/main/stApp.py)

see [intro post](https://ljmartin.github.io/sideprojects/dockviz.html) for motivation.

## blurb
This is a streamlit app to visualize docking hits for the purposes of manual curation before moving to _in vitro_ testing. You might use this to prioritize a ranked list of docking hits so that you only spend money on buying the ligands with the best chance of successfully binding the target.

In this test case, the 5 example ligands have been docked against the D4 receptor, crystal structure 5WIU bound to nemonapride. The ligands were pulled from a public dataset in [doi](https://doi.org/10.1038/s41586-019-0917-9)). Smina was used for docking, and obabel for file conversion to/from pdbqt/pdb. The docking can be done within `run_smina.ipynb` assuming you have smina installed.

## how to use
Load the app, and select whether or not to use Annotations that describe the pharmacophores (these can be helpful but sometimes messy). Write a number, from 0 - 4 inclusive, and press enter to load one of the docked ligand poses. How you judge the likelihood of a hit is up to you ;)

## how to adapt
Have a set of docked ligands named `0.pdb`, `1.pdb`, etc... in the `files` dir, up to as many as you want. You'll also need to have a receptor pdb file in the `files` dir, and change the name in `stApp.py` to the name of your receptor. Same goes for the co-crystallized ligand.

## credits
This pulls together some amazing libraries written by others:

- `py3dmol` for 3d visualization: https://3dmol.csb.pitt.edu/ , Rego & Koes: https://doi.org/10.1093/bioinformatics/btu829
- molecular properties are calculated with the `RDKit`: https://www.rdkit.org/
- docking used `smina`: https://sourceforge.net/projects/smina/
- docking pipeline was described by `Esben Jannik Bjerrum` at https://www.cheminformania.com/ligand-docking-with-smina/
- ligands are from a publicly available dataset from `Lyu et al`: https://doi.org/10.1038/s41586-019-0917-9
- structure file conversions by obabel: https://openbabel.org/docs/dev/Command-line_tools/babel.html
- and the pieces were put together by lewis martin: https://ljmartin.github.io/