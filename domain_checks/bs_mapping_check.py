#!/bin/python

import glob
import csv
import re
import pandas as pd
import itertools
import numpy as np
from tqdm.notebook import tqdm
import csv
import matplotlib.pyplot as plt 
import seaborn as sns 
plot_kwds = {'alpha' : 0.5, 's' : 80, 'linewidths':0}
from collections import defaultdict

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')
from rdkit.Chem import Draw, rdFMCS

from rdkit.Chem import rdRGroupDecomposition, AllChem, rdmolops
from rdkit.Chem import rdqueries
from rdkit.Chem import rdDepictor, rdmolfiles
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import IPythonConsole

quints_infos = pd.read_csv("output/quints_infos.csv", names=["set", "pertname", "pertsmarts", "num_ha", "sem"])
quints_infos.sort_values(by="sem")

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def plotInfo(ligname1, ligname2):
    """Plots pert mols with atom indices and atom mappings with and without sanitisation."""
    
    
    mol1 = BSS.IO.readPDB("./quintup_ligands/ligand_files/{}.pdb".format(ligname1))[0]
    mol2 = BSS.IO.readPDB("./quintup_ligands/ligand_files/{}.pdb".format(ligname2))[0]

    mol1_r = Chem.SDMolSupplier("./quintup_ligands/sdffiles/{}.sdf".format(ligname1))[0]
    mol2_r = Chem.SDMolSupplier("./quintup_ligands/sdffiles/{}.sdf".format(ligname2))[0]

    mol1_r = rdmolops.AddHs(mol1_r)
    mol2_r = rdmolops.AddHs(mol2_r)
    
    AllChem.Compute2DCoords(mol1_r)
    AllChem.Compute2DCoords(mol2_r)
    


    mapping = BSS.Align.matchAtoms(mol1, mol2)
    mapping_s = BSS.Align.matchAtoms(mol1, mol2, sanitize=True)
    
    if mapping != mapping_s:
    	return True
    else:
    	return False


#plot_perts(quints_infos.sort_values(by="sem", ascending=False).head(200))
handled = []
with open("to_do_sims/mapping_mismatches", "r") as writefile:
	writer = csv.writer(writefile)

	for pert in tqdm(quints_infos["pertname"].values):
	    
	    if not pert in handled:

	    	mismatch = plotInfo(pert.split("~")[0], pert.split("~")[1])

	    	if mismatch:
	    		
	    		# write pert to file.
	    		writer.writerow(pert)

	    		
	    
	    # make sure we don't do the same pert again (or the inverse of it).
	    handled.append(pert)
	    handled.append(pert.split("~")[1]+"~"+pert.split("~")[0])