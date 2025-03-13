#Here stores all the functionalities that redfish's project need
from collections import Counter
from rdkit.Chem import AllChem as Chem
from MEACRNG.MolEncoder.C1MolLib import C1MolLib
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation


#this is the main script, using to detect the C adding reaction, so we can get that C containing specie, and then back generate it
def detect_reactions(pathway):
	mols_list = [ Chem.MolFromSmiles(specie) for specie in pathway ]
	adding_carbon_containing_nodes = []
	for mol_num in range(len(mols_list)-1):
		current_mol = mols_list[mol_num]
		next_mol = mols_list[mol_num + 1]
		current_C_count = Counter(pathway[mol_num])['C']
		next_C_count = Counter(pathway[mol_num + 1])["C"]
		if next_C_count - current_C_count == 1:
			#need to get the C containing IM, which requires lots of reaction examination
			adding_carbon_containing_nodes.append(get_other_reactant(pathway[mol_num], pathway[mol_num + 1]))
	assert len(adding_carbon_containing_nodes) == 2, f"did't find two carbon adding nodes, pathways is {pathway}"
	return adding_carbon_containing_nodes


def get_other_reactant(mol1, mol2):
	mol1 = Chem.MolFromSmarts(mol1)
	mol2 = Chem.MolFromSmarts(mol2)
	carbon_containing_building_blocks = [
            Chem.RWMol(Chem.MolFromSmarts("C")),
            C1MolLib.CH,
            C1MolLib.CH2,
            C1MolLib.CH3,
            C1MolLib.CO,
            C1MolLib.CHO,
            C1MolLib.COH,
            C1MolLib.CH2O,
            C1MolLib.CHOH
            ]
	#try mol1
	for mol in carbon_containing_building_blocks:
		reaction1 = OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(mol, 1)
		reaction2 = OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(mol, 2)
		next_mol_1 = reaction1(mol1)
		next_mol_2 = reaction2(mol1)
		next_mol = [Chem.MolToSmiles(next_mol_2), Chem.MolToSmiles(next_mol_1)]
		if Chem.MolToSmiles(mol2) in next_mol:
			print(f"found the reaction!")
			print(f'it is {Chem.MolToSmiles(mol)}')
			return mol
	for mol in carbon_containing_building_blocks:
		reaction1 = OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(mol, 1)
		reaction2 = OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(mol, 2)
		next_mol_1 = reaction1(mol2)
		next_mol_2 = reaction2(mol2)
		next_mol = [Chem.MolToSmiles(next_mol_2), Chem.MolToSmiles(next_mol_1)]
		if Chem.MolToSmiles(mol1) in next_mol:
			print(f"found the reaction!")
			print(f'it is {Chem.MolToSmiles(mol)}')
			return mol
	print(f"didn't find the C adding block, and the current molecules are {Chem.MolToSmiles(mol1)}, {Chem.MolToSmiles(mol2)}")



