#Here stores all the functionalities that redfish's project need
from collections import Counter
from rdkit.Chem import AllChem as Chem
from MEACRNG.MolEncoder.C1MolLib import C1MolLib
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation
from MEACRNG.Tools.utils import ReactionLibUtils
from MEACRNG.CRNG.ReactionLibrary.C3ReactionLib import C3OneReactionOneProductReactionLib
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import ReactionType, Reaction


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


def get_other_reactant(specie1, specie2):
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

#a way to simplify the main script
def extract_pathway(pg_instance, reactions, building_blocks, total_visits):
    pathway_dic, activated_IM, total_visits, reactions = pg_instance.find_a_pathway(reactions, building_blocks, total_visits)
    pathway_list = [pathway_dic[num + 1] for num in range(len(pathway_dic))]
    return pathway_list, activated_IM, total_visits, reactions

#get the lost reactant, to get the whole reaction, and then make the mkm input file
def get_the_lost_reactant(specie1, specie2):
	#all the reactions
    C = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(Chem.RWMol(Chem.MolFromSmarts("C")), j),
            "Add %s on carbon" % ('C'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CH = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CH, j),
            "Add %s on carbon" % ('CH'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CH2 = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CH2, j),
            "Add %s on carbon" % ('CH2'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CH3 = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CH3, j),
            "Add %s on carbon" % ('CH3'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CO = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CO, j),
            "Add %s on carbon" % ('CO'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CHO = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CHO, j),
            "Add %s on carbon" % ('CHO'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    COH = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.COH, j),
            "Add %s on carbon" % ('COH'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CH2O = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CH2O, j),
            "Add %s on carbon" % ('CH2O'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    CHOH = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CHOH, j),
            "Add %s on carbon" % ('CHOH'),
            ReactionType.OneReactantOneProduct)
        for j in [1, 2]
    ]

    reaction_add_H_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, i, "C", 0, 3),
        "Add H %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]]

    reaction_add_H_on_the_Os = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, i, "O", 0, 1),
        "Add H on %sth O" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]]

    reaction_add_OH_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, i, "C", 1, 3),
        "Add OH on %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]]

    reaction_add_O_single_bond_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, i, "C", 0, 3, Chem.BondType.SINGLE),
        "Add O on %s C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]]

    reaction_remove_H_on_Os = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            j, "H", i, "O"),
        "Remove H on %sth O" % i,
        ReactionType.OneReactantOneProduct) for i in range(1, 7)
    for j in [1,2]]

    reaction_remove_H_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            j, "H", i, "C"),
        "Remove H on %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]
    for j in [1,2,3,4]]

    reaction_remove_O_or_OH_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            j, "O", i, "C", True),
        "Remove 1st O/OH on %s C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2, 3]
    for j in [1,2,3,4]]

    H = ReactionLibUtils.merge_reactions([reaction_add_H_on_carbons, reaction_add_H_on_the_Os])
    OH = ReactionLibUtils.merge_reactions([reaction_add_OH_on_carbons])
    O = ReactionLibUtils.merge_reactions([reaction_add_O_single_bond_on_carbons])
    reaction_group = {'H': H, 'OH': OH, 'O': O, 'C': C, 'CH': CH, 'CH2': CH2, 'CH3': CH3, 'CO': CO, 'CHO': CHO, 'COH': COH, 'CH2O': CH2O, 'CHOH': CHOH}
    for reaction_string, reaction in reaction_group.items():
    	print(reaction_string)
    	flag_of_find_reaction = find_reaction(specie1, specie2, reaction)
    	if flag_of_find_reaction == True:
    		return reaction_string
    print(f"doesn't find the lost reactant, and the species are {specie1} and {specie2}")


def find_reaction(specie1, specie2, reactions):
	mol1 = Chem.RWMol(Chem.MolFromSmarts(specie1))
	mol2 = Chem.RWMol(Chem.MolFromSmarts(specie2))
	for index, reaction in enumerate(reactions):
		mol1 = Chem.RWMol(Chem.MolFromSmarts(specie1))
		mol2 = Chem.RWMol(Chem.MolFromSmarts(specie2))
		print(index+1)
		print(Chem.MolToSmiles(mol1))
		product = reaction.reaction_func(mol1)
		if type(product) != bool:
			product_specie = Chem.MolToSmiles(product)
			print(product_specie)
			print('mol1')
			if product_specie == specie2:
				return True
		product = reaction.reaction_func(mol2)
		if type(product) != bool:
			product_specie = Chem.MolToSmiles(product)
			print(product_specie)
			print('mol2')
			if product_specie == specie1:
				return True

