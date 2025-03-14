#some experiments for funcionality that i want to implement in the main scripts
#same setting in the main Pathway_of_...
from rdkit.Chem import AllChem as Chem
from MEACRNG.CRNG.ReactionNetworkGenerator import SpeciesNode, PathwayGenerator #the most important one
from MEACRNG.CRNG.ReactionLibrary.C1ReactionLib import C1OneReactantOneProductReactionLib
from MEACRNG.CRNG.ReactionLibrary.C2ReactionLib import C2OneReactionOneProductReactionLib
from MEACRNG.CRNG.ReactionLibrary.C3ReactionLib import C3OneReactionOneProductReactionLib #make the reaction
from MEACRNG.CRNG.ValidMolCheckFuncLibrary.ValidMolCheckFuncLib import CarbonValidMolCheckFunc #don't care
from MEACRNG.Tools.utils import ReactionLibUtils
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation, ReactionType, Reaction
from MEACRNG.Tools.Smiles2formula import expression_transform #in the final part, you want to make it better looking

from MEACRNG.MolEncoder.C1MolLib import C1MolLib #to make the carbon_containing_building_blocks
import json
import ast
from collections import Counter

# allowed steps,for multi carbon products，use functions in MetaReation or write a new lib
reactions = C3OneReactionOneProductReactionLib.reaction_add_carbon_species_on_jth_carbon

# boundary condition check
check_func = [
    CarbonValidMolCheckFunc.no_more_than_four_bond_in_carbon_check,
    CarbonValidMolCheckFunc.no_more_than_M_bond_in_atom_N_check(n='O',m=2),
    CarbonValidMolCheckFunc.valid_C_O_bond_structure_check_when_only_use_single_CO_bond_simple,
    CarbonValidMolCheckFunc.no_more_than_three_carbon_check,
    CarbonValidMolCheckFunc.num_of_C_and_O_should_smaller_equal_than_4,
    CarbonValidMolCheckFunc.num_of_CHO_should_bigger_equal_than_1
]
#[H]C([H])C([H])([H])C([H])([H])O and [H]OC([H])([H])C([H])([H])C([H])[H]

smiles = "[H]C([H])C([H])([H])C([H])([H])O"
from MEACRNG.Tools.Smiles2formula import expression_transform# 编码转换的工具
formula = expression_transform(smiles)
print(formula)




from MEACRNG.Tools.redfish_functionalities import get_the_lost_reactant

with open("network_with_energy.json", "r") as file:
    network_with_energy = json.load(file)
reactants_set_1 = set()
reactants_set_2 = set()
reactants_list_1 = []
reactants_list_2 = []
for reactants_couple in network_with_energy.keys():
    reactants_couple = ast.literal_eval(reactants_couple)
    #print(reactants_couple[0])
    specie1 = expression_transform(reactants_couple[0])
    specie2 = expression_transform(reactants_couple[1])
    reactants_set_1.add(specie1)
    reactants_set_1.add(specie2)
print(len(reactants_set_1))





"""print(len(reactants_set_1))
print(len(reactants_set_2))
print(len(reactants_list_1))
print(reactants_set_1)

duplicates = {item for item, count in Counter(reactants_list_2).items() if count > 1}
print('hhh')
print(len(reactants_list_2))
print(len(duplicates))"""

"""final_list = list(reactants_set_1)
smooth_list = [specie.replace('[','').replace(']','').replace('(','').replace(')','') for specie in final_list]
print(smooth_list)
duplicates = {item for item, count in Counter(smooth_list).items() if count > 1}
print(len(duplicates))
print(duplicates)"""


"""specie1 = "[H]CC([H])([H])C([H])([H])O"
specie2 = "[H]C([H])C([H])([H])C([H])([H])O"
lost_reactant = get_the_lost_reactant(specie1, specie2)
print('hhh')
print(lost_reactant)
mol1 = Chem.RWMol(Chem.MolFromSmarts(specie1))
mol2 = Chem.RWMol(Chem.MolFromSmarts(specie2))

specie1_copy = specie1.replace('[','').replace(']','').replace('(','').replace(')','')
specie2_copy = specie2.replace('[','').replace(']','').replace('(','').replace(')','')  

longer_one = max(specie1_copy, specie2_copy, key=len)
shorter_one = min(specie1_copy, specie2_copy, key=len)
side1 = longer_one
side2 = [shorter_one, lost_reactant]
print(side1, side2)


#reaction1 = OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(1, 2, "C", 0, 3)
reaction1 = OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(1, 1, "O", 0, 1)



print(Chem.MolToSmiles(mol1))
print(type(reaction1))
hhh = reaction1(mol1)
print(type(hhh))
print(Chem.MolToSmiles(hhh))"""



"""#reactants and results in Smiles，
CH3CH2CH3 = SpeciesNode(Chem.MolFromSmarts('[H]C([H])([H])C([H])([H])C([H])([H])[H]'))# 'AddHs' will fill H，eg. C+4H==CH4
CO = SpeciesNode(Chem.MolFromSmiles('CO'))
#CH3 = SpeciesNode(Chem.MolFromSmiles('[H]C([H])[H]',sanitize=False))# sanitize=False，input smile which keeps H directly.

#from this place it's the experiment part
import copy

CH3CH2CH3 = Chem.Mol(CH3CH2CH3.mol)
print(type(CH3CH2CH3))
check = Chem.MolFromSmarts("[H]C[H](C([H])([H])[H])C([H])([H])[H]")
print(type(check))

check = Chem.MolToSmiles(check)
print(check)
CH3CH2CH3 = Chem.MolToSmiles(CH3CH2CH3)
print(CH3CH2CH3)
print(type(check))
print(type(CH3CH2CH3))

mol2 = Chem.CanonSmiles(CH3CH2CH3)
print(mol2)
mol1 = Chem.MolFromSmarts("[H]C[H](C([H])([H])[H])C([H])([H])[H]")
print(mol1)
print(mol2)
if mol1 == mol2:
    print("good")


reaction = OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(1, 1, "C", 0, 3)
molecule = Chem.MolFromSmarts("C-O-[H]")
molecule = Chem.RWMol(molecule)
print(type(reaction))
result = reaction(molecule)
print(type(result))
string = Chem.MolToSmiles(result)
print(string)

#test
reaction = OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(2, "H", 1, "O")
molecule = Chem.MolFromSmarts("C-O-[H]")
molecule = Chem.RWMol(molecule)
print(type(molecule))
print(type(reaction))

r1 = OneReactantOneProductMetaReation.reaction_remove_the_kth_bond_with_certain_type_and_connect_to_certain_atom_in_base_atom(
                k=2,
                base_atom_idx=1,
                bond_type=None,
                bonded_atom="H",
                brake_all_bonds_in_bonded_atom=False
            )

molecule1 = r1(molecule)
print(type(molecule1))
string = Chem.MolToSmiles(molecule1)
print(string)

r2 = OneReactantOneProductMetaReation.reaction_remove_atoms_without_bonds(
                remove_ignore="C"
            )

molecule2 = r2(molecule1)

print(type(molecule2))
string = Chem.MolToSmiles(molecule2)
print(string)


print('*'*100)
string = Chem.MolToSmiles(molecule)
print(string)
result = reaction(molecule)
print(type(result))
string = Chem.MolToSmiles(result)
print(string)


reaction = OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            3, "O", 2, "C", True)
molecule = Chem.MolFromSmarts("[H]CC([H])(O[H])C([H])[H]")
if molecule is None:
    print("damn")
#molecule = Chem.AddHs(molecule)
molecule = Chem.RWMol(molecule)
result = reaction(molecule)
string = Chem.MolToSmiles(result)
print(string)

#[H]CC([H])(O[H])C([H])[H]"""
