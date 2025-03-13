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

# allowed steps,for multi carbon products，use functions in MetaReation or write a new lib
reactions = ReactionLibUtils.merge_reactions([
    #C3OneReactionOneProductReactionLib.reaction_add_carbon_species_on_jth_carbon,
    C3OneReactionOneProductReactionLib.reaction_add_H_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_OH_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_O_single_bond_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_H_on_the_Os,
    C3OneReactionOneProductReactionLib.reaction_remove_H_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_remove_H_on_Os,
    C3OneReactionOneProductReactionLib.reaction_remove_O_or_OH_on_carbons,
])

# boundary condition check
check_func = [
    CarbonValidMolCheckFunc.no_more_than_four_bond_in_carbon_check,
    CarbonValidMolCheckFunc.no_more_than_M_bond_in_atom_N_check(n='O',m=2),
    CarbonValidMolCheckFunc.valid_C_O_bond_structure_check_when_only_use_single_CO_bond_simple,
    CarbonValidMolCheckFunc.no_more_than_three_carbon_check,
    CarbonValidMolCheckFunc.num_of_C_and_O_should_smaller_equal_than_4,
    CarbonValidMolCheckFunc.num_of_CHO_should_bigger_equal_than_1
    
]


mol1 = Chem.MolFromSmarts("[H]CC([H])O[H]")
mol2 = Chem.MolFromSmarts("[H]CC([H])(O[H])C([H])[H]")
reaction = OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(C1MolLib.CH2, 1)
next_mol_1 = reaction(mol1)
if Chem.MolToSmiles(next_mol_1) == Chem.MolToSmiles(mol2):
    print('good')
print(Chem.MolToSmiles(next_mol_1))
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
