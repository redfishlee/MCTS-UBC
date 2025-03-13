from rdkit.Chem import AllChem as Chem
from MEACRNG.CRNG.ReactionNetworkGenerator import SpeciesNode, PathwayGenerator #the most important one
from MEACRNG.CRNG.ReactionLibrary.C1ReactionLib import C1OneReactantOneProductReactionLib
from MEACRNG.CRNG.ReactionLibrary.C2ReactionLib import C2OneReactionOneProductReactionLib
from MEACRNG.CRNG.ReactionLibrary.C3ReactionLib import C3OneReactionOneProductReactionLib #make the reaction
from MEACRNG.CRNG.ValidMolCheckFuncLibrary.ValidMolCheckFuncLib import CarbonValidMolCheckFunc #don't care
from MEACRNG.Tools.utils import ReactionLibUtils
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation, ReactionType, Reaction
from MEACRNG.Tools.Smiles2formula import expression_transform #in the final part, you want to make it better looking
from MEACRNG.Tools.redfish_functionalities import detect_reactions

from MEACRNG.MolEncoder.C1MolLib import C1MolLib #to make the carbon_containing_building_blocks

import json

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

#reactants and results in Smiles，
CH3CH2CH3 = SpeciesNode(Chem.MolFromSmarts('[H]C([H])([H])C([H])([H])C([H])([H])[H]'))# 'AddHs' will fill H，eg. C+4H==CH4
CO = SpeciesNode(Chem.MolFromSmiles('CO'))
#CH3 = SpeciesNode(Chem.MolFromSmiles('[H]C([H])[H]',sanitize=False))# sanitize=False，input smile which keeps H directly.

# set reactant and product
PG = PathwayGenerator.Generator(
    reactant=CO,
    product=CH3CH2CH3,
    valid_mol_check_func_list=check_func)

#some configuration
activated_species = [] #store all the intermediates that have been visited
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
            ] #this is all rdkit mol now


#get the main pathway
all_pathways = {}
total_visits = 0
while True:
    pathways_single_iteration = []
    total_visits += 1
    pathway_dic, activated_IM, total_visits, reactions = PG.find_a_pathway(reactions, carbon_containing_building_blocks, total_visits)
    pathway_list = [ pathway_dic[num+1] for num in range(len(pathway_dic)) ]
    main_pathway = PG.simplify_pathway(pathway_list)
    pathways_single_iteration.append(main_pathway)

    #identify the carbon containing building blocks
    adding_carbon_containing_nodes = detect_reactions(main_pathway)
    specie_smiles = [Chem.MolToSmiles(node) for node in adding_carbon_containing_nodes]

    #generating pathway for those two C containing species
    c_building_nodes = [SpeciesNode(Chem.MolFromSmarts(specie_smile)) for specie_smile in specie_smiles]
    for c_building_node in c_building_nodes:
        PG_c_building = PathwayGenerator.Generator(
        reactant=CO,
        product=c_building_node,
        valid_mol_check_func_list=check_func)

        pathway_dic, activated_IM, total_visits, reactions = PG_c_building.find_a_pathway(reactions, carbon_containing_building_blocks, total_visits)
        pathway_list = [ pathway_dic[num+1] for num in range(len(pathway_dic)) ]

        c_building_block_pathway = PG.simplify_pathway(pathway_list)
        pathways_single_iteration.append(c_building_block_pathway)

    assert isinstance(pathways_single_iteration, list), "this pathways is not a list"
    assert len(pathways_single_iteration) == 3, f"this pathways doesn't have 3 pathways, current pathways is {pathways_single_iteration}"
    all_pathways[total_visits] = pathways_single_iteration
    if total_visits >= 100:
        break

with open("All_pathways.json", "w") as f:
    json.dump(all_pathways, f, indent=4)


