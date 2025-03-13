from MEACRNG.CRNG.ReactionLibrary.MetaReactions import ReactionType, Reaction
from rdkit.Chem import AllChem as Chem
from MEACRNG.MolEncoder.C1MolLib import C1MolLib
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation


# C2的Lib是基于C1建立的，而并非从0开始
class C2OneReactionOneProductReactionLib:



    reaction_add_carbon_species_on_1st_carbon = [
        Reaction(
            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_with_O_and_make_CC_single_bond(i, 1),
            "Add %s on carbon" % (Chem.MolToSmiles(i)),
            ReactionType.OneReactantOneProduct)
        for i in [
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
    ]

    reaction_add_H_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, i, "C", 0, 3),
        "Add H %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_add_H_on_the_Os = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, i, "O", 0, 1),
        "Add H on %sth O" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_add_OH_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, i, "C", 1, 3),
        "Add OH on %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_add_O_single_bond_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, i, "C", 0, 3, Chem.BondType.SINGLE),
        "Add O on %s C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_remove_H_on_Os = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "H", i, "O"),
        "Remove H on %sth O" % i,
        ReactionType.OneReactantOneProduct) for i in range(1, 5)]

    reaction_remove_H_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "H", i, "C"),
        "Remove H on %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_remove_first_O_or_OH_on_carbons = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "O", i, "C", True),
        "Remove 1st O/OH on %s C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]

    reaction_remove_second_O_or_OH_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            2, "O", i, "C", True),
        "Remove 2nd O/OH on %sth C" % i,
        ReactionType.OneReactantOneProduct) for i in [1, 2]]
