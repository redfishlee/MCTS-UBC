from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation, ReactionType, Reaction
from rdkit.Chem import AllChem as Chem




class C1OneReactantOneProductReactionLib:
    reaction_remove_H_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "H", 1, "C"),
        "Remove H on C",
        ReactionType.OneReactantOneProduct)]

    reaction_add_H_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, 1, "C", 0, 3),
        "Add H on C",
        ReactionType.OneReactantOneProduct)]

    reaction_add_H_on_the_O = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            1, 1, "O", 0, 1),
        "Add H on O",
        ReactionType.OneReactantOneProduct)]

    reaction_remove_H_on_O = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "H", 1, "O"),
        "Remove H on O",
        ReactionType.OneReactantOneProduct)]

    reaction_add_OH_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, 1, "C", 1, 3),
        "Add OH on C",
        ReactionType.OneReactantOneProduct)]

    reaction_add_O_double_bond_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, 1, "C", 0, 3, Chem.BondType.DOUBLE),
        "Add O on C",
        ReactionType.OneReactantOneProduct)]

    reaction_add_O_single_bond_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_add_atom_A_on_nth_atom_B_with_m_Hs_to_B_if_B_have_bonds_fewer_equal_than_l(
            8, 1, "C", 0, 3, Chem.BondType.SINGLE),
        "Add O on C",
        ReactionType.OneReactantOneProduct)]

    reaction_remove_first_O_or_OH_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            1, "O", 1, "C", True),
        "Remove 1st O/OH on C",
        ReactionType.OneReactantOneProduct)]

    reaction_remove_second_O_or_OH_on_carbon = [Reaction(
        OneReactantOneProductMetaReation.reaction_remove_kth_atom_A_on_nth_atom_B(
            2, "O", 1, "C", True),
        "Remove 2nd O/OH on C",
        ReactionType.OneReactantOneProduct)]


