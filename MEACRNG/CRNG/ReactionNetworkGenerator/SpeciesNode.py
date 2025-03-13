from rdkit.Chem import AllChem as Chem
from MEACRNG.Tools.utils import ReactFlags
from MEACRNG.Tools import MolViz as mv
import copy
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import Reaction,ReactionType
import math

class SpeciesNode(object):
    '''
    play as a node in reaction network
    '''

    def __init__(self,rdkit_mol):
        self.mol = rdkit_mol
        self.mol = Chem.RWMol(self.mol)# read/wirte
        self.smiles = Chem.MolToSmiles(self.mol)
        # if the species did the reaction, then set to True, orop = one reactant one product
        self.orop_reaction_to_bReacted_dict = {}
        self.orop_reaction_list = []
        #redfish
        self.activated_flag = False #important signal, one usage is to tell whether stop or keep search which require C-containing intermediate in pathway searching process
        #redfish
        #following code will add attributes to each speciesnode for further grading in the pathway searching process
        self.visited_times = 0
        self.score = 0
        self.reaction_energy = 0
        self.barrier_energy = 0
        self.grade = 0

    #redfish, give node grade, so that the search could begin
    def node_grading(self, total_visits):
        exploration_param = 1.7
        if self.visited_times == 0:
            self.grade == float('inf')
        else:
            self.grade = self.score + exploration_param * math.sqrt(math.log(total_visits + 1) / self.visited_times) 

    #redfish, set all flags to false
    def sign_one_reactant_one_product_reactions(self, reaction_list):#将所有的反应的flag变成false，这样相当于重启所有反应，所有反应都没有发生过
        self.reaction_list = reaction_list
        for i in reaction_list:
            assert isinstance(i,Reaction)
            assert i.reaction_type == ReactionType.OneReactantOneProduct# 一反应物一产物类型的反应
            self.orop_reaction_to_bReacted_dict[i] = False

    #redfish, to generate all leaves from one leave, just one step, won't go down forever
    def react_one_reactant_one_product(self, new_oto_reactions_to_sign_for_product=None):#对一个中间体去发生所有可能的反应，然后生成可能的产物
        results = []
        used_reactions_names = []
        for reaction in self.orop_reaction_to_bReacted_dict:
            if self.orop_reaction_to_bReacted_dict[reaction] == True:
                continue
            self.orop_reaction_to_bReacted_dict[reaction] = True
            mol = copy.deepcopy(self.mol)
            result = reaction.reaction_func(mol)

            if result:
                assert isinstance(result, Chem.RWMol)
                result = SpeciesNode(copy.deepcopy(result))
                if new_oto_reactions_to_sign_for_product is not None:#set all the flag false for all newly generated mol
                    assert isinstance(new_oto_reactions_to_sign_for_product,list)
                    result.sign_one_reactant_one_product_reactions(new_oto_reactions_to_sign_for_product)
                results.append(result)
                used_reactions_names.append(reaction.name)

        if len(results) == 0:
            return ReactFlags.NoProduct,None
        return results,used_reactions_names



"""

    def show(self):
        mv.show_mol(self.mol)

    def sign_one_reactant_one_product_reactions(self, reaction_list):#将所有的反应的flag变成false，这样相当于重启所有反应，所有反应都没有发生过
        self.reaction_list = reaction_list
        for i in reaction_list:
            assert isinstance(i,Reaction)
            assert i.reaction_type == ReactionType.OneReactantOneProduct# 一反应物一产物类型的反应
            self.orop_reaction_to_bReacted_dict[i] = False

    def react_one_reactant_one_product(self, new_oto_reactions_to_sign_for_product=None):#对一个中间体去发生所有可能的反应，然后生成可能的产物
        results = []
        used_reactions_names = []
        for reaction in self.orop_reaction_to_bReacted_dict:
            if self.orop_reaction_to_bReacted_dict[reaction] == True:
                continue
            self.orop_reaction_to_bReacted_dict[reaction] = True
            mol = copy.deepcopy(self.mol)
            result = reaction.reaction_func(mol)

            if result:
                assert isinstance(result, Chem.RWMol)
                result = SpeciesNode(copy.deepcopy(result))
                if new_oto_reactions_to_sign_for_product is not None:
                    assert isinstance(new_oto_reactions_to_sign_for_product,list)
                    result.sign_one_reactant_one_product_reactions(new_oto_reactions_to_sign_for_product)
                results.append(result)
                used_reactions_names.append(reaction.name)

        if len(results) == 0:
            return ReactFlags.NoProduct,None
        return results,used_reactions_names

    def __str__(self):
        return Chem.MolToSmiles(self.mol)




if __name__ == '__main__':
    # TODO: move following code into the documentation of the method!
    # mol = Chem.MolFromSmiles('OCC')
    # mol = Chem.AddHs(mol)
    # s = SpeciesNode(mol)
    # #print(s.get_index_of_the_nth_atom_of_element_with_total_m_bonds(1,"O",2))
    # # status = s.remove_the_bond_with_certain_type_and_connect_to_certain_atom_in_base_atom(
    # #     base_atom_idx=0,
    # #     bond_type=Chem.BondType.SINGLE,
    # #     bonded_atom="H")
    # #print(status)
    # s.bonding_an_new_atom_to_the_atom_with_certain_index_with_certain_bond_type(
    #     atom_index=2,
    #     atomic_number_of_the_atom_to_add=6,
    #     bond_type=Chem.BondType.SINGLE,
    #     num_of_H_add_to_the_new_atom=2
    # )
    # s.show()
    # print(s.mol.GetNumBonds())
    #s.mol = Chem.AddHs(s.mol)
    #s.show()



    # CH4 test generate
    # mol = Chem.MolFromSmiles("C")
    # mol = Chem.AddHs(mol)
    # s = SpeciesNode(mol)
    # carbon_index = s.get_index_of_the_nth_atom_of_element_with_total_m_bonds(1,"C",-1)
    # new_s = []
    # # move Hs
    # for _ in range(3):
    #     s.remove_the_bond_with_certain_type_and_connect_to_certain_atom_in_base_atom(
    #         bond_type=Chem.BondType.SINGLE,
    #         bonded_atom="H",
    #         base_atom_idx=carbon_index)
    #     s.remove_atoms_without_bonds()
    #     new_s.append(copy.deepcopy(s))
    # new_new_s = []
    # for i in new_s:
    #     i.bonding_an_new_atom_to_the_atom_with_certain_index_with_certain_bond_type(
    #         atom_index=carbon_index,
    #         atomic_number_of_the_atom_to_add=8,
    #         bond_type=Chem.BondType.SINGLE,
    #         num_of_H_add_to_the_new_atom=1 # OH
    #     )
    #     print(i.mol.GetNumAtoms())
    #     print(Chem.MolToSmiles(i.mol))
    #     new_new_s.append(copy.deepcopy(i))
    # for i in new_new_s:
    #     i.show()
    def remove_H_reaction():
        pass


"""





