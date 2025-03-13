import copy
import os
import random
from itertools import product# 迭代功能的内置库
from MEACRNG.Tools.Smiles2formula import expression_transform# 编码转换的工具
from MEACRNG.MolEncoder.C1MolLib import C1MolLib
from rdkit.Chem import AllChem as Chem
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import ReactionType, Reaction
from MEACRNG.CRNG.ReactionLibrary.MetaReactions import OneReactantOneProductMetaReation

#for testing
import time


from MEACRNG.CRNG.ReactionNetworkGenerator.SpeciesNode import SpeciesNode
from MEACRNG.Tools.utils import ReactFlags, color_string

class ReactionStep(object):

    def __init__(self, reactant_list, product_list, reaction_name):
        assert isinstance(reactant_list, list)# 判断输入变量类型
        assert isinstance(product_list, list)
        assert isinstance(reaction_name, str)
        self.reactant_list = reactant_list
        self.product_list = product_list
        self.reaction_name = reaction_name
        self.color_string= False

    def __str__(self):
        rs = " + ".join([(Chem.MolToSmiles(i.mol)) for i in self.reactant_list])
        ps = " + ".join([(Chem.MolToSmiles(i.mol)) for i in self.product_list])

        # 编码转换，暂时不用，先依然print难以看懂的形式
        # rs = expression_transform('C',rs)
        # ps = expression_transform('C',ps)

        if self.color_string == True:

            return color_string(rs, "RED") + \
               " > " + color_string(self.reaction_name, "BLUE") + \
               " > " + color_string(ps, "GREEN")
        else:
            return rs + \
                   " > " + self.reaction_name + \
                   " > " + ps

    @property
    def hash_regardless_direction(self):
        '''
        返回反应物和产物的hash值
        这个值会随着内存位置变化
        :return:
        '''
        rs = " + ".join([str(i) for i in self.reactant_list])
        ps = " + ".join([str(i) for i in self.product_list])
        return hash(rs) + hash(ps)



class Generator(object):

    def __init__(self, reactant, product=None, valid_mol_check_func_list=None):
        self.reactant = reactant# 传入初始物种
        self.product = product  # 传入产物
        self.valid_mol_check_func_list = valid_mol_check_func_list# 传入check_function
        if product is not None:# 暂时产物不会影响结果
            print("目前final results用于判断是否在网络内，网络会生成到没有新物种为止")
        assert isinstance(reactant, SpeciesNode)
        assert isinstance(product, SpeciesNode)
        self.steps = []


    def find_a_pathway(self, reactions, carbon_containing_building_blocks, total_visits):
        #configuration
        activated_IM = [] #IM: intermediate
        activated_IM.append(self.reactant) #attention, the list stores nodes
        current_node = self.reactant
        iteration_num = 0
        while True:#need to get a reliable pathway
            current_node = self.reactant
            pathway = {}
            pathway_found = False
            step_number = 0
            iteration_num += 1 #add iteration number if the pathway is not found
            current_node.sign_one_reactant_one_product_reactions(reactions)#set the flag to false, so we will consider all the reactions
            while True:#break when pathway is too long or succeed
                step_number += 1
                pathway[step_number] = Chem.MolToSmiles(current_node.mol)
                
                #get next node contain all the next level leaves from the former leave, and you don't need to care about reaction names
                next_nodes, reaction_names = current_node.react_one_reactant_one_product()
                
                #check is the new nodes are valid
                next_nodes = [node for node in next_nodes if self.check_if_mol_valid(node)]
                
                #grade all the next level nodes
                smiles_of_IM = {node.smiles: node for node in activated_IM}

                for node in next_nodes:
                    node_smiles = Chem.MolToSmiles(node.mol)
                    if node_smiles in smiles_of_IM:
                            node = smiles_of_IM[node_smiles]
                    node.node_grading(total_visits)


                #randomly pick one node with highest grade
                max_grade = max(node.grade for node in next_nodes)
                top_nodes = [node for node in next_nodes if node.grade == max_grade]
                random_top_node = random.choice(top_nodes)
                
                #upgrade items in activated_IM if specienode are already in it. if so, change the activated_IM instead. if not, change node, and then add it to activated_IM
                smiles_of_IM = {node.smiles: node for node in activated_IM}
                if random_top_node.smiles in smiles_of_IM.keys():
                    smiles_of_IM[random_top_node.smiles].visited_times += 1
                else:   
                    random_top_node.visited_times += 1
                    random_top_node.activated_flag = True
                    activated_IM.append(random_top_node)


                #update the status
                current_step = ReactionStep([copy.deepcopy(current_node)], [copy.deepcopy(random_top_node)], 'donno, and doncare')
                current_node = copy.deepcopy(random_top_node)
            
                #update the reactions(if we have IM of CH, we should be able able to have a CH adding reaction)
                if carbon_containing_building_blocks != None:
                    for rdkit_mol in carbon_containing_building_blocks:
                        if current_node.smiles == Chem.MolToSmiles(rdkit_mol):
                            new_reaction = Reaction(
                            OneReactantOneProductMetaReation.reaction_add_mol_on_nth_carbon_and_make_CC_single_bond(current_node.mol, 1),
                            "Add %s on carbon" % (Chem.MolToSmiles(current_node.mol)),
                            ReactionType.OneReactantOneProduct)
                            carbon_containing_building_blocks.remove(rdkit_mol)
                            reactions.append(new_reaction)

                #set actionable reactions to false for the select node
                current_node.sign_one_reactant_one_product_reactions(reactions)

                #check: 1.is pathway too long? 2.have we reached product yet?
                if len(pathway) > 50:
                    iteration_num += 1
                    break
                if Chem.MolToSmiles(current_node.mol) == Chem.MolToSmiles(self.product.mol):
                    iteration_num += 1
                    pathway_found = True
                    step_number += 1
                    pathway[step_number] = Chem.MolToSmiles(current_node.mol)
                    #print("congratulations! A reliable pathway has been found!")
                    return pathway, activated_IM, total_visits, reactions


    def check_if_mol_valid(self, species_node):
        if self.valid_mol_check_func_list is None:
            return True
        for func in self.valid_mol_check_func_list:
            if func(species_node.mol) == False:
                return False
        return True

    #simplify the pathway
    def simplify_pathway(self, full_pathway):
        full_pathway_copy = full_pathway[:]
        while True:
            found_duplicate = False  # Flag to control breaking out of both loops
            for main_number in range(len(full_pathway_copy)):
                for find_number in range(len(full_pathway_copy)):
                    if full_pathway_copy[main_number] == full_pathway_copy[find_number] and find_number != main_number:
                        #print(full_pathway_copy[main_number:find_number])
                        del full_pathway_copy[main_number:find_number]
                        found_duplicate = True  # Set flag to True to indicate a break
                        break  # Break the inner loop
                if found_duplicate:
                    break  # Break the outer loop if a duplicate was found
            if not found_duplicate:
                # No duplicates found in this iteration, exit the while loop
                break
                
        return full_pathway_copy

