import copy
import os
import random
from itertools import product# 迭代功能的内置库
from MEACRNG.Tools.Smiles2formula import expression_transform# 编码转换的工具

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
        rs = " + ".join([str(i) for i in self.reactant_list])
        ps = " + ".join([str(i) for i in self.product_list])

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



class NetworkGenerator(object):

    def __init__(self, initial_reactants, final_results=None, valid_mol_check_func_list=None):
        self.species_pool = initial_reactants# 传入初始物种
        self.product_pool = final_results  # 传入产物
        self.species_name_showed = [str(i) for i in self.species_pool]# 获得初始物种的字符串
        self.product_name_showed = [str(i) for i in self.product_pool]  # 获得初始物种的字符串
        self.valid_mol_check_func_list = valid_mol_check_func_list# 传入check_function
        self.target_pool = final_results# 传入目标产物
        if final_results is not None:# 暂时产物不会影响结果
            print("目前final results用于判断是否在网络内，网络会生成到没有新物种为止")
        for i in initial_reactants:# 输入的初始反应物必须是SpeciesNode类
            assert isinstance(i, SpeciesNode)
        if final_results is not None:
            for i in final_results:
                assert isinstance(i, SpeciesNode)
        self.steps = []#i think this is used to generate readable reaction information

    def get_all_involved_mols_and_mol_name_to_mol_dict(self):

        mol_name_to_mol_dict = {}
        all_mol_name_set = set()
        for i in self.steps:
            for s in i.reactant_list:
                mol_name_to_mol_dict[str(s)] = s
                all_mol_name_set.add(str(s))
            for s in i.product_list:
                mol_name_to_mol_dict[str(s)] = s
                all_mol_name_set.add(str(s))
        all_mol = [mol_name_to_mol_dict[i] for i in all_mol_name_set]
        return all_mol, all_mol_name_set, mol_name_to_mol_dict

    def check_product_in_dict(self,check=False):
        if check:
            product_not_involed=[]
            all_mol, all_mol_name,name_to_mol_dict=self.get_all_involved_mols_and_mol_name_to_mol_dict()
            for i in self.product_name_showed:
                if i in all_mol_name:
                    pass
                else:
                    product_not_involed.append(i)
            if len(product_not_involed)==0:
                print('All products in final results are involed')
            else:
                print('%s is not involed in the network ,check if the condition is true' % product_not_involed)
        return

    def species_do_one_reactant_one_product_reaction_until_no_new_product(self, reactions):
        for s in self.species_pool:
            s.sign_one_reactant_one_product_reactions(reactions)#set the flag to false, so we will consider all the reactions
        while 1: #这里不断的迭代，去让中间体不断的往下走，然后不断扩充species_pool，直到产生不出new product
            new_product = []
            for s in self.species_pool:
                result, reaction_names = s.react_one_reactant_one_product(
                    new_oto_reactions_to_sign_for_product=reactions)
                if result != ReactFlags.NoProduct:
                    for j in range(len(result)):#
                        r = result[j]
                        # if showed, only add reaction, but not add to species pool!
                        if str(r) in self.species_name_showed:
                            name = reaction_names[j]
                            self.steps.append(ReactionStep([copy.deepcopy(s)], [copy.deepcopy(r)], name))
                            continue#as the mol have already been generated, we don't add it to pool, but its reaction is new, so we add it
                        if self.check_if_mol_valid(r) == False:
                            continue
                        self.species_pool.append(r)
                        self.species_name_showed.append(str(r))
                        name = reaction_names[j]
                        self.steps.append(ReactionStep([copy.deepcopy(s)], [copy.deepcopy(r)], name))
                    new_product.extend(result)

            if len(new_product) == 0:
                break

    def __len__(self):
        return len(self.steps)

    def remove_repeated_reaction_regardless_direction(self, debug=False):
        if len(self.steps) == 0:
            print("No reaction!")
            return
        print("Reaction num before remove repeated reaction %s" % len(self.steps))

        unique_steps = {}
        for i in self.steps:
            if debug:
                if i.hash_regardless_direction in unique_steps.keys():
                    print(
                        "Get repeated reaction %s, will keep %s" % (unique_steps[i.hash_regardless_direction], str(i)))
            unique_steps[i.hash_regardless_direction] = i
        self.steps = list(unique_steps.values())
        print("Reaction num after remove repeated reaction %s" % len(self.steps))

    def remove_repeated_reaction_regardless_direction_change(self):
        if len(self.steps) == 0:
            print("No reaction!")
            return
        print("Reaction num before remove repeated reaction %s" % len(self.steps))

        unique_steps = {}
        for i in self.steps:
            if str(i) not in [
                # '[H]C([H])[H] > Add C on first carbon > [H]C([H])([H])C',
                # '[H]C([H])[H] > Add C on last carbon > [H]C([H])([H])C',
                # '[H]C([H])[H] > Add CO on first carbon > [H]C([H])([H])CO',
                # '[H]C([H])[H] > Add CO on last carbon > [H]C([H])([H])CO',
                # '[H]C([H])[H] > Add [H]C on first carbon > [H]CC([H])([H])[H]',
                # '[H]C([H])[H] > Add [H]C on last carbon > [H]CC([H])([H])[H]',
                # '[H]C([H])[H] > Add [H]C[H] on first carbon > [H]C([H])C([H])([H])[H]',
                # '[H]C([H])[H] > Add [H]C[H] on last carbon > [H]C([H])C([H])([H])[H]',
                # '[H]C[H] > Add C on first carbon > [H]C([H])C',
                # '[H]C[H] > Add C on last carbon > [H]C([H])C',
                # '[H]C[H] > Add CO on first carbon > [H]C([H])CO',
                # '[H]C[H] > Add CO on last carbon > [H]C([H])CO',
                # '[H]C[H] > Add [H]C on first carbon > [H]CC([H])[H]',
                # '[H]C[H] > Add [H]C on last carbon > [H]CC([H])[H]',
                # '[H]C > Add C on first carbon > [H]CC',
                # '[H]C > Add C on last carbon > [H]CC',
                # '[H]C > Add CO on first carbon > [H]CCO',
                # '[H]C > Add CO on last carbon > [H]CCO',
                # 'CO > Add C on first carbon > CCO',
                # 'CO > Add C on last carbon > CCO',
            ] :
                if i.hash_regardless_direction in unique_steps.keys():
                    print(
                        "Get repeated reaction %s, will keep %s" % (unique_steps[i.hash_regardless_direction], str(i)))
                unique_steps[i.hash_regardless_direction] = i
        self.steps = list(unique_steps.values())
        print("Reaction num after remove repeated reaction %s" % len(self.steps))




    def get_total_readable_reaction_string(self,color = 1):

        string = ""
        for s in self.steps:
            s.color_string=color
            string += str(s) + "\n"
        return string


    def group_reaction_by_reaction_name(self):
        grouped_steps = []
        rx_name_to_rx_dict = {}
        for s in self.steps:
            try:
                rx_name_to_rx_dict[s.reaction_name].append(s)
            except KeyError:
                rx_name_to_rx_dict[s.reaction_name] = [s]
        for i in rx_name_to_rx_dict:
            grouped_steps.extend(rx_name_to_rx_dict[i])
        self.steps = grouped_steps

    def sort_reaction_by_reaction_str(self):
        self.steps.sort(key=lambda i: str(i))

    def add_custom_step(self, reactant, product, reaction_name):
        assert isinstance(reactant, SpeciesNode)
        assert isinstance(product, SpeciesNode)
        assert isinstance(reaction_name, str)
        self.steps.append(ReactionStep([reactant], [product], reaction_name))

    def check_if_mol_valid(self, species_node):
        if self.valid_mol_check_func_list is None:
            return True
        for func in self.valid_mol_check_func_list:
            if func(species_node.mol) == False:
                return False
        return True


    def reaction_to_catmap_usable(self):
        raise NotImplementedError()

    #从这里开始放新的
    def reaction_to_allmol_in_reactionlist(self):
        if self.add_custom_step() != None:
            pass
            # return .reactionlist
    # 只要指定两种产物即可 和反应物的道理类似
    # mdzz



    def plot_network_to_pdf_file(self, filepath="network.pdf", ):
        def make_gviz_valid_name(specie_name):
            return specie_name.replace("[", "a").replace("]", "b").replace("(", "c").replace(")", "d").replace("=", "e")

        def make_specie_str(specie):
            output = make_gviz_valid_name(specie) + "[label=<" + specie
            output += "> fontsize=80 penwidth=1 arrowsize=10 , fontname=\"arial\" ];\n"
            # the width of edge below
            return output

        print("Plotting Network with Smiles coding...")
        start_str = "digraph html {"
        define_str = ""
        node_str = ""

        species = []
        new_el_rxns = []

        for step in self.steps:
            species.extend([str(i) for i in step.reactant_list])
            species.extend([str(i) for i in step.product_list])

            reaction = [[str(i) for i in step.reactant_list],
                        [str(i) for i in step.product_list]]

            new_el_rxns.append(reaction)

        species = set(species)
        for specie in species:
            define_str += make_specie_str(specie)

        for reaction in new_el_rxns:
            for n in range(len(reaction) - 1):
                for path in product(reaction[n], reaction[n + 1]):
                    if set(path) < set(species):
                        if random.random() > 0.5:
                            path1 = path[0]
                            path2 = path[1]
                        else:
                            path1 = path[1]
                            path2 = path[0]

                        node_str += make_gviz_valid_name(path1) + "->" + make_gviz_valid_name(path2)
                        node_str += "[penwidth=5 arrowType=open arrowsize=2 arrowhead=\"none\"]" + "\n"

        total_str = start_str + define_str + node_str + "}"
        textfilepath = filepath.split(".")[0] + ".txt"

        with open(textfilepath, "w") as f:
            f.write(total_str)

        os.system("dot -Tpdf %s -o %s" % (textfilepath, filepath))

    def plot_network_to_pdf_file_labeltransed(self, filepath="network.pdf", ):
            def make_gviz_valid_name(specie_name):
                return specie_name.replace("[", "a").replace("]", "b").replace("(", "c").replace(")", "d").replace("=", "e")

            def make_specie_str(specie):
                output = make_gviz_valid_name(specie) + "[label=<" + expression_transform(specie)
                output += "> fontsize=80 penwidth=1 arrowsize=10 , fontname=\"arial\" ];\n"
                # the width of edge below
                return output

            print("Plotting Network with transformed expression...")
            start_str = "digraph html {"
            define_str = ""
            node_str = ""

            species = []
            new_el_rxns = []

            for step in self.steps:
                species.extend([str(i) for i in step.reactant_list])
                species.extend([str(i) for i in step.product_list])

                reaction = [[str(i) for i in step.reactant_list],
                            [str(i) for i in step.product_list]]

                new_el_rxns.append(reaction)

            species = set(species)
            for specie in species:
                define_str += make_specie_str(specie)

            for reaction in new_el_rxns:
                for n in range(len(reaction) - 1):
                    for path in product(reaction[n], reaction[n + 1]):
                        if set(path) < set(species):
                            if random.random() > 0.5:
                                path1 = path[0]
                                path2 = path[1]
                            else:
                                path1 = path[1]
                                path2 = path[0]

                            node_str += make_gviz_valid_name(path1) + "->" + make_gviz_valid_name(path2)
                            node_str += "[penwidth=5 arrowType=open arrowsize=2 arrowhead=\"none\"]" + "\n"

            total_str = start_str + define_str + node_str + "}"
            textfilepath = filepath.split(".")[0] + ".txt"

            with open(textfilepath, "w") as f:
                f.write(total_str)

            os.system("dot -Tpdf %s -o %s" % (textfilepath, filepath))