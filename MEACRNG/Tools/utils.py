'''
reference: http://www.rdkit.org/docs/GettingStartedInPython.html
when coding using RDKit, search for document is very frequent!
To find the class of the instance you create, and find its method,
and to fine the class of the instance what its method returns.
'''
import matplotlib.pyplot as plt
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import math
import matplotlib.pyplot as plt




def get_index_and_total_bond_num_of_the_nth_atom_of_element_with_total_m_bonds(# TOUNDER
        mol,
        n=1,
        element="C",
        m=-1):
    '''
    
    获得某个原子的index以及总共有多少个键，

    注意输入的mol并不是字符串，而是编码过的分子
    这个原子是mol中的第n个： 元素为element，并且有m个bond  的原子
    e.g.
    n=1,element="C",m=4: get first C with 4 bonds
    n=3,element="O",m=2: get third O with 2 bonds
    n=3,element="O",m=-1: get third O with any bonds
    return the index of that atom
    '''

    if n <= 0:
        raise ValueError("n >= 1")# 注意valueerror和typeerror的区别，前者是指接受了可行的参数类型，但参数不符合要求
    atoms = mol.GetAtoms()# 从rdkit编码的分子中使用getatoms得到其中原子，原子仍然使用的是rdkit的编码
    nth_atom = 0# 初始化一个标记
    for i in atoms:
        if i.GetSymbol() == element:# 想要atom能够和字符串原字符号对应，必须对rdkit编码的原子再使用getsymbol
            total_bond_num = get_total_bond_num_regarding_triple_double_single_as_3_2_1(i)# 获取该原子上总键数
            if m == -1:# m=-1代表任何键，前面的注释中有解释
                nth_atom += 1
                if nth_atom >= n:# 如果找不到符合条件的情况，比如在CH4中想找到第2个连接任意键数的C，最后就会return none
                    return i.GetIdx(), total_bond_num# 获取原子在原分子的index位置
            else:
                if len(i.GetBonds()) == m:# getbonds返回的是该原子上每一根键，不管所谓的键型是否相同
                    nth_atom += 1
                    if nth_atom >= n:
                        return i.GetIdx(), total_bond_num
    return None, None# 如果找不到符合条件的情况（可能由于很多问题导致），就返回index none，总键数 none


class ReactionLibUtils:# 合并可能发生步骤的class
    @staticmethod
    def merge_reactions(list_of_reaction_list):
        all_reactions = []
        for i in list_of_reaction_list:
            all_reactions.extend(i)
        return all_reactions


class ReactFlags:# TOUNDER
    NoProduct = "NoProduct"


def get_total_bond_num_regarding_triple_double_single_as_3_2_1(atom):# 统计rdkit编码的某个原子身上的键数，三键认为是3根键
    total_bond_num = 0# 初始化
    for b in atom.GetBonds():# 对rdkit编码的原子进行getbonds，一个原子上不止一根，甚至不止一种键
        if b.GetBondType() == Chem.BondType.SINGLE:# 检查键类型，具体有哪些类型可以到rdkit官方文档上查找
            total_bond_num += 1
        elif b.GetBondType() == Chem.BondType.DOUBLE:
            total_bond_num += 2
        elif b.GetBondType() == Chem.BondType.TRIPLE:
            total_bond_num += 3
        else:
            raise ValueError("Do not support such bond type!")# 暂时只需要统计单双三键
    return total_bond_num


def get_num_of_single_double_triple_bonds_and_number_of_atom_B_with_centerd_atom_A(B_symbol, A_atom):# 统计某个选定的中心原子上指定的某种原子个数，之间化学键的数量和类型
    '''
    注意输入的B是字符串symbol，A则是rdkit编码的原子
    e.g.: B="O", A=Chem.Mol("C=O(-O)"),
    will return:
    2(number of O),
    1(number of single bond),
    1(number of double bond),
    0(number of triple bond)
    '''
    num_B = 0
    num_single = 0
    num_double = 0
    num_triple = 0

    for b in A_atom.GetBonds():
        if b.GetBeginAtom().GetSymbol() == B_symbol or b.GetEndAtom().GetSymbol() == B_symbol:# 所谓的键头键尾原子与encode时的方式有关(所以要用or)
            num_B += 1
            if b.GetBondType() == Chem.BondType.SINGLE:
                num_single += 1
            elif b.GetBondType() == Chem.BondType.DOUBLE:
                num_double += 1
            elif b.GetBondType() == Chem.BondType.TRIPLE:
                num_triple += 1
            else:
                raise ValueError("Can not cover other bonds except single double triple!")
    return num_B, num_single, num_double, num_triple







def color_string(strPrint, front_color="", back_color=""):# 使得输出字符串带有一定颜色以及格式。注意如果使用的控制台无法识别终端转义码，则会直接当做一般的字符串print，比如sublime有这个问题，pycharm则能正常运作
    '''
    linux中
    echo -e "\033[字背景颜色;字体颜色m字符串\033[0m"  

    '''

    fd = {
        "BLACK": 30,
        'RED': 31,
        'GREEN': 32,
        'BROWN': 33,
        'BLUE': 34,
        'PURPLE': 35,
        'CYAN': 36,
        'WHITE': 37,
        'UNDERLINE_ON': 38,
        'UNDERLINE_OFF': 39,
        "": ""}

    bd = {
        "BLACK": 40,
        'RED': 41,
        'GREEN': 42,
        'BROWN': 43,
        'BLUE': 44,
        'PURPLE': 45,
        'CYAN': 46,
        'WHITE': 47,
        "": ""

    }

    strPrefix = "\033[0;%s;%sm" % (fd[str.upper(front_color)], bd[str.upper(back_color)])

    strSuffix = "\033[0m"
    strMsg = strPrefix + strPrint + strSuffix

    return strMsg


def get_all_bond_as_edge_in_mol(mol):# 得出一个分子的原子位置以及之间的连键信息
    edges = set()
    for a in mol.GetAtoms():
        for b in a.GetBonds():
            edges.add((b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx(), b.GetBondType()))
    return edges

