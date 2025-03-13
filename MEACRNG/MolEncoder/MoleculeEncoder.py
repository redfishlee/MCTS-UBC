from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from MEACRNG.Tools.MolViz import show_mol
import string


class OrganicSurfaceSpeciesEncoder:
    def __init__(self,exclude_COC=True):
        '''
        :param mol_string:
        :param exclude_COC: # 为True时COC当成OCC，否则当做醚处理
        '''
        '''
        Encoder暂时只能以C作为骨架原子，若遇到能作骨架的其他原子比如N, S
        则需要注意，一般会认为是支链，但有时候会出现bug
        '''
        self.exclude_COC = exclude_COC

    @staticmethod
    def add_atom_from_left_to_right(mol_string):

        Head_Carbon = mol_string.find('C')
        l = list(mol_string)
        # 建立结构
        structure = Chem.MolFromSmarts("C")
        structure = Chem.RWMol(structure)

        carbon_position = 0
        position = 0
        # 检查左边是否有原子
        if Head_Carbon != 0:
            nl = l[:Head_Carbon]
            nl.reverse()
            l = l[Head_Carbon:]

            # 对左边原子进行处理
            while len(nl) != 0:  # 当首字母不是C，那么左边有原子，执行左边加原子。若首字母已经C，不执行此循环。过程中头碳位置始终为0.
                if nl[0] in string.digits:
                    atomstring = nl.pop(1)
                    atomstring = '[' + atomstring + ']'  # smarts不支持单个H，需要[H]
                    atom = Chem.MolFromSmarts(atomstring)
                    atom = atom.GetAtomWithIdx(0)
                    for count in range(int(nl[0])):
                        structure.AddAtom(atom)
                        structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                    position = structure.GetNumAtoms() - 1
                    nl.pop(0)
                elif nl[0] in string.ascii_uppercase:
                    atomstring = nl.pop(0)
                    atomstring = '[' + atomstring + ']'
                    atom = Chem.MolFromSmarts(atomstring)
                    atom = atom.GetAtomWithIdx(0)
                    structure.AddAtom(atom)
                    structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                    if atomstring == '[H]':
                        position = carbon_position
                    else:
                        position = structure.GetNumAtoms() - 1
                    if nl != []:
                        if nl[0] == 'O':
                            position = carbon_position
                elif nl[0] in string.punctuation:
                    if nl[0] == '-':
                        nl.pop(0)
                        atomstring = nl.pop(0)
                        atomstring = '[' + atomstring + ']'
                        atom = Chem.MolFromSmarts(atomstring)
                        atom = atom.GetAtomWithIdx(0)
                        structure.AddAtom(atom)
                        structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.AROMATIC)
                        position = structure.GetNumAtoms() - 1
            position = 0
            carbon_position = 0
        # 开始加右边
        # if mol.find('COC') != -1:
        #     inter_mol = ''.join(l)
        #     inter_mol = inter_mol.replace('COC','OCC')
        #     l = list(inter_mol)
        # #为了符合文献中COC==OCC.不过又得右边加了，这样不行
        l.pop(0)  # 创建初始结构后，敲掉首字母位置的碳，该位置为0
        while len(l) != 0:  # 该循环最大注意点是position变化
            if l[0] in string.digits:  # 如果是数字，就代表将前一个原子add到前一个原子位置n次，因此输入position无需变化
                numberstring = l.pop(0)  # 获得该数字具体数值
                for count in range(int(numberstring) - 1):  # 加n次
                    structure.AddAtom(atom)
                    structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
            elif l[0] in string.punctuation:  # 如果是括号，过渡态或者支链
                origin_position = carbon_position  # 由于加完括号内东西后，position应该回到原先的C上，因此记录下先前点
                if mol_string.find('(') != -1 and mol_string.find(')') != -1:
                    while l[0] != ')':  # 该环节运行到)为止
                        if l[0] == '(':  # (直接敲掉
                            l.pop(0)
                        elif l[0] == '-':  # 因为-后面是要加上去的原子，敲掉-后直接开始加原子，键使用
                            l.pop(0)  # 敲掉-
                            atomstring = l.pop(0)  # 得到原子，并且敲掉
                            atomstring = '[' + atomstring + ']'  # smarts不支持单个H，需要[H]
                            atom = Chem.MolFromSmarts(atomstring)
                            atom = atom.GetAtomWithIdx(0)
                            structure.AddAtom(atom)
                            structure.AddBond(position, structure.GetNumAtoms() - 1,
                                              Chem.BondType.AROMATIC)  # 此处position应是之前的C
                            previous_atom = atom
                            position = structure.GetNumAtoms() - 1  # 由于加完这个原子后，如果后面还有原子，应加在这个原子后面。
                            if atomstring == '[C]':
                                position = structure.GetNumAtoms() - 1
                                carbon_position = structure.GetNumAtoms() - 1
                            elif atomstring == '[H]':  # 加完H后回到C
                                position = carbon_position
                            else:
                                if l != []:
                                    if l[0] in string.ascii_uppercase:
                                        if l[0] == 'O':
                                            position = structure.GetNumAtoms() - 2  # 专门为羧基写的，存在问题，是否还有其他基团需要特别照顾
                                        else:
                                            position = structure.GetNumAtoms() - 1
                        else:
                            if l[0] in string.ascii_uppercase:  # 加完过渡态原子后，接着加之后的原子，如果没有的话，会摸到)，直接结束循环
                                atomstring = l.pop(0)
                                atomstring = '[' + atomstring + ']'  # smarts不支持单个H，需要[H]
                                atom = Chem.MolFromSmarts(atomstring)
                                atom = atom.GetAtomWithIdx(0)
                                structure.AddAtom(atom)
                                structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                                previous_atom = atom
                                if atomstring == '[C]':  # 如果下一个原子是C，则需要调整carbon_position
                                    position = structure.GetNumAtoms() - 1
                                    carbon_position = structure.GetNumAtoms() - 1  # 将下一个原子的连接位置调整到现在这个C上
                                elif atomstring == '[H]':  # H之后不可能有东西，所以加完H后回到C。#####################可能问题：如果有其他原子具有类似C的功能，比如Si，那么将无法编译
                                    position = carbon_position
                                else:
                                    if l != []:
                                        if l[0] in string.ascii_uppercase:
                                            if l[0] == 'O':
                                                position = structure.GetNumAtoms() - 2  # 专门为羧基写的。########################存在问题，是否还有其他基团需要特别照顾
                                            else:
                                                position = structure.GetNumAtoms() - 1
                            elif l[0] in string.digits:
                                numberstring = l.pop(0)
                                for count in range(int(numberstring) - 1):
                                    structure.AddAtom(atom)
                                    structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                    l.pop(0)
                    position = origin_position
                else:
                    if l[0] == '-':
                        l.pop(0)
                        atomstring = l.pop(0)
                        atomstring = '[' + atomstring + ']'
                        atom = Chem.MolFromSmarts(atomstring)
                        atom = atom.GetAtomWithIdx(0)
                        structure.AddAtom(atom)
                        structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.AROMATIC)
                        previous_atom = atom
                        position = structure.GetNumAtoms() - 1
                        if atomstring == '[C]':
                            position = structure.GetNumAtoms() - 1
                            carbon_position = structure.GetNumAtoms() - 1
                        elif atomstring == '[H]':  # 加完H后回到C
                            position = carbon_position
                        else:
                            if l != []:
                                if l[0] in string.ascii_uppercase:
                                    if l[0] == 'O':
                                        if l[1] == 'O':
                                            position = carbon_position  # 专门为羧基写的，存在问题，是否还有其他基团需要特别照顾
                                        else:
                                            position = structure.GetNumAtoms() - 1
                                    else:
                                        position = structure.GetNumAtoms() - 1

            elif l[0] in string.ascii_uppercase:  # 正常加原子
                atomstring = l.pop(0)
                atomstring = '[' + atomstring + ']'  # smarts不支持单个H，需要[H]
                atom = Chem.MolFromSmarts(atomstring)
                atom = atom.GetAtomWithIdx(0)
                structure.AddAtom(atom)
                structure.AddBond(position, structure.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                previous_atom = atom
                if atomstring == '[C]':  # 如果下一个原子是C，则需要调整carbon_position
                    position = structure.GetNumAtoms() - 1
                    carbon_position = structure.GetNumAtoms() - 1  # 将下一个原子的连接位置调整到现在这个C上
                elif atomstring == '[H]':  # H之后不可能有东西，所以加完H后回到C。#####################可能问题：如果有其他原子具有类似C的功能，比如Si，那么将无法编译
                    position = carbon_position
                else:
                    if l != []:
                        if l[0] in string.ascii_uppercase:
                            if l[0] == 'O':
                                position = carbon_position  # 专门为羧基写的。########################存在问题，是否还有其他基团需要特别照顾
                            else:
                                position = structure.GetNumAtoms() - 1
                        elif l[0] == '-':
                            if atomstring == '[C]':
                                position = carbon_position
                            elif atomstring == '[H]':
                                position = carbon_position
                            elif atomstring == '[O]' and l[1] == 'H':
                                position = structure.GetNumAtoms() - 1
                            elif atomstring == '[O]' and l[1] == 'O':
                                position = carbon_position

        return structure

    def encode(self,mol_string):

        # 排除非有机物
        if mol_string.find('C') == -1:
            error = "Got error mol string %s, only support organics with carbon structure!" % mol_string
            raise ValueError(error)

        if self.exclude_COC == True:
            if mol_string.find('COC') != -1:
                raise ValueError(
                    "Got mol string %s may misunderstanding, ether is not allowed, COC will be recognized as ether"
                    "\nif you want to encode ether, please set 'exclude_COC' to False." % mol_string)
        mol_string = mol_string.replace('*', '')

        return self.add_atom_from_left_to_right(mol_string)


    def show_graph(self,mol_string):
        structure = self.add_atom_from_left_to_right(mol_string)
        # 纠正显示错误
        structure = Chem.MolFromSmiles(Chem.MolToSmiles(structure))
        structure = Draw.MolToImage(structure)
        structure.show()

if __name__ == '__main__':

    encoder = OrganicSurfaceSpeciesEncoder(True)
    mol = encoder.encode("CH2(CHCH3)OH")
    # mol = encoder.encode("CH2(CHCH3)NH2OH")这种就无法正确显示
    show_mol(mol)
