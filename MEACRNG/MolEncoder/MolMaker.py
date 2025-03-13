from rdkit.Chem import AllChem as Chem
from MEACRNG.Tools import utils as u


class MMConnectConfig:# 建立键类型信息，another_atom_atomic_number代表连键原子的原子序数，Hs_to_add_to_another_atom代表有多少个H加在连键原子上，比如OH就是1
    def __init__(self, bond_type=Chem.BondType.SINGLE, another_atom_atomic_number=1, Hs_to_add_to_another_atom=0):
        self.bond_type = bond_type
        self.another_atom = another_atom_atomic_number
        self.Hs_to_add_to_another_atom = Hs_to_add_to_another_atom



class MolMaker:

    @staticmethod
    def connect_bonds_from_config_to_first_atom_A_in_initial_mol(
            initial_mol,
            mmcc_list,
            A_symbol="C"):
        '''
        给一个initial mol分子，在这个分子的第一个symbol为A的原子上按照mmcc的设置添加原子
        '''
        assert isinstance(mmcc_list, list)
        for a in initial_mol.GetAtoms():
            if a.GetSymbol() == A_symbol:
                for c in mmcc_list:
                    assert isinstance(c, MMConnectConfig)
                    initial_mol.AddAtom(Chem.Atom(c.another_atom))#typical workflow of adding atom, add atom, find index, add bond
                    index_of_add_atom = initial_mol.GetNumAtoms() - 1
                    initial_mol.AddBond(a.GetIdx(), index_of_add_atom, c.bond_type)
                    for _ in range(c.Hs_to_add_to_another_atom):#add H to atoms, doesn't matter that much
                        initial_mol.AddAtom(Chem.Atom(1))
                        initial_mol.AddBond(index_of_add_atom, initial_mol.GetNumAtoms() - 1, Chem.BondType.SINGLE)
                return initial_mol

        return False


    @staticmethod
    def add_mol_B_to_mol_A_and_make_bonds_in_ith_atom_in_A_with_jth_atom_in_B(
            atomA,
            atomB,
            idx_of_atom_in_A_to_make_bond,
            idx_of_atom_in_B_to_make_bond,
            bond_type=Chem.BondType.SINGLE
    ):
        '''
        将AB两个分子加在一起，同时A中的第i个原子和B中的第j个原子建立bond

        算法：
        先获得B中所有的bond的信息，得到边，然后add atom，把所有index都加上A原子的个数（因为加上过后按照顺序排是+1）
        最后把bond加上
        '''
        atomA = Chem.RWMol(atomA)
        atomB = Chem.RWMol(atomB)
        all_bonds_in_B = u.get_all_bond_as_edge_in_mol(atomB)
        num_offset_when_added_B = atomA.GetNumAtoms()
        # 添加B的atom
        for a in atomB.GetAtoms():
            atomA.AddAtom(a)

        # 添加B的键
        for i in all_bonds_in_B:
            atomA.AddBond(i[0]+num_offset_when_added_B,i[1]+num_offset_when_added_B,i[2])
        # 添加AB atom的键
        atomA.AddBond(idx_of_atom_in_A_to_make_bond, idx_of_atom_in_B_to_make_bond+num_offset_when_added_B,bond_type)
        return atomA





def test_add_two_mol():
    import MEACRNG.MolEncoder.C1MolLib as c

    mol1 = c.C1MolLib.COOH_single_CO_bond
    # u.show_mol(mol1)

    mol2 = c.C1MolLib.HCOO_double_CO_bond

    C_index_in_CHO, _ = u.get_index_and_total_bond_num_of_the_nth_atom_of_element_with_total_m_bonds(mol1, 1, "C")
    C_index_in_CHO_2, _ = u.get_index_and_total_bond_num_of_the_nth_atom_of_element_with_total_m_bonds(mol2, 1, "C")

    new_mol = MolMaker.add_mol_B_to_mol_A_and_make_bonds_in_ith_atom_in_A_with_jth_atom_in_B(
        mol1,
        mol2,
        C_index_in_CHO,
        C_index_in_CHO_2)
    u.show_mol(new_mol)









if __name__ == '__main__':
    test_add_two_mol()


