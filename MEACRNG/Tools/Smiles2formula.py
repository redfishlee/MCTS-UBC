import string

def txt_wrap_by(start_str, end, target):
    start = target.find(start_str)
    if start >= 0:
        start += len(start_str)
        end = target.find(end, start)
        if end >= 0:
            return target[start:end].strip()

def expression_transform(immaformula):
    '''
    TODO
    This function is still in test, and it supports formula which only contain C H O at present.
    '''

    available_atom = ['C','H','O']# skeletal atom list. There is an order in this atom list, Carbon usually enjoys the first order.
    # check if there are atoms in immaformula beyond available atom.
    for element in immaformula:
        if element in string.punctuation:
            continue
        if element in string.digits:
            continue
        if element in string.ascii_letters:
            if element not in available_atom:
                return immaformula
    '''
    this function only support Carbon contained linear molecule temporary. 
    '''
    if 'C' not in immaformula:
        return immaformula
    All_Ts_dic = {'C':{'HCH3':'CH4', 'HCH2':'CH3', 'HCH':'CH2', 'HC':'CH'}, 'O':{'HO':'OH', 'HOH':'H2O'}, }
    Ts_dic = All_Ts_dic['C']

    i = immaformula
    i = i.replace('([H])','H').replace('[H]','H')
    i = i.replace('HHH','H3').replace('HH','H2')
    # print(i)
    if i.count('C') == 1:# C1单独处理
        for sp in Ts_dic.keys():
            i = i.replace(sp,Ts_dic[sp])
        return i
    else:# Cn情况
        result = ''
        while i.find('C') != -1:
            current_first_carbon = i.find('C') 
            content = txt_wrap_by('C','C',i)
            current_second_carbon = i.find('C',i.find('C')+1)

            if current_second_carbon == -1:
                result += i
                i = ''
            else:
                mid = i[:current_first_carbon+1] + content
                if string.punctuation not in mid == True:# 无支链情况
                    if mid in Ts_dic.keys():
                        mid = Ts_dic[mid]
                    result += mid
                    i = i[current_second_carbon:]
                else:
                    left_br_num = mid.count('(')
                    right_br_num = mid.count(')')
                    if left_br_num == right_br_num:# C(O)(O)(O)C情况，将()内容取出后，进行变换，然后再加回去
                        backup = ''
                        for iter_cyc in range(left_br_num):
                            mid_mid = '(' + txt_wrap_by('(',')',mid) + ')'
                            mid = mid.replace(mid_mid,'',1)
                            backup += mid_mid
                        if mid in Ts_dic.keys():
                            mid = Ts_dic[mid]
                        result += mid + backup
                        i = i[current_second_carbon:]
                    if left_br_num - 1 == right_br_num:# C(O)(O)(O)(C)CCH3情况，也就是content会变成(O)(O)(O)(，去掉(后同C(O)(O)(O)C情况
                        mid = mid.rstrip('(')
                        backup = ''
                        for iter_cyc in range(right_br_num):
                            mid_mid = '(' + txt_wrap_by('(',')',mid) + ')'
                            mid = mid.replace(mid_mid,'',1)
                            backup += mid_mid
                        if mid in Ts_dic.keys():
                            mid = Ts_dic[mid]
                        result += mid + backup + '('
                        i = i[current_second_carbon:]
                    if left_br_num + 1 == right_br_num:# C)CH3情况，也就是C(O)(O)(O)(C)CCH3第一次循环处理完后的结果，直接正常运行即可
                        if mid in Ts_dic.keys():
                            mid = Ts_dic[mid]
                        result += mid
                        i = i[current_second_carbon:]
        # print(result)
        return result
 


    

if __name__ == '__main__':


    with open(r'C:\Users\OldKuroCat\Desktop\private\me\T2\C4mollist.txt','r') as f:
        test = f.readlines()

    atom_key_list = ['C']

    # test = ['[H]C([H])(O)(O)(C)C([H])([H])[H]',]
    test = ['[F]']

    species = []
    for i in test:
        species.append(i.strip('\n'))


    final_state = []
    for i in species:
        # print(i)
        for atom in atom_key_list:
            final_state.append(expression_transform(i))
        print(expression_transform(i))
        # print('-----------------------------------------------')

    # for i in final_state:
    #     print(i)