from MEACRNG.Tools.redfish_functionalities import get_the_lost_reactant
from MEACRNG.Tools.Smiles2formula import expression_transform# 编码转换的工具

import json
import ast

with open("network_with_energy.json", "r") as file:
	network_with_energy = json.load(file)
#print(network_with_energy)
mkm_file = "mkm_input.m"

def generating_mkm_input(network_with_energy, mkm_file):
	step_label = 0
	with open(mkm_file, 'w') as f:
		f.write("%(x): for note    Energy input:[ Ea | H0 ]\n")
		for reaction, energy in network_with_energy.items():
			step_label += 1
			#constructe the step description
			species_list = ast.literal_eval(reaction)
			specie1 = species_list[0]
			specie2 = species_list[1]
			lost_reactant = get_the_lost_reactant(specie1, specie2)

			specie1_length = len(specie1.replace('[','').replace(']','').replace('(','').replace(')',''))
			specie2_length = len(specie2.replace('[','').replace(']','').replace('(','').replace(')',''))
			formula1 = expression_transform(specie1)
			formula2 = expression_transform(specie2)

			if specie1_length > specie2_length:
				longer_one = formula1 
				shorter_one = formula2
			if specie1_length < specie2_length:
				longer_one = formula2
				shorter_one = formula1
			side1 = longer_one
			side2 = [shorter_one, lost_reactant]
			print(side2, side1)

			#get the step description of each step
			step_description = get_step_description(side2, side1)
			whole_description = f"({step_label})" + ":" + step_description + f"[ {energy[0]}  {energy[1]} ]" + '\n'
			f.write(whole_description)

	i = step_label + 1
	gas_reaction_description = []
	gas_reaction_description.append(f"({i}): CO+# <-> CO#    [ 0.000 -0.922]\n")
	i += 1
	gas_reaction_description.append(f"({i}): H2+#+# <-> H#+H#    [0.000 -0.694]\n")
	i += 1
	gas_reaction_description.append(f"({i}): CH3CH2CH3# <-> CH3CH2CH3+#    [0.321 -0.159]\n")
	i += 1
	gas_reaction_description.append(f"({i}): H2O+O#+# <-> HO#+HO#    [0.221 -0.143]\n")
	i += 1
	gas_reaction_description.append(f"({i}): H2O+#+# <-> HO#+H#    [0.000 -0.542]\n")
	with open(mkm_file, 'a') as f:
		for i in range(len(gas_reaction_description)):
			f.write(gas_reaction_description[i])
	other_setting = []
	other_setting.append('CalcDRC = 1;')
	other_setting.append('PlotType = 0')
	other_setting.append('PlotMode = 0')
	other_setting.append('T = 500;')
	other_setting.append('Q_v_INIT = 1;')
	other_setting.append('npar = 8;')
	other_setting.append('P_H2_FROZ = 6.67;')
	other_setting.append('P_CO_FROZ = 3.33;')
	other_setting.append('P_CH3CH2CH3_FROZ = 0.457;')
	other_setting.append('P_H2O_FROZ = 0.01;')
	with open(mkm_file, 'a') as f:
		for i in range(len(other_setting)):
			f.write(f"{other_setting[i]}\n")

#get the step description of each step
def get_step_description(specie1, specie2):
	step_description = specie1[0] + '#+' + specie1[1] + '#' + '<->' + specie2 + '#+#' + '  '
	return step_description

generating_mkm_input(network_with_energy, mkm_file)