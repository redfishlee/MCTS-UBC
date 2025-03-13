import json
import ast

with open("network_with_energy.json", "r") as file:
	network_with_energy = json.load(file)
print(network_with_energy)
for key, value in network_with_energy.items():
	print(key)
	key_list = ast.literal_eval(key)
	print(type(key_list))
	print(key_list[0])
	print(key_list[1])