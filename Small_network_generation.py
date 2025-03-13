import json
import random

def reactions_in_pathway(pathway):
	reactions = []
	for index_of_pathways in range(3):
		for index_of_specie in range(len(pathway[index_of_pathways])-1):
			current_reaction = (pathway[index_of_pathways][index_of_specie], pathway[index_of_pathways][index_of_specie+1])
			if current_reaction in reactions:
				continue
			reactions.append(current_reaction)
	return reactions

with open("All_pathways.json", "r") as file:
	all_pathways = json.load(file)
assert isinstance(all_pathways, dict), "all_pathway file doesn't contain a dictionary"
assert len(all_pathways) == 100, "all_pathway file has pathways less that 100"
length_of_pahtway = {}
for index_of_whole_pathway, key_of_whole_pathway in enumerate(all_pathways):
	length_of_pahtway[index_of_whole_pathway] = len(reactions_in_pathway(all_pathways[key_of_whole_pathway]))

sorted_keys = sorted(length_of_pahtway, key=length_of_pahtway.get)[:10]
unique_reactions = set()
for index in sorted_keys:
	print(index, length_of_pahtway[index])
	print(reactions_in_pathway(all_pathways[f"{index+1}"]))
	unique_reactions.update(reactions_in_pathway(all_pathways[f"{index+1}"]))
print(unique_reactions)
print(len(unique_reactions))

network_with_energy = {}
barrier_energy = 0
for reaction in unique_reactions:
	reaction_energy = round(random.uniform(0, 2), 3)
	network_with_energy[f"{reaction}"] = [reaction_energy, barrier_energy]  

print(network_with_energy)
print(len(network_with_energy))
with open("network_with_energy.json", "w") as file:
	json.dump(network_with_energy, file, indent=4)



























