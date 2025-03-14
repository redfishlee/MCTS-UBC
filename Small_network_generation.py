import json
import random
from MEACRNG.Tools.utils import ReactionLibUtils



reactions = ReactionLibUtils.merge_reactions([
    C3OneReactionOneProductReactionLib.reaction_add_carbon_species_on_jth_carbon,
    C3OneReactionOneProductReactionLib.reaction_add_H_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_OH_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_O_single_bond_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_add_H_on_the_Os,
    C3OneReactionOneProductReactionLib.reaction_remove_H_on_carbons,
    C3OneReactionOneProductReactionLib.reaction_remove_H_on_Os,
    C3OneReactionOneProductReactionLib.reaction_remove_O_or_OH_on_carbons,
])



#get all teh reactions (specie1, specie2) from All_pathways, this will make it easier to generate mkm file
def reactions_in_pathway(pathway):
	reactions = []
	for index_of_pathways in range(3):
		for index_of_specie in range(len(pathway[index_of_pathways])-1):
			current_reaction = (pathway[index_of_pathways][index_of_specie], pathway[index_of_pathways][index_of_specie+1])
			current_other_speice = get_other_specie(current_reaction)
			current_whole_reaction = (pathway[index_of_pathways][index_of_specie], pathway[index_of_pathways][index_of_specie+1], current_other_speice)
			#avoid duplicates
			if current_reaction in reactions:
				continue
			reactions.append(current_reaction)
	return reactions

#determine the amount of unique reactions in each pathway
with open("All_pathways.json", "r") as file:
	all_pathways = json.load(file)
assert isinstance(all_pathways, dict), "all_pathway file doesn't contain a dictionary"
assert len(all_pathways) == 100, "all_pathway file has pathways less that 100"
length_of_pahtway = {}
for index_of_whole_pathway, key_of_whole_pathway in enumerate(all_pathways):
	length_of_pahtway[index_of_whole_pathway] = len(reactions_in_pathway(all_pathways[key_of_whole_pathway]))

#get the top 10 shortest pathway
sorted_keys = sorted(length_of_pahtway, key=length_of_pahtway.get)[:10]

#generate a small network with these 10 shortest pathways
unique_reactions = set()#to prevent duplicate
for index in sorted_keys:
	unique_reactions.update(reactions_in_pathway(all_pathways[f"{index+1}"]))

#generate network with reaction and barrier energy
network_with_energy = {}
barrier_energy = 0
for reaction in unique_reactions:
	reaction_energy = round(random.uniform(0, 2), 3)
	network_with_energy[f"{reaction}"] = [reaction_energy, barrier_energy]  

#put the final network with energy into json file
with open("network_with_energy.json", "w") as file:
	json.dump(network_with_energy, file, indent=4)



























