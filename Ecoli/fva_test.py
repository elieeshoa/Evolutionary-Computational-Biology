from __future__ import division
import ast
from os import write
from cobra.flux_analysis import flux_variability_analysis
import cobra
cobra_config = cobra.Configuration()
import cobra.test
from os.path import join
data_dir = cobra.test.data_dir

# The path for to your .xml SBML model, e.g.
# /Users/elieeshoa/Desktop/Elie-Zomorrodi/Scerevisiae/Scerevisiae_iAZ900.xml 
path_to_xml_model = "/Users/elieeshoa/Desktop/Elie-Zomorrodi/Ecoli/Ecoli_iJO1366_updated.xml"

# Define the organism's metabolic model
met_model = cobra.io.read_sbml_model(join(data_dir, path_to_xml_model))

# Inspect the model: How many reactions? Exchange Reactions? Metabolites? 
# Genes? What are the names of the exchange reactions?
print("Number of Reactions", len(met_model.reactions))
print("Number of Exchange Reactions", len(met_model.exchanges))
print("Number of Metabolites", len(met_model.metabolites))
print("Number of Genes", len(met_model.genes))
print("Names of Exchange Reactions", met_model.exchanges)
print("Names of All Reactions", met_model.reactions)

# Setting reversible and irreversible bounds
for k in met_model.reactions:
    if met_model.reactions.get_by_id(k.id).reversibility == True:
        met_model.reactions.get_by_id(k.id).bounds = (-1000, 1000)
    else:
        met_model.reactions.get_by_id(k.id).bounds = (0, 1000)

# Set lower bounds on all exchange reactions to 0.
for c in met_model.exchanges:
    met_model.reactions.get_by_id(c.id).lower_bound = 0

# Import the .py file (in this example it is iAZ_minimal.py) that has any 
# subcollection of these exact variables: `excess_nutrients` & `ngam_atp` 
# & `regulation` & `others`
# Make sure this file is in the same directory you're working in
import iJO1366_minimal_glucose_aerobic_edited as growth_medium

#setting bounds based on excess nutrients
try:
    for key, value in growth_medium.excess_nutrients.items():
        met_model.reactions.get_by_id(key).bounds = value
except:
    pass        

#ATPM bounds.
try:
    for key, value in growth_medium.ngam_atp.items():
        met_model.reactions.get_by_id(key).bounds = value
except:
    pass

# Setting bounds based on others
try:
    for key, value in growth_medium.others.items():
        met_model.reactions.get_by_id(key).bounds = value
except:
    pass

# Setting bounds based on regulation
try:
    for key, value in growth_medium.regulation.items():
        met_model.reactions.get_by_id(key).bounds = value
except:
    pass


#setting the sugar and oxygen
met_model.reactions.get_by_id('EX_glc(e)').bounds = (-10, 1000)
# Aerobic
met_model.reactions.get_by_id('EX_o2(e)').bounds = (-20, 1000)
# Anaerobic
# Don't set the o2
# print(met_model.reactions.get_by_id('EX_fru(e)').bounds)
# met_model.reactions.get_by_id('EX_fru(e)').bounds = (0, 1000)
# met_model.reactions.get_by_id('EX_sucr(e)').bounds = (0, 1000)

if __name__ == '__main__':  
    # met_model.optimize()
    # # print("FVA summary: ", met_model.summary(fva=0.95))
    # summary = met_model.summary(fva=0.95)
    # print("hi")
    # with open("fva_html", "w") as f:
    #     f.write(summary._repr_html_())
    # # summary.to_string()
    # with open("fva_out.txt", "w") as f:
    #     f.write(summary.to_string())











    # dict = {}
    # for rxn in cobra.flux_analysis.flux_variability_analysis(met_model, met_model.reactions, fraction_of_optimum=0.95):
    #         dict[rxn.id] = (rxn.minimum, rxn.maximum)
    # print(dict)






    # final = cobra.flux_analysis.flux_variability_analysis(met_model, met_model.reactions, fraction_of_optimum=0.95)
    # # # final.to_dict()
    # dict = {}
    # final = final.to_dict('index')
    # # print(final)
    # for key, value in final.items(): 
    #     dict[key] = [value["minimum"], value["maximum"]]
    # print(dict)
    # with open("fva_dict.txt", "a+") as f:
    #     print(dict, file=f)
    #     f.write("{\n")
    #     for key in dict.keys(): 
    #         f.write(f"{k}:{dict[k]}\n")
    #     f.write("}")
    #     #     f.write('%s:%s\n' % (key, value))


    # Turn into desired dictionary format
    with open("fva_dict.txt", "r") as f:
        content = f.read()
        d = ast.literal_eval(content)

    with open("fva_dict1.txt", "a+") as f:   
        f.write("rxns_flux_ranges = {\n")
        for item in d.items(): 
            f.write(f"'{str(item[0])}' : {str(item[1])}")
            f.write("\n")
        f.write("}")


        # f.write(final.to_dict('index').__str__)


    # with open('fva_out.txt', 'wb') as file:
    #     pickle.dump(summary, file)
    # with open("fva_out.txt", "w") as f:
    #     write(f, met_model.summary(fva=0.95))
    # flux_variability_analysis(met_model, met_model.reactions[:10])