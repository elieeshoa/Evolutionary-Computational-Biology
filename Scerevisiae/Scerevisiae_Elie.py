#This python file is a demo on how to load a model and use basic cobra toolbox functions. I hope you find it helpful!
#By Maya Rayle, Dr Ali Zomorrodi Lab Research Intern 2020. mayarayle@college.harvard.edu

#I am doing this in an IDE (An Integrated Development Environment) called PyCharm.
#IDEs make it easier to run and format code, and often catch mistakes
#You can pick a good IDE that works for you, or code directly in the terminal, if you would like to.
#Because there are small formatting differences between IDEs, some of this code might not work for you, but should be easily adaptable.

#Getting Started:

#This line is useful for later functions and needs to be at the top of the document.
from __future__ import division
from os import write
import cobra
import cobra.test
from os.path import join

cobra_config = cobra.Configuration()
data_dir = cobra.test.data_dir
#Import your model
#When you have your model picked, you should save it as an SBML file into the folder where the rest of your python files are
#Then, you can find what the path is to get to it. This is called the absolute path.
#The path is essentially a list of all the folders from broadest to most specific that your file is nested inside
#Think of it like a road map to your file
Scerevisiae_model = cobra.io.read_sbml_model(join(data_dir, '/Users/elieeshoa/Desktop/Elie-Zomorrodi/Scerevisiae/Scerevisiae_iAZ900_noCycles_03_25_2015.xml'))

#Inspect the model: How many reactions? Exchange Reactions? Metabolites? Genes? What are the names of the exchange reactions?
print("Number of Reactions", len(Scerevisiae_model.reactions))
print("Number of Exchange Reactions", len(Scerevisiae_model.exchanges))
print("Number of Metabolites", len(Scerevisiae_model.metabolites))
print("Number of Genes", len(Scerevisiae_model.genes))
# print("Names of Exchange Reactions", Scerevisiae_model.exchanges)
# print("Names of All Reactions", Scerevisiae_model.reactions)


#Setting bounds:
#for a model, it is important to set the bounds on the fluxes before you can ask questions about the model or modify it further
#typically, a reversible reaction will have a lower bound of -1000, unless specified otherwise.
#an irreversible reaction will have a lower bound of 0, unless specified otherwise.
#all reactions have upper bounds of 1000
#exchange reactions should have lower bounds of 0, unless specified otherwise.
    #note that exchange reactions are a type of reversible reaction, but they have a different LB than most other reversible reactions.
    #a negative flux on an exchange reaction means that the microbe is taking in that nutrient from the environment
    #a positive flux on an exchange reaction means that the microbe is secreting that nutrient into the environment
    #if you want a nutrient to be present in the environment and available for the microbe to pick up, set the lower bound to a negative value
        #for instance, if you want to grow your microbe on a glucose medium with oxygen available, set bounds as follows:
        #EX_glc_D_e bounds: [-10, 1000]
        #EX_o2_e bounds: [-2, 1000]
    #if the nutrient is not available in the environment, set the lower bound to 0.

#Setting reversible and irreversible bounds
for k in Scerevisiae_model.reactions:
    if Scerevisiae_model.reactions.get_by_id(k.id).reversibility == True:
        Scerevisiae_model.reactions.get_by_id(k.id).bounds = (-1000, 1000)
    else:
        Scerevisiae_model.reactions.get_by_id(k.id).bounds = (0, 1000)

#Set lower bounds on all exchange reactions to 0.
for c in Scerevisiae_model.exchanges:
    Scerevisiae_model.reactions.get_by_id(c.id).lower_bound = 0

#Setting up the growth medium
#In this example, my growth medium is comprised of the following nutrients
#this means my microbe can take in these nutrients to help it grow
#this could also be referenced in a separate python file, but I have it here so everything is in one place
excess_nutrients = {
#'EX_ca2_e_':[-1000.00000000,1000.00000000],
#'EX_cbl1_e_':[-1000.00000000,1000.00000000],
#'EX_cl_e_':[-1000.00000000,1000.00000000],
#'EX_co2_e_':[-1000.00000000,1000.00000000],
#'EX_cobalt2_e_':[-1000.00000000,1000.00000000],
#'EX_cu2_e_':[-1000.00000000,1000.00000000],
#'EX_fe2_e_':[-1000.00000000,1000.00000000],
#'EX_fe3_e_':[-1000.00000000,1000.00000000],
'EX_h_e_':[-1000.00000000,1000.00000000],
'EX_h2o_e_':[-1000.00000000,1000.00000000],
'EX_k_e_':[-1000.00000000,1000.00000000],
#'EX_mg2_e_':[-1000.00000000,1000.00000000],
#'EX_mn2_e_':[-1000.00000000,1000.00000000],
#'EX_mobd_e_':[-1000.00000000,1000.00000000],
'EX_na1_e_':[-1000.00000000,1000.00000000],
'EX_nh4_e_':[-1000.00000000,1000.00000000],
'EX_pi_e_':[-1000.00000000,1000.00000000],
'EX_so4_e_':[-1000.00000000,1000.00000000],
#'EX_tungs_e_':[-1000.00000000,1000.00000000],
#'EX_zn2_e_':[-1000.00000000,1000.00000000],

# Trace amount of essential nutrients that are present in experimental minimal medium (see iMM904 paper)
'EX_4abz_e_':[-0.5,1000],    # 4-aminobenzoate
'EX_btn_e_':[-0.5,1000],     # biotin
'EX_inost_e_':[-0.5,1000],   # inositol
'EX_nac_e_':[-0.5,1000],     # nicotinate
'EX_pnto_R_e_':[-0.5,1000],  # pantothenate
'EX_thm_e_':[-0.5,1000],     # thiamin
}

regulation = {
}

ngam_atp = {
'ATPM':[1,1],
}

others = {
}



#setting bounds based on growth medium: excess nutrients
for key, value in excess_nutrients.items():
   Scerevisiae_model.reactions.get_by_id(key).bounds = value

#setting bounds based on regulation
# for key, value in regulation.items():
#    Scerevisiae_model.reactions.get_by_id(key).bounds = value

#setting atp bounds
#The reaction below is a sink for ATP (non growth associated ATP). It's ATP for basic life sustaining processes not associated with growth.
#housekeeping.
#GAM = growth associated maintenance ATP
for key, value in ngam_atp.items():
   Scerevisiae_model.reactions.get_by_id(key).bounds = value

#setting bounds based on others
# for key, value in others.items():
#    Scerevisiae_model.reactions.get_by_id(key).bounds = value

#setting the sugar and oxygen
Scerevisiae_model.reactions.get_by_id('EX_glc_e_').bounds = (-10, 1000)
# Aerobic
Scerevisiae_model.reactions.get_by_id('EX_o2_e_').bounds = (-2, 1000)
# Anaerobic
# Scerevisiae_model.reactions.get_by_id('EX_o2_e_').bounds = (0, 1000)
# Scerevisiae_model.reactions.get_by_id('EX_fru_e_').bounds = (-10, 1000)
# Scerevisiae_model.reactions.get_by_id('EX_sucr_e_').bounds = (-10, 1000)

#Now that we've set up the bounds properly, we can begin to analyze the model.

#FLUX BALANCE ANALYSIS
#How to run Flux Balance Analysis (FBA)
#The default is for FBA to optimize biomass
# import time
# start = time.time()
# # fba_solution = Scerevisiae_model.optimize()
# end = time.time()

Scerevisiae_model.objective = 'biomass_wildType'
fba_solution = Scerevisiae_model.optimize()
print("FBA Solution Biomass", fba_solution)
# Scerevisiae_model.optimize()
print('hi')
print(Scerevisiae_model.reactions.get_by_id('EX_glc_e_').bounds) 
print(Scerevisiae_model.reactions.get_by_id('EX_o2_e_').bounds) 
print(Scerevisiae_model.reactions.get_by_id('EX_fru_e_').bounds)
print(Scerevisiae_model.reactions.get_by_id('EX_sucr_e_').bounds)

# for r in Scerevisiae_model.exchanges:
#     print(r.id, Scerevisiae_model.reactions.get_by_id(r.id).flux)
#         # print(key, value)
# print("FBA Solution", fba_solution)
# print("FBA summary", Scerevisiae_model)
# with open("outputs.txt", "w") as f:
#     f. write(str(fba_solution.objective_value))

# print()




# for ex in Scerevisiae_model.exchanges:
#     if Scerevisiae_model.reactions.get_by_id(ex.id).flux > 0:
#         print(ex.id, Scerevisiae_model.reactions.get_by_id(ex.id).flux)


print(" EX_co2_e_ : " + fba_solution.fluxes.get("EX_co2_e_").__str__())
print(" EX_etoh_e_ : " + fba_solution.fluxes.get("EX_etoh_e_").__str__())
print(" EX_for_e_ : " + fba_solution.fluxes.get("EX_for_e_").__str__())
print(" EX_h2o_e_ : " + fba_solution.fluxes.get("EX_h2o_e_").__str__())
print()
# print("Time elapsed to solve the FBA problem: " + str(end - start))
print()

print()


# for x in Scerevisiae_model.exchanges:
#     y = x.id
#     if fba_solution.fluxes.get(key=y):
#         if fba_solution.fluxes.get(key=y) > 0:
#             print('exchange flux metabolites secreted', x.name, fba_solution.fluxes.get(key=y))


# Aerobic
aerobic = ["EX_co2_e_", "EX_etoh_e_", "EX_for_e_", "EX_h2o_e_", "EX_h_e_"]
print()
with open("scerevisiae-aerobic-output.txt", "w") as f:
    print("----- iAZ900 S. cerevisiae under aerobic glucose minimal condition  -----", file=f)
    print("", file=f)
    for x in aerobic:
        print(f"{x} : " + fba_solution.fluxes.get(x).__str__(), file=f)
    print("", file=f)
    # print("Time elapsed to solve the FBA problem: " + str(end - start), file=f)
    print("Objective value: " + str(fba_solution.objective_value), file=f)
    print("", file=f)

    for x in Scerevisiae_model.exchanges:
        y = x.id
        if fba_solution.fluxes.get(key=y):
            if fba_solution.fluxes.get(key=y) > 0:
                print('exchange flux metabolites secreted', x.name, \
                      fba_solution.fluxes.get(key=y), file=f)

# Anaerobic
# anaerobic = ["EX_co2_e_", "EX_etoh_e_"]
# print()
# with open("scerevisiae-anaerobic-output.txt", "w") as f:
#     for x in anaerobic:
#         print(f"{x} : " + fba_solution.fluxes.get(x).__str__(), file=f)
#     print("", file=f)
#     print("Time elapsed to solve the FBA problem: " + str(end - start), file=f)
#     print("Objective value: " + str(fba_solution.objective_value), file=f)
#     print("", file=f)

#     for x in Scerevisiae_model.exchanges:
#         y = x.id
#         if fba_solution.fluxes.get(key=y):
#             if fba_solution.fluxes.get(key=y) > 0:
#                 print('exchange flux metabolites secreted', x.name, \
#                       fba_solution.fluxes.get(key=y), file=f)

#The default is for FBA to optimize biomass. That's why this function will produce the same answer as above
# Scerevisiae_model.objective = 'biomass_wildType'
# fba_solution_biomass = Scerevisiae_model.optimize()
# print("FBA Solution Biomass", fba_solution_biomass)

#If we want to run FBA to optimize a different reaction, let's say we want to see what the minimal glucose uptake is, we can have it optimize for glucose
#note that FBA will maximize the specified reaction.
#Since EX_glc(e) is going to take on a negative flux since it's being taken in by the model, this is asking for the minimum amount of glucose needed for the microbe to live (but not grow).


# Scerevisiae_model.objective = 'EX_glc(e)'
# fba_solution_glc_biomass_growth = Scerevisiae_model.optimize().fluxes.Ec_biomass_iJO1366_core_53p95M
# print("FBA Solution Glucose Biomass Growth", fba_solution_glc_biomass_growth)
# fba_solution_glc_minimal_uptake = Scerevisiae_model.optimize().fluxes.get('EX_glc(e)')
# print("FBA Solution Glucose Minimal Uptake", fba_solution_glc_minimal_uptake)







# #For more on FBA, read the following:
# #https://cobrapy.readthedocs.io/en/latest/simulating.html

# #Let's say I want to find out what the biomass growth of an E coli would be on a sucrose medium instead of a glucose medium:
# #And I want to have a separate model with this condition set
# #How to Copy a Model
# import copy
# Scerevisiae_model_sucrose = copy.deepcopy(Scerevisiae_model)

# #Setting new bounds on sugars and less oxygen, just to see what happens
# Scerevisiae_model_sucrose.reactions.get_by_id('EX_sucr(e)').bounds = (-10, 1000)
# Scerevisiae_model_sucrose.reactions.get_by_id('EX_fru(e)').bounds = (0, 1000)
# Scerevisiae_model_sucrose.reactions.get_by_id('EX_glc(e)').bounds = (0, 1000)
# Scerevisiae_model_sucrose.reactions.get_by_id('EX_o2(e)').bounds = (-2, 1000)

# #new FBA problem
# Scerevisiae_model_sucrose.objective = 'Ec_biomass_iJO1366_core_53p95M'
# print("Sucrose model biomass FBA Solution", Scerevisiae_model_sucrose.optimize())

# #Going back to our original model, let's learn a little more.

# #let's say I want to know what metabolites are secreted out of the model if the model is optimizing growth (as we can assume a microbe will do)
# #I'll use a for loop to run through all the exchange reactions, and print any with a positive flux value
# #Remember that a positive flux means secretion
# for x in Scerevisiae_model.exchanges:
#     y = x.id
#     if fba_solution_biomass.fluxes.get(key=y):
#         if fba_solution_biomass.fluxes.get(key=y) > 0:
#             print('exchange flux metabolites secreted', x.name, fba_solution_biomass.fluxes.get(key=y))

# #How to edit a reaction
# #If I want to add metabolites to a reaction, I can go in and edit the reaction
# #In this example, I will be adding ATP to the D-amino acid dehydrogenase reaction
# DAAD_rxn = Scerevisiae_model.reactions.DAAD
# print("Original DAAD Reaction", DAAD_rxn)
# print(Scerevisiae_model.reactions.DAAD.check_mass_balance())
# DAAD_rxn.add_metabolites({Scerevisiae_model.metabolites.get_by_id("atp_c"): -1,
#     Scerevisiae_model.metabolites.get_by_id("adp_c"): 1,
#     Scerevisiae_model.metabolites.get_by_id("pi_c"): 1,
#     Scerevisiae_model.metabolites.get_by_id("h_e"): 1
#     })
# print("Modified DAAD Reaction", DAAD_rxn)
# print(Scerevisiae_model.reactions.DAAD.check_mass_balance())

# #Now, I can take away the metabolites.
# DAAD_rxn.subtract_metabolites({Scerevisiae_model.metabolites.get_by_id("atp_c"): -1,
#     Scerevisiae_model.metabolites.get_by_id("adp_c"): 1,
#     Scerevisiae_model.metabolites.get_by_id("pi_c"): 1,
#     Scerevisiae_model.metabolites.get_by_id("h_e"): 1
#     })
# print("Back to Original DAAD Reaction", DAAD_rxn)
# print(Scerevisiae_model.reactions.DAAD.check_mass_balance())

# #Gene Knockouts
# # Knock out all reactions in the model and report the reaction deletions that lead to a zero max biomass flux
# # (these are called essential reactions).
# # A for loop to run through all of the reactions deleting one at a time and extract the solution. Print if it is 0.
# Scerevisiae_model.objective = 'Ec_biomass_iJO1366_core_53p95M'
# fba_solution_biomass = Scerevisiae_model.optimize()
# print('complete model: ', fba_solution_biomass)
# for m in Scerevisiae_model.reactions:
#      name_rxn = m.id
#      with Scerevisiae_model:
#          Scerevisiae_model.reactions.get_by_id(name_rxn).knock_out()
#          optimize_solution_knockout = Scerevisiae_model.optimize()
#          if optimize_solution_knockout.objective_value < .001 and optimize_solution_knockout.objective_value > -.001:
#             print('Model Solutions', name_rxn, Scerevisiae_model.optimize())

# #I hope this has been a useful demo! If you are looking for more information, just look at the following website:
# #https://cobrapy.readthedocs.io/en/latest/
