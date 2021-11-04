#This python file is a demo on how to load a model and use basic cobra toolbox functions. I hope you find it helpful!
#By Maya Rayle, Dr Ali Zomorrodi Lab Research Intern 2020. mayarayle@college.harvard.edu

#How To Connect to ERISone?
# Setup a VPN.
# Use Cisco AnyConnect Secure Mobility Client.
# Follow these instructions to download it: https://www.partners.org/vpn/
# After you’ve downloaded it, when it says VPN: ready to connect, enter “pvc.partners.org/saml” in the box.
# Enter your MGH credentials afterwards when prompted.
# You will have to re-enter every 12 hours, because it logs you out automatically.

# Connecting to ERISone Linux Cluster:
# Make sure that your VPN is activated.
# Follow the instructions here: https://rc.partners.org/kb/article/2814
# Make sure that you enter your partners account in lower case letters. For example: mmm11@erisone.partners.org
# For the password, enter your regular MGH login password.

#Before you start coding, you will need to install cobra through your terminal.
#If you have a Mac, do the following:
#open terminal
#type "ls" to get a list of your folders in your library
#type "cd [name of the folder with python files in it]"
#type "source venv/bin/activate" to access the virtual environment so you can download packages
#type "pip install cobra"
#If you have a Windows operating system or something else, follow instructions here:
#https://github.com/opencobra/cobrapy/blob/devel/INSTALL.rst

#I am doing this in an IDE (An Integrated Development Environment) called PyCharm.
#IDEs make it easier to run and format code, and often catch mistakes
#You can pick a good IDE that works for you, or code directly in the terminal, if you would like to.
#Because there are small formatting differences between IDEs, some of this code might not work for you, but should be easily adaptable.

#Getting Started:

#This line is useful for later functions and needs to be at the top of the document.
from __future__ import division
import cobra
cobra_config = cobra.Configuration()
import cobra.test
from os.path import join
data_dir = cobra.test.data_dir

# Set your output file and your objective function
output_file = "Sc_out.txt"
objective = "biomass_core"
#Import your model
#When you have your model picked, you should save it as an SBML file into the folder where the rest of your python files are
#Then, you can find what the path is to get to it. This is called the absolute path.
#The path is essentially a list of all the folders from broadest to most specific that your file is nested inside
#Think of it like a road map to your file

# The path for to your .xml SBML model, e.g.
# /Users/elieeshoa/Desktop/Elie-Zomorrodi/Scerevisiae/Scerevisiae_iAZ900.xml 
path_to_xml_model = "/Users/elieeshoa/Desktop/Elie-Zomorrodi/Scerevisiae/" + \
                    "Scerevisiae_iAZ900_noCycles_03_25_2015.xml"

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

# Setting reversible and irreversible bounds
for k in met_model.reactions:
    if met_model.reactions.get_by_id(k.id).reversibility == True:
        met_model.reactions.get_by_id(k.id).bounds = (-1000, 1000)
    else:
        met_model.reactions.get_by_id(k.id).bounds = (0, 1000)

# Set lower bounds on all exchange reactions to 0.
for c in met_model.exchanges:
    met_model.reactions.get_by_id(c.id).lower_bound = 0

#Setting up the growth medium
#In this example, my growth medium is comprised of the following nutrients
#this means my microbe can take in these nutrients to help it grow
#this could also be referenced in a separate python file, but I have it here so everything is in one place

# Import the .py file (in this example it is iAZ_minimal.py) that has any 
# subcollection of these exact variables: `excess_nutrients` & `ngam_atp` 
# & `regulation` & `others`
# Make sure this file is in the same directory you're working in
import iAZ900_minimal as growth_medium

rxn = met_model.reactions.get_by_id("ZYMSTESTH_SC")
# for rxn in model.reactions:
#     coefficients[rxn.forward_variable] = 1.
#     coefficients[rxn.reverse_variable] = 1.
print()
print()
rxn.objective_coefficient = 1.
print(rxn.objective_coefficient)
print()
print()

#setting bounds based on excess nutrients
try:
    for key, value in growth_medium.excess_nutrients.items():
        met_model.reactions.get_by_id(key).bounds = value
except:
    pass        

#ATPM bounds.
#met_model.reactions.ATPM.bounds = (1,1)
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

# the aerobic condition.
met_model.reactions.EX_o2_b.bounds = (-2, 1000)
# # the anaerobic condition.
# met_model.reactions.EX_o2_b.bounds = (0, 1000)

#Set all "_e_" reactions to -1000, 1000
for x in met_model.reactions:
    if "EX_" in x.id and "_e_" in x.id:
        met_model.reactions.get_by_id(x.id).bounds = (-1000, 1000)

#Setting Sugar Bounds
met_model.reactions.EX_glc_D_b.bounds = (0, 1000)
met_model.reactions.EX_fru_b.bounds = (0, 1000)
met_model.reactions.EX_sucr_b.bounds = (-10, 1000)

#Now that we've set up the bounds properly, we can begin to analyze the model.

#FLUX BALANCE ANALYSIS
#How to run Flux Balance Analysis (FBA)
#The default is for FBA to optimize biomass
try:
    fba_solution = met_model.optimize()
except:
    print('model is infeasible')

print("FBA Solution", fba_solution)
print("FBA summary", met_model.summary)

#The default is for FBA to optimize biomass. That's why this function will produce the same answer as above
met_model.objective = objective
try:
    fba_solution_biomass = met_model.optimize()
except:
    print('model is infeasible')
print("FBA Solution Biomass", fba_solution_biomass)

with open(output_file, "a+") as f:
    print("FBA Solution Biomass", fba_solution_biomass.objective_value, file=f)

rxn = met_model.reactions.get_by_id("ZYMSTESTH_SC")
# for rxn in model.reactions:
#     coefficients[rxn.forward_variable] = 1.
#     coefficients[rxn.reverse_variable] = 1.
print()
print()
print(rxn.objective_coefficient)
print()
print()

#If we want to run FBA to optimize a different reaction, let's say we want to see what the minimal glucose uptake is, we can have it optimize for glucose
#note that FBA will maximize the specified reaction.
#Since EX_glc_D_b is going to take on a negative flux since it's being taken in by the model, this is asking for the minimum amount of glucose needed for the microbe to live (but not grow).
met_model.objective = 'EX_glc_D_b'
fba_solution_glc_biomass_growth = met_model.optimize().fluxes.biomass_core
print("FBA Solution Glucose Biomass Growth", fba_solution_glc_biomass_growth)
fba_solution_glc_minimal_uptake = met_model.optimize().fluxes.get('EX_glc_D_b')
print("FBA Solution Glucose Minimal Uptake", fba_solution_glc_minimal_uptake)

#For more on FBA, read the following:
#https://cobrapy.readthedocs.io/en/latest/simulating.html

#Let's say I want to find out what the biomass growth of an E coli would be on a sucrose medium instead of a glucose medium:
#And I want to have a separate model with this condition set
#How to Copy a Model
import copy
met_model_sucrose = copy.deepcopy(met_model)

#Setting new bounds on sugars and less oxygen, just to see what happens
met_model_sucrose.reactions.get_by_id('EX_sucr_b').bounds = (-10, 1000)
met_model_sucrose.reactions.get_by_id('EX_fru_b').bounds = (0, 1000)
met_model_sucrose.reactions.get_by_id('EX_glc_D_b').bounds = (0, 1000)
met_model_sucrose.reactions.get_by_id('EX_o2_b').bounds = (-2, 1000)

#new FBA problem
met_model_sucrose.objective = objective
print("Sucrose model biomass FBA Solution", met_model_sucrose.optimize())

#Going back to our original model, let's learn a little more. 

#let's say I want to know what metabolites are secreted out of the model if the model is optimizing growth (as we can assume a microbe will do)
#I'll use a for loop to run through all the exchange reactions, and print any with a positive flux value
#Remember that a positive flux means secretion
for x in met_model.exchanges:
    y = x.id
    if fba_solution_biomass.fluxes.get(key=y):
        if fba_solution_biomass.fluxes.get(key=y) > 0:
            print('exchange flux metabolites secreted', x.name, fba_solution_biomass.fluxes.get(key=y))
            with open(output_file, "a+") as f:
                print('exchange flux metabolites secreted', x.name, fba_solution_biomass.fluxes.get(key=y), file=f)

# with open("Sc_out_aerobic.txt", "w") as f:
#     print("FBA Solution: \n", file=f)
#     print(fba_solution.fluxes.biomass_core, file=f)
#     print("\n\n\n FBA summary: \n", file=f)
#     print(met_model.summary.__str__(), file=f)

#How to edit a reaction
#If I want to add metabolites to a reaction, I can go in and edit the reaction
#In this example, I will be adding ATP to the EX_2mbtoh_b reaction
EX_2mbtoh_b_rxn = met_model.reactions.EX_2mbtoh_b
print("Original EX_2mbtoh_b Reaction", EX_2mbtoh_b_rxn)
print(met_model.reactions.EX_2mbtoh_b.check_mass_balance())
EX_2mbtoh_b_rxn.add_metabolites({met_model.metabolites.get_by_id("atp_c"): -1,
    met_model.metabolites.get_by_id("adp_c"): 1,
    met_model.metabolites.get_by_id("pi_c"): 1,
    met_model.metabolites.get_by_id("h_e"): 1
    })
print("Modified EX_2mbtoh_b Reaction", EX_2mbtoh_b_rxn)
print(met_model.reactions.EX_2mbtoh_b.check_mass_balance())

#Now, I can take away the metabolites.
EX_2mbtoh_b_rxn.subtract_metabolites({met_model.metabolites.get_by_id("atp_c"): -1,
    met_model.metabolites.get_by_id("adp_c"): 1,
    met_model.metabolites.get_by_id("pi_c"): 1,
    met_model.metabolites.get_by_id("h_e"): 1
    })
print("Back to Original EX_2mbtoh_b Reaction", EX_2mbtoh_b_rxn)
print(met_model.reactions.EX_2mbtoh_b.check_mass_balance())

#Gene Knockouts
# Knock out all reactions in the model and report the reaction deletions that lead to a zero max biomass flux
# (these are called essential reactions).
# A for loop to run through all of the reactions deleting one at a time and extract the solution. Print if it is 0.
met_model.objective = objective
fba_solution_biomass = met_model.optimize()
print('complete model: ', fba_solution_biomass)
for m in met_model.reactions:
     name_rxn = m.id
     with met_model:
         met_model.reactions.get_by_id(name_rxn).knock_out()
         optimize_solution_knockout = met_model.optimize()
         if optimize_solution_knockout.objective_value < .001 and optimize_solution_knockout.objective_value > -.001:
            print('Model Solutions', name_rxn, met_model.optimize())

#I hope this has been a useful demo! If you are looking for more information, just look at the following website:
#https://cobrapy.readthedocs.io/en/latest/
