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
from os import write

#Importing Cobra and preparing to import your model
import cobra
cobra_config = cobra.Configuration()
import cobra.test
from os.path import join
data_dir = cobra.test.data_dir
#Import your model
#When you have your model picked, you should save it as an SBML file into the folder where the rest of your python files are
#Then, you can find what the path is to get to it. This is called the absolute path.
#The path is essentially a list of all the folders from broadest to most specific that your file is nested inside
#Think of it like a road map to your file
E_coli_model = cobra.io.read_sbml_model(join(data_dir, '/Users/elieeshoa/Desktop/Elie-Zomorrodi/Ecoli/Ecoli_iJO1366_updated.xml'))




#Inspect the model: How many reactions? Exchange Reactions? Metabolites? Genes? What are the names of the exchange reactions?
print("Number of Reactions", len(E_coli_model.reactions))
print("Number of Exchange Reactions", len(E_coli_model.exchanges))
print("Number of Metabolites", len(E_coli_model.metabolites))
print("Number of Genes", len(E_coli_model.genes))
# print("Names of Exchange Reactions", E_coli_model.exchanges)
# print("Names of All Reactions", E_coli_model.reactions)


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
for k in E_coli_model.reactions:
    if E_coli_model.reactions.get_by_id(k.id).reversibility == True:
        E_coli_model.reactions.get_by_id(k.id).bounds = (-1000, 1000)
    else:
        E_coli_model.reactions.get_by_id(k.id).bounds = (0, 1000)

#Set lower bounds on all exchange reactions to 0.
for c in E_coli_model.exchanges:
    E_coli_model.reactions.get_by_id(c.id).lower_bound = 0

#Setting up the growth medium
#In this example, my growth medium is comprised of the following nutrients
#this means my microbe can take in these nutrients to help it grow
#this could also be referenced in a separate python file, but I have it here so everything is in one place


# Aerobic
excess_nutrients = {
'EX_ca2(e)':[-1000,1000],
'EX_cbl1(e)':[-1000,1000],
'EX_cl(e)':[-1000,1000],
'EX_co2(e)':[-1000,1000],
'EX_cobalt2(e)':[-1000,1000],
'EX_cu2(e)':[-1000,1000],
'EX_fe2(e)':[-1000,1000],
'EX_fe3(e)':[-1000,1000],
'EX_h(e)':[-1000,1000],
'EX_h2o(e)':[-1000,1000],
'EX_k(e)':[-1000,1000],
'EX_mg2(e)':[-1000,1000],
'EX_mn2(e)':[-1000,1000],
'EX_mobd(e)':[-1000,1000],
'EX_na1(e)':[-1000,1000],
'EX_nh4(e)':[-1000,1000],
'EX_ni2(e)':[-1000,1000],
'EX_pi(e)':[-1000,1000],
'EX_sel(e)':[-1000,1000],
'EX_slnt(e)':[-1000,1000],
'EX_so4(e)':[-1000,1000],
'EX_tungs(e)':[-1000,1000],
'EX_zn2(e)':[-1000,1000],
}

regulation = {
'14GLUCANabcpp':[0,0],
'14GLUCANtexi':[0,0],
'ACOAD1f':[0,0],
'ACOAD2f':[0,0],
'ACOAD3f':[0,0],
'ACOAD4f':[0,0],
'ACOAD5f':[0,0],
'ACOAD6f':[0,0],
'ACOAD7f':[0,0],
'ACOAD8f':[0,0],
'ACS':[0,0],
'ALDD2y':[0,0],
'ALDD3y':[0,0],
'ARAI':[0,0],
'ARBabcpp':[0,0],
'ARBt2rpp':[0,0],
'ASNNpp':[0,0],
'ASPT':[0,0],
'ASPt2_2pp':[0,0],
'CADVtpp':[0,0],
'CRNBTCT':[0,0],
'CRNCAR':[0,0],
'CRNCBCT':[0,0],
'CRNCDH':[0,0],
'CRNt7pp':[0,0],
'CRNt8pp':[0,0],
'CYANST':[0,0],
'CYSDS':[0,0],
'CYTD':[0,0],
'DAAD':[0,0],
'DCYTD':[0,0],
'DMSOR1':[0,0],
'DMSOR2':[0,0],
'DOGULNR':[0,0],
'DXYLK':[0,0],
'ECOAH1':[0,0],
'ECOAH2':[0,0],
'ECOAH3':[0,0],
'ECOAH4':[0,0],
'ECOAH5':[0,0],
'ECOAH6':[0,0],
'ECOAH7':[0,0],
'ECOAH8':[0,0],
'F6Pt6_2pp':[0,0],
'FCI':[0,0],
'FCLK':[0,0],
'FORt2pp':[0,0],
'FORtppi':[0,0],
'FRD2':[0,0],
'FRD3':[0,0],
'FRUURt2rpp':[0,0],
'FUCtpp':[0,0],
'FUMt2_2pp':[0,0],
'G3PCabcpp':[0,0],
'G3PD5':[0,0],
'G3PD6':[0,0],
'G3PD7':[0,0],
'G3PEabcpp':[0,0],
'G3PGabcpp':[0,0],
'G3PIabcpp':[0,0],
'G3PSabcpp':[0,0],
'G6Pt6_2pp':[0,0],
'GALabcpp':[0,0],
'GALS3':[0,0],
'GAM6Pt6_2pp':[0,0],
'GLCabcpp':[0,0],
'GLCtexi':[0,0],
'GLUNpp':[0,0],
'GLYC3Pabcpp':[0,0],
'GLYC3Pt6pp':[0,0],
'GLYK':[0,0],
'GPDDA1pp':[0,0],
'GPDDA2pp':[0,0],
'GPDDA3pp':[0,0],
'GPDDA4pp':[0,0],
'GPDDA5pp':[0,0],
'HACD1':[0,0], # Name was HACD1i in iAF1260 
'HACD2':[0,0], # Name was HACD2i in iAF1260
'HACD3':[0,0], # Name was HACD3i in iAF1260
'HACD4':[0,0], # Name was HACD4i in iAF1260
'HACD5':[0,0], # Name was HACD5i in iAF1260
'HACD6':[0,0], # Name was HACD6i in iAF1260
'HACD7':[0,0], # Name was HACD7i in iAF1260
'HACD8':[0,0], # Name was HACD8i in iAF1260
'HYD1pp':[0,0],
'HYD2pp':[0,0],
'HYD3pp':[0,0],
'ICL':[0,0],
#'KAT1':[0,0],  # Not in the iJO1366 model
#'KAT2':[0,0],  # Not in the iJO1366 model
#'KAT3':[0,0],  # Not in the iJO1366 model
#'KAT4':[0,0],  # Not in the iJO1366 model
#'KAT5':[0,0],  # Not in the iJO1366 model
#'KAT6':[0,0],  # Not in the iJO1366 model
#'KAT7':[0,0],  # Not in the iJO1366 model
#'KAT8':[0,0],  # Not in the iJO1366 model
'LACZ':[0,0],
'LCARS':[0,0],
'LCTStpp':[0,0],
'LYXI':[0,0],
'LYXt2pp':[0,0],
'MALDt2_2pp':[0,0],
'MALS':[0,0],
'MALt2_2pp':[0,0],
'MALTabcpp':[0,0],
'MALTHXabcpp':[0,0],
'MALTHXtexi':[0,0],
'MALTPTabcpp':[0,0],
'MALTptspp':[0,0],
'MALTPTtexi':[0,0],
'MALTtexi':[0,0],
'MALTTRabcpp':[0,0],
'MALTTRtexi':[0,0],
'MALTTTRabcpp':[0,0],
'MALTTTRtexi':[0,0],
'MAN6Pt6_2pp':[0,0],
'MELIBt2pp':[0,0],
'NO2t2rpp':[0,0],
'NTRIR2x':[0,0],
'OBTFL':[0,0],
'OROTt2_2pp':[0,0],
'P5CD':[0,0],
'PFL':[0,0],
'PPAt4pp':[0,0],
'PROD2':[0,0],
'PROt4pp':[0,0],
'RBK':[0,0],
'RBK_L1':[0,0],
'RBP4E':[0,0],
'RMI':[0,0],
'RMK':[0,0],
'RMNtpp':[0,0],
'RMPA':[0,0],
'RNTR1c2':[0,0],  # Name was RNTR1c in iAF1260
'RNTR2c2':[0,0],  # Name was RNTR2c in iAF1260
'RNTR3c2':[0,0],  # Name was RNTR3c in iAF1260
'RNTR4c2':[0,0],  # Name was RNTR4c in iAF1260
'SBTPD':[0,0],
'SBTptspp':[0,0],
'SERD_D':[0,0],
'SUCCt2_2pp':[0,0],
'TMAOR1':[0,0],
'TMAOR2':[0,0],
'TMDPP':[0,0],
'TRE6PH':[0,0],
'TREptspp':[0,0],
'TRPAS2':[0,0],
#'TRPt2rpp':[0,0],
'XYLabcpp':[0,0],
'XYLI1':[0,0],
'XYLI2':[0,0],
'XYLK':[0,0],
}

ngam_atp = {
'ATPM':[8.39,8.39],
}

others = {
'ALAt2pp':[0,0],   # Creates a thermodynamically infeasible loop with ALAt2rpp 
'GLYt2pp':[0,0]    # Creates a thermodynamically infeasible loop with GLYt2rpp
}



# Anaerobic
# excess_nutrients = {
# 'EX_ca2(e)':[-1000,1000],
# 'EX_cbl1(e)':[-1000,1000],
# 'EX_cl(e)':[-1000,1000],
# 'EX_co2(e)':[-1000,1000],
# 'EX_cobalt2(e)':[-1000,1000],
# 'EX_cu2(e)':[-1000,1000],
# 'EX_fe2(e)':[-1000,1000],
# 'EX_fe3(e)':[-1000,1000],
# 'EX_h(e)':[-1000,1000],
# 'EX_h2o(e)':[-1000,1000],
# 'EX_k(e)':[-1000,1000],
# 'EX_mg2(e)':[-1000,1000],
# 'EX_mn2(e)':[-1000,1000],
# 'EX_mobd(e)':[-1000,1000],
# 'EX_na1(e)':[-1000,1000],
# 'EX_nh4(e)':[-1000,1000],
# 'EX_ni2(e)':[-1000,1000],
# 'EX_pi(e)':[-1000,1000],
# 'EX_sel(e)':[-1000,1000],
# 'EX_slnt(e)':[-1000,1000],
# 'EX_so4(e)':[-1000,1000],
# 'EX_tungs(e)':[-1000,1000],
# 'EX_zn2(e)':[-1000,1000],
# }

# regulation = {
# '14GLUCANabcpp':[0,0],
# '14GLUCANtexi':[0,0],
# 'ALDD2y':[0,0],
# 'ALDD3y':[0,0],
# 'ARAI':[0,0],
# 'ARBabcpp':[0,0],
# 'ARBt2rpp':[0,0],
# 'CRNBTCT':[0,0],
# 'CRNCAR':[0,0],
# 'CRNCBCT':[0,0],
# 'CRNCDH':[0,0],
# 'CRNt7pp':[0,0],
# 'CRNt8pp':[0,0],
# 'CYANST':[0,0],
# 'CYSDS':[0,0],
# 'CYTD':[0,0],
# 'DAAD':[0,0],
# 'DCYTD':[0,0],
# 'DOGULNR':[0,0],
# 'DXYLK':[0,0],
# 'ECOAH1':[0,0],
# 'ECOAH2':[0,0],
# 'ECOAH3':[0,0],
# 'ECOAH4':[0,0],
# 'ECOAH5':[0,0],
# 'ECOAH6':[0,0],
# 'ECOAH7':[0,0],
# 'ECOAH8':[0,0],
# 'F6Pt6_2pp':[0,0],
# 'FCI':[0,0],
# 'FCLK':[0,0],
# 'FRUURt2rpp':[0,0],
# 'FUCtpp':[0,0],
# 'G3PCabcpp':[0,0],
# 'G3PEabcpp':[0,0],
# 'G3PGabcpp':[0,0],
# 'G3PIabcpp':[0,0],
# 'G3PSabcpp':[0,0],
# 'G6Pt6_2pp':[0,0],
# 'GALabcpp':[0,0],
# 'GALS3':[0,0],
# 'GAM6Pt6_2pp':[0,0],
# 'GLCabcpp':[0,0],
# 'GLCtexi':[0,0],
# 'GLYC3Pabcpp':[0,0],
# 'GLYC3Pt6pp':[0,0],
# 'GLYK':[0,0],
# 'HACD1':[0,0], # Name was HACD1i in iAF1260 
# 'HACD2':[0,0], # Name was HACD2i in iAF1260
# 'HACD3':[0,0], # Name was HACD3i in iAF1260
# 'HACD4':[0,0], # Name was HACD4i in iAF1260
# 'HACD5':[0,0], # Name was HACD5i in iAF1260
# 'HACD6':[0,0], # Name was HACD6i in iAF1260
# 'HACD7':[0,0], # Name was HACD7i in iAF1260
# 'HACD8':[0,0], # Name was HACD8i in iAF1260
# #'KAT1':[0,0],  # Not in iJO1366 model
# #'KAT2':[0,0],  # Not in iJO1366 model
# #'KAT3':[0,0],  # Not in iJO1366 model
# #'KAT4':[0,0],  # Not in iJO1366 model
# #'KAT5':[0,0],  # Not in iJO1366 model
# #'KAT6':[0,0],  # Not in iJO1366 model
# #'KAT7':[0,0],  # Not in iJO1366 model
# #'KAT8':[0,0],  # Not in iJO1366 model
# 'LACZ':[0,0],
# 'LCARS':[0,0],
# 'LCTStpp':[0,0],
# 'LYXI':[0,0],
# 'LYXt2pp':[0,0],
# 'MALTabcpp':[0,0],
# 'MALTHXabcpp':[0,0],
# 'MALTHXtexi':[0,0],
# 'MALTPTabcpp':[0,0],
# 'MALTptspp':[0,0],
# 'MALTPTtexi':[0,0],
# 'MALTtexi':[0,0],
# 'MALTTRabcpp':[0,0],
# 'MALTTRtexi':[0,0],
# 'MALTTTRabcpp':[0,0],
# 'MALTTTRtexi':[0,0],
# 'MAN6Pt6_2pp':[0,0],
# 'MELIBt2pp':[0,0],
# 'P5CD':[0,0],
# 'PPAt4pp':[0,0],
# 'PROD2':[0,0],
# 'PROt4pp':[0,0],
# 'RBK':[0,0],
# 'RBK_L1':[0,0],
# 'RBP4E':[0,0],
# 'RMI':[0,0],
# 'RMK':[0,0],
# 'RMNtpp':[0,0],
# 'RMPA':[0,0],
# 'SBTPD':[0,0],
# 'SBTptspp':[0,0],
# 'SERD_D':[0,0],
# 'TMDPP':[0,0],
# 'TRE6PH':[0,0],
# 'TREptspp':[0,0],
# 'TRPAS2':[0,0],
# #'TRPt2rpp':[0,0],
# 'XYLabcpp':[0,0],
# 'XYLI1':[0,0],
# 'XYLI2':[0,0],
# 'XYLK':[0,0],
# }

# ngam_atp = {
# 'ATPM':[8.39000000,8.39000000],
# }

# others = {
# 'ALAt2pp':[0,0],   # Creates a thermodynamically infeasible loop with ALAt2rpp 
# 'GLYt2pp':[0,0]    # Creates a thermodynamically infeasible loop with GLYt2rpp
# }

#setting bounds based on growth medium: excess nutrients
for key, value in excess_nutrients.items():
   E_coli_model.reactions.get_by_id(key).bounds = value

#setting bounds based on regulation
for key, value in regulation.items():
   E_coli_model.reactions.get_by_id(key).bounds = value

#setting atp bounds
#The reaction below is a sink for ATP (non growth associated ATP). It's ATP for basic life sustaining processes not associated with growth.
#housekeeping.
#GAM = growth associated maintenance ATP
for key, value in ngam_atp.items():
   E_coli_model.reactions.get_by_id(key).bounds = value

#setting bounds based on others
for key, value in others.items():
   E_coli_model.reactions.get_by_id(key).bounds = value

#setting the sugar and oxygen
E_coli_model.reactions.get_by_id('EX_glc(e)').bounds = (-10, 1000)
# Aerobic
E_coli_model.reactions.get_by_id('EX_o2(e)').bounds = (-20, 1000)
# Anaerobic
# Don't set the o2
E_coli_model.reactions.get_by_id('EX_fru(e)').bounds = (0, 1000)
E_coli_model.reactions.get_by_id('EX_sucr(e)').bounds = (0, 1000)

#Now that we've set up the bounds properly, we can begin to analyze the model.

# FVA
cobra.flux_analysis.flux_variability_analysis(E_coli_model, E_coli_model.reactions[:10], fraction_of_optimum=0.9)
fva_solution = E_coli_model.optimize()
# print(E_coli_model.summary(fva=0.95))
print("FVA Solution", fva_solution)
print("FVA summary", E_coli_model.summary(fva=0.95))

with open("fva_out.txt", "w") as f:
    write(f, E_coli_model.summary(fva=0.95))




# #FLUX BALANCE ANALYSIS
# #How to run Flux Balance Analysis (FBA)
# #The default is for FBA to optimize biomass
# import time
# start = time.time()
# fba_solution = E_coli_model.optimize()
# end = time.time()
# print("FBA Solution", fba_solution)
# print("FBA summary", E_coli_model.summary)
# # with open("ecoli-aerobic-output.txt", "w") as f:
# #     f. write(str(fba_solution.objective_value) + " EX_5mtr(e) : " + fba_solution.fluxes.get("EX_5mtr(e)").__str__())

# # Aerobic
# aerobic = ["EX_5mtr(e)", "EX_co2(e)", "EX_h(e)", "EX_h2o(e)"]
# print()
# with open("ecoli-aerobic-output.txt", "w") as f:
#     for x in aerobic:
#         print(f"{x} : " + fba_solution.fluxes.get(x).__str__(), file=f)
#     print("", file=f)
#     print("Time elapsed to solve the FBA problem: " + str(end - start), file=f)
#     print("Objective value: " + str(fba_solution.objective_value), file=f)
#     print("", file=f)

#     for x in E_coli_model.exchanges:
#         y = x.id
#         if fba_solution.fluxes.get(key=y):
#             if fba_solution.fluxes.get(key=y) > 0:
#                 print('exchange flux metabolites secreted', x.name, \
#                       fba_solution.fluxes.get(key=y), file=f)

# # Anaerobic
# # anaerobic = ["EX_5mtr(e)", "EX_ac(e)", "EX_co2(e)", "EX_glyclt(e)", "EX_h(e)",\
# #              "EX_h2(e)", "EX_succ(e)"]
# # print()
# # with open("ecoli-anaerobic-output.txt", "w") as f:
# #     for x in anaerobic:
# #         print(f"{x} : " + fba_solution.fluxes.get(x).__str__(), file=f)
# #     print("", file=f)
# #     print("Time elapsed to solve the FBA problem: " + str(end - start), file=f)
# #     print("Objective value: " + str(fba_solution.objective_value), file=f)
# #     print("", file=f)

# #     for x in E_coli_model.exchanges:
# #         y = x.id
# #         if fba_solution.fluxes.get(key=y):
# #             if fba_solution.fluxes.get(key=y) > 0:
# #                 print('exchange flux metabolites secreted', x.name, \
# #                       fba_solution.fluxes.get(key=y), file=f)

# # #The default is for FBA to optimize biomass. That's why this function will produce the same answer as above
# # E_coli_model.objective = 'Ec_biomass_iJO1366_core_53p95M'
# # fba_solution_biomass = E_coli_model.optimize()
# # print("FBA Solution Biomass", fba_solution_biomass)

# #If we want to run FBA to optimize a different reaction, let's say we want to see what the minimal glucose uptake is, we can have it optimize for glucose
# #note that FBA will maximize the specified reaction.
# #Since EX_glc(e) is going to take on a negative flux since it's being taken in by the model, this is asking for the minimum amount of glucose needed for the microbe to live (but not grow).
# E_coli_model.objective = 'EX_glc(e)'
# fba_solution_glc_biomass_growth = E_coli_model.optimize().fluxes.Ec_biomass_iJO1366_core_53p95M
# print("FBA Solution Glucose Biomass Growth", fba_solution_glc_biomass_growth)
# fba_solution_glc_minimal_uptake = E_coli_model.optimize().fluxes.get('EX_glc(e)')
# print("FBA Solution Glucose Minimal Uptake", fba_solution_glc_minimal_uptake)

# # #For more on FBA, read the following:
# # #https://cobrapy.readthedocs.io/en/latest/simulating.html

# # #Let's say I want to find out what the biomass growth of an E coli would be on a sucrose medium instead of a glucose medium:
# # #And I want to have a separate model with this condition set
# # #How to Copy a Model
# # import copy
# # E_coli_model_sucrose = copy.deepcopy(E_coli_model)

# # #Setting new bounds on sugars and less oxygen, just to see what happens
# # E_coli_model_sucrose.reactions.get_by_id('EX_sucr(e)').bounds = (-10, 1000)
# # E_coli_model_sucrose.reactions.get_by_id('EX_fru(e)').bounds = (0, 1000)
# # E_coli_model_sucrose.reactions.get_by_id('EX_glc(e)').bounds = (0, 1000)
# # E_coli_model_sucrose.reactions.get_by_id('EX_o2(e)').bounds = (-2, 1000)

# # #new FBA problem
# # E_coli_model_sucrose.objective = 'Ec_biomass_iJO1366_core_53p95M'
# # print("Sucrose model biomass FBA Solution", E_coli_model_sucrose.optimize())

# # #Going back to our original model, let's learn a little more.

# # #let's say I want to know what metabolites are secreted out of the model if the model is optimizing growth (as we can assume a microbe will do)
# # #I'll use a for loop to run through all the exchange reactions, and print any with a positive flux value
# # #Remember that a positive flux means secretion
# # for x in E_coli_model.exchanges:
# #     y = x.id
# #     if fba_solution_biomass.fluxes.get(key=y):
# #         if fba_solution_biomass.fluxes.get(key=y) > 0:
# #             print('exchange flux metabolites secreted', x.name, fba_solution_biomass.fluxes.get(key=y))

# # #How to edit a reaction
# # #If I want to add metabolites to a reaction, I can go in and edit the reaction
# # #In this example, I will be adding ATP to the D-amino acid dehydrogenase reaction
# # DAAD_rxn = E_coli_model.reactions.DAAD
# # print("Original DAAD Reaction", DAAD_rxn)
# # print(E_coli_model.reactions.DAAD.check_mass_balance())
# # DAAD_rxn.add_metabolites({E_coli_model.metabolites.get_by_id("atp_c"): -1,
# #     E_coli_model.metabolites.get_by_id("adp_c"): 1,
# #     E_coli_model.metabolites.get_by_id("pi_c"): 1,
# #     E_coli_model.metabolites.get_by_id("h_e"): 1
# #     })
# # print("Modified DAAD Reaction", DAAD_rxn)
# # print(E_coli_model.reactions.DAAD.check_mass_balance())

# # #Now, I can take away the metabolites.
# # DAAD_rxn.subtract_metabolites({E_coli_model.metabolites.get_by_id("atp_c"): -1,
# #     E_coli_model.metabolites.get_by_id("adp_c"): 1,
# #     E_coli_model.metabolites.get_by_id("pi_c"): 1,
# #     E_coli_model.metabolites.get_by_id("h_e"): 1
# #     })
# # print("Back to Original DAAD Reaction", DAAD_rxn)
# # print(E_coli_model.reactions.DAAD.check_mass_balance())

# # #Gene Knockouts
# # # Knock out all reactions in the model and report the reaction deletions that lead to a zero max biomass flux
# # # (these are called essential reactions).
# # # A for loop to run through all of the reactions deleting one at a time and extract the solution. Print if it is 0.
# # E_coli_model.objective = 'Ec_biomass_iJO1366_core_53p95M'
# # fba_solution_biomass = E_coli_model.optimize()
# # print('complete model: ', fba_solution_biomass)
# # for m in E_coli_model.reactions:
# #      name_rxn = m.id
# #      with E_coli_model:
# #          E_coli_model.reactions.get_by_id(name_rxn).knock_out()
# #          optimize_solution_knockout = E_coli_model.optimize()
# #          if optimize_solution_knockout.objective_value < .001 and optimize_solution_knockout.objective_value > -.001:
# #             print('Model Solutions', name_rxn, E_coli_model.optimize())

# # #I hope this has been a useful demo! If you are looking for more information, just look at the following website:
# # #https://cobrapy.readthedocs.io/en/latest/
