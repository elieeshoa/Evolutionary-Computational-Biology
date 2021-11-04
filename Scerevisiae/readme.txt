--------------- Model ------------------
Scerevisiae_iAZ900_noCycles_03_25_2015.xml
SBML file of the model (PMID: 21190580)

The iAZ900 model has two biomass reactions: biomass_core and biomass_wildType. For simulations, use the former.

iAZ900.txt
The list of all reactions in the iAZ900 model

--------------- Growth media -------------
iAZ900_minimal.py:
Minimal media bounds 

Different carbon sources that you can try include:
Glucose EX_glc_e: [LB, UB] = [-10, 1000]
Sucrose: EX_surc_e: [LB< UB] = [-10, 1000]
Fructose: EX_fru_e: [LB, UB] = [-10, 1000]

For aerobic condition set the UB on oxygen uptake to -2: EX_o2_e: [LB, UB] = [-2, 1000].
For anaerobic condition set the LB for oxygen uptake to zero.

