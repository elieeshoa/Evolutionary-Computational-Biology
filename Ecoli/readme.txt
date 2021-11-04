-------------------------- SBML model file --------------------------- 
Ecoli_iJO1366_updated.xml:
E. coli iJ1366 model (PMID:  21988831)

---------------- Growth media files ---------------------------------
Minimal media flux bounds:
iJO1366_minimal_aerobic.py		
Minimal media, aerobic condition, any carbon source

iJO1366_minimal_glucose_aerobic.py
Minimal media, aerobic condition, glucose as the carbon source

iJO1366_minimal_anaerobic.py		
Minimal media, anaerobic condition, any carbon source

iJO1366_minimal_glucose_anaerobic.py
Minimal media, anaerobic condition, glucose as the carbon source

Set the flux of hte carbon source to -10
Glucose: EX_glc(e) = [LB, UB] = :[-10,1000]

For aerobic condition, set the LB on the flux of exchange reaction for oxygen to two times that for the limiting carbon source
EX_o2(e) = [LB, UB] = [-20, 1000]   




