from __future__ import division
import sys, time
sys.path.append('../../')
from tools.globalVariables import *
from fbaTools import fbaTools
from fba import fba
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from pyomo import environ # It was previously: "from coopr import pyomo" for versions of pyomo older than 4.X

# The following lines change the temporary directory for pyomo
from pyutilib.services import TempfileManager
TempfileManager.tempdir = pyomo_tmp_dir

def fva(model, selected_rxns = [], optimization_solver = default_optim_solver, save_to_model = False, results_filename = '', simulation_conditions = '', warmstart = False, warnings = True, stdout_msgs = True):
    """
    Performs flux variability analysis

    INPUTS:
    -------
    model: 
    An instance of class model containing the information
    about the metabolic model

    selected_rxns:
    A list (or tuple) of selected reactions for which FVA should be performed. If no
    input is provided FVA is performed for all reactios in the model

    optimization_solver: 
    Name of the LP solver to be used to solve the LP. Current 
    allowable choices are cplex and gurobi

    save_to_model: 
    If True, it stores the identified bounds on reaciton fluxes in fva_flux_bounds. Otherwise
    they are stored in a dictionary whose keys are ids and values are a list of two elements 
    in the form [fva_LB,fva_UB], whith fva_LB and fva_UB being the FVA LB and UB on fluxes 

    results_filename: 
    A string containing the name of the file to save the results in. If an empty string 
    is provided the results are not saved to a file

    simulation_conditions: 
    A string describing simulation conditions


    OUTPUTS:
    --------
    fva_flux_bounds: 
    A dictionary with keys being reactions ids and values beiing a list ot two elements 
    containing the fva flux bounds in the form [LB, UB]

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 08-22-2017 
    """   
    # save_to_model
    if not isinstance(save_to_model,bool):
        raise TypeError('save_to_model must be either True or False')

    # selected_rxns
    if not isinstance(selected_rxns,list) and not isinstance(selected_rxns, tuple):
        raise userError('selected_rxns must be a list or tuple of reaction objects')
 
    # optimization_solver
    if not isinstance(optimization_solver,str):
        raise TypeError('optimization_solver must be a string')
    elif optimization_solver.lower() not in ['gurobi','cplex','gurobi_ampl','cplexamp']:
        raise ValueError('Invalid value for optimization_solver. Allowed choices are gurobi and cplex')

    # simulation_conditions 
    if not isinstance(simulation_conditions,str):
        raise TypeError('simulation_conditions must be a string')

    # warmstart 
    if not isinstance(warmstart,bool):
        raise TypeError('use_warmsart must be either True or False')

    # warnings and stdout_msgs 
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True or False')
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True or False')

    # If warmstart is True use gurobi_ampl
    if warmstart and optimization_solver in ['gurobi', 'gurobi_ampl']:
        optimization_solver = 'gurobi_ampl'
    elif warmstart and optimization_solver not in ['gurobi', 'gurobi_ampl']:
        # If the solver is not gurobi or gurobi_ampl warmstart should be turned off
        # because other solvers such as gurobi will return an error (see fbaTools.py)
        warmstart = False
        print '**WARNING (fva.py)! warmstart was turned off becasue it can be used only with gurobi or gurobi_ampl as the solver. The specified solver ({}) does not support warmstart'.format(optimization_solver)

    # A dictionary holding the FVA flux bounds
    fva_flux_bounds = dict([(r.id,[None,None]) for r in model.reactions]) 

    #--- Minimize rxn fluxes ---
    for rxn in model.reactions:
        rxn.objective_coefficient = 0

    # Reactions to consider
    if len(selected_rxns) == 0:
        rxns_to_consider = model.reactions
    else:
        rxns_to_consider = selected_rxns

    counter = 0
    for rxn in rxns_to_consider: 
        counter += 1
       
        rxn.objective_coefficient = 1

        if counter == 1:
            fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = False, save_to_model = False, simulation_conditions = simulation_conditions, warmstart = warmstart, warmings = warnings, stdout_msgs = False, show_solver_output = False)

        # From counter 2 on, turn off build_new_optModel and preprocessing and turn on warmstart 
        elif counter == 2:
            fba_model.build_new_optModel = False

        # Redefine the objective function if counter > 1
        if counter > 1:
            fba_model.optModel.del_component('objectiveFunc')
            fba_model.optModel.objectiveFunc = environ.Objective(rule = fba_model.objectiveFunc_rule, sense = environ.minimize)

            # Supply the current solution as the warm start
            for j in fba_model.optModel.J:
                fba_model.optModel.v[j] = fba_model.solution['opt_rxnFluxes'][j]

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            LB = fba_model.solution['objective_value']
        else:
            raise userError('fba problem to find LB for rxn {} in fva did not end with an optimal solution: exit_flag = {}'.format(rxn.id, fba_model.solution['exit_flag']))

        # Store the results
        if save_to_model:
            rxn.fva_flux_bounds[0] = LB
        else:
            fva_flux_bounds[rxn.id][0] = LB

        rxn.objective_coefficient = 0

    #--- Maximize rxn flux ---
    for rxn in model.reactions:
        rxn.objective_coefficient = 0

    counter = 0
    for rxn in rxns_to_consider: 
        counter += 1
       
        rxn.objective_coefficient = 1

        if counter == 1:
            fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = True, save_to_model = False, simulation_conditions = simulation_conditions, warmstart = warmstart, warnings = warnings, stdout_msgs = False, show_solver_output = False)

        # From counter 2 on, turn off build_new_optModel and preprocessing and turn on warmstart 
        elif counter == 2:
            fba_model.build_new_optModel = False

        # Redefine the objective function if counter > 1
        if counter > 1:
            fba_model.optModel.del_component('objectiveFunc')
            fba_model.optModel.objectiveFunc = environ.Objective(rule = fba_model.objectiveFunc_rule, sense = environ.maximize)

            # Supply the current solution as the warm start
            for j in fba_model.optModel.J:
                fba_model.optModel.v[j] = fba_model.solution['opt_rxnFluxes'][j]

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            UB = fba_model.solution['objective_value']
        else:
            raise userError('fba problem to find UB for rxn {} in fva ended with a non-optimal solution: exit_flag = {}'.format(rxn.id, fba_model.solution['exit_flag']))

        # Store the results
        if save_to_model:
            rxn.fva_flux_bounds[1] = UB
        else:
            fva_flux_bounds[rxn.id][1] = UB

        rxn.objective_coefficient = 0

        # Save results into a file 
        if results_filename != '':
            with open(results_filename,'w') as f:
                f.write('fva_flux_bounds = {\n')
                for rxn in fva_flux_bounds.keys():
                    f.write("'{}':{},\n".format(rxn, fva_flux_bounds[rxn]))
                f.write('}')

    # Return fva_flux_bounds if save_to_model is not True
    return fva_flux_bounds        

