"""
-----------------------------------------------
--------------- NashEq Finder -----------------
-----------------------------------------------
This script is a python implementation of the NashEq Finder algorithm presented in

Zomorrodi, AR, Segre, D, "Microbial games at genomic resolution: understanding the 
evolution of intercellular metabolic interactions in microbial communities", Nat Comm (2017)

This code can identify all pure strategy Nash equilibria of a game with any number of players
and strategies inone shot.

NOTE:
1. This code requires installing pyomor, which is a python-based optimization modeling software package
   Check out the followingn link for details:
   http://www.pyomo.org

2. This code also requires an optimizaton solver such as gurobo or IBM cplex. Consult the respected 
   website for further details.   

Ali R. Zomorrodi, Segre Lab @ Boston University
Last updated: July 06th, 2017

Please contact Ali Zomorrodi at ali.r.zomorrodi@gmail.com for questions and updates

"""

from __future__ import division
import datetime
from operator import concat
import re, sys, math, copy, time, random
from datetime import timedelta
from tracemalloc import start
from uuid import uuid4
from numpy import nonzero  # To convert elapsed time to hh:mm:ss format
from sympy.sets.sets import FiniteSet
sys.path.append('/content/drive/MyDrive/Manuscript/Ali_codes/')
import optlang # 1.5.2
# from optlang.cplex_interface import Model, Variable, Constraint, Objective
from optlang.gurobi_interface import Model, Variable, Constraint, Objective
import sympy # 1.8
import matplotlib.pyplot as mpl

# import pretty print
from pprint import pprint


class NashEqFinder(object):
    """
    General class for NashEq Finder. Sample usage is provided at the end 
    """   

    def __init__(self, game, NashEq_type = 'pure', optimization_solver = 'gurobi', warnings = True, stdout_msgs=False, output_file = '', stdout_timing=True):
        """
        INPUTS 
        ------
        game: 
        An instance of the class game (see game.py for details) 

        NashEq_type:
        Type of the Nash equilibrium to find (currently only pure strategy Nash equilibrium)

        optimization_solver: 
        Name of the LP solver to be used to solve the LP. Current 
        allowable choices are cplex and gurobi

        warnings: 
        Can be True or False indicating whether warnings should be written 
        in the standard output

        stdout_msgs: 
        By default (True) writes a summary including the solve 
        status, optimality status (if not optimal), objective 
        function value and the elapsed time on the screen.
        if set to a value of False no resuults are written on 
        the screen, in which case The user can instead specifiy 
        an output file using the option output_file, or store 
        them in a variable (see the 'run' method for details)

        output_file: 
        Optional input. It is a string containg the path to a 
        file and its name (e.g., 'results1/fbaResults.txt'), where
        the results should be written to. 
        """
       
        # Metabolic model
        self.game = game

        # Type of the Nash equilibrium to find
        if NashEq_type.lower() not in ['pure','mixed']:
            raise ValueError("Invalid NashEq_type (allowed choices are 'pure' or 'mixed')")
        else:
            self.NashEq_type = NashEq_type

        # Solver name
        if optimization_solver == None:
            self.optimization_solver = 'gurobi'
        else:
            if optimization_solver.lower() in ['cplex','gurobi']:
                self.optimization_solver = optimization_solver
            else:
                raise ValueError('Invalid solver name (eligible choices are cplex and gurobi)\n')          
               
        # Output to the screen 
        if not isinstance(warnings,bool):
            raise TypeError("Error! warnings should be True or False")
        else:
             self.warnings = warnings

        if not isinstance(stdout_msgs,bool):
            raise TypeError("Error! stdout_msgs should be True or False")
        else:
            self.stdout_msgs = stdout_msgs

        if not isinstance(stdout_timing,bool):
            raise TypeError("Error! stdout_timing should be True or False")
        else:
            self.stdout_timing = stdout_timing

        # Output file
        if not isinstance(output_file, str):
            raise TypeError('output_file must be a string')
        else:
            self.output_file = output_file

        # Lower bound on payoff values according to the payoff matrix
        payoffMin = min([k for sublist in self.game.payoff_matrix.values() for k in sublist.values()]) 

        # Sometimes we run into problems when both LB and max payoff of a plyaer given
        # the fixed strategies of other players are zero (i.e., we arrive at a trivial
        # solution of e.g., 2 >= 0 as both terms in the RHS of constraint NashCond are
        # Cancelled out. This happens for problem 9 of Homework 1 of Game Theory I
        # for example). Therefore, it is better to always avoid a LB of zero. 
        if payoffMin - 1 > 0: 
            self.payoffLB = payoffMin - 1 
        else:
            self.payoffLB = payoffMin - 2 


    def convert_to_payoffMatrix_key(self,i):
        """
        This function converts the elements of the set I in the pyomo model
        (or elements of gameStatesForI) to the format of keys of the payoff matrix
        of the game, i.e., ('p1','s1','p2','s2') is converted to the format of
        (('p1','s1'), ('p2','s2')). See optModel.I for how the set I is defined
        """
        gameState = []
        done = 0
        k1 = list(i)
        while done == 0:
            gameState.append(tuple(k1[0:2]))
            del k1[0:2]
            if len(k1) == 0:
                done = 1
        return tuple(gameState)

        
    # Elie
    def createOptlangModel(self):
        """
        This creates a optlang optimization model 
        """   

        def tuple_to_var_name(t):
            return str(t).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',','_')
        
        def add_optlang_NashCond_rule(optlangOptModel,p,i):

            # Convert the game state to the format of keys of the payoff matrix
            i = self.convert_to_payoffMatrix_key(i)
            
            responseP = [
                k for k in self.game.payoff_matrix.keys() if False not in \
                [dict(k)[pp] == dict(i)[pp] for  pp in dict(i).keys() if pp != p]
                ]

            # Find the payoff of the best response of player P 
            bestResP = max([self.game.payoff_matrix[k][p] for k in responseP])

            model.add(Constraint(
                self.game.payoff_matrix[i][p] - \
                    bestResP*optlangOptModel.variables[tuple_to_var_name(i)] - \
                    self.payoffLB * (1 - optlangOptModel.variables[tuple_to_var_name(i)]),
                lb=0)) 

        # set optlang solver as cplex
        # optlang.cplex_interface

        model = Model(name='Original Optlang Model')
        # print("hi from createOptlangModel")
        # exit()
        model.players = self.game.players_names
        model.indices = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()]

        variables_names = []

        # Add the variables to the model
        for index in model.indices:
            var_str = tuple_to_var_name(index)
            var = Variable(var_str, type='binary', problem=model)
            variables_names.append(var_str)
            model.add(var)

        # Add the objective function
        model.objective = \
            Objective(
                expression=sympy.Add(*sympy.symbols(variables_names)), \
                direction='max')
   
        # Add the constraints
        for player in model.players:
            for index in model.indices:
                add_optlang_NashCond_rule(model, player, index)
        
        self.optModel = model

    def optlangFindPure(self):
        """ 
        This method runs the optimization problem finding the pure strategy Nash
        equilbirium. 

        OUTPUTS:
        -------
        Nash_equilibria: 
        Is a list containing the labels of the cells of the payoff matrix
        that were found to be a pure strategy Nash equilibrium. For example, in a two-player 
        game if the set of strategies for players 1 and 2 are {s11,s12} and {s21,s22},
        respectively, the optimal values of binary varaibles for each cell can be as follows 
        {('s11','s21'):0,('s11','s21'):1,('s12','s21'):0,('s21','s22'):0}
        and additionally we may have an alternative solution as:
        {('s11','s21'):0,('s11','s21'):0,('s12','s21'):0,('s21','s22'):1}
        Nash_equilibria would be then be a list [('s11','s21'),('s21','s22')] 

        exit_flag: 
        Shows the condition the termination condition of the code (this is different from 
        optimExitflag for solving the optimization problem). exit_flag can take 
        either of the following values:
        - 'objIsZero': The objective function is zero
        - 'solverError': There was an error in both optimization solvers (cplex & guorobi)
        - 'objNotZeroNotOne': An erroneous case where the objective function is neither
                              zero nor one
        - A string showing a non-optimal solution for the optimization problem     
        """
        # Helper function to convert a tuple to a string representing a variable name
        def tuple_to_var_name(t):
            return str(t).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',','_')
        
        # Processing and wall time
        start_run_pt = time.process_time()
        start_run_wt = time.time()

        #---- Creating and instantiating the optModel ----
        start_optlang_pt = time.process_time()
        start_optlang_wt = time.time()

        

        # Create the self.optModel model    
        self.createOptlangModel()
        # print("hi from optlangFindPure")
        # exit()

        #---- Solve the model ----
        # Create a solver and set the options
        elapsed_optlang_pt = \
            str(timedelta(seconds = time.process_time() - start_optlang_pt))
        elapsed_optlang_wt = \
            str(timedelta(seconds = time.time() - start_optlang_wt))

        #- Solve the optModel
        start_solver_pt = time.process_time()
        start_solver_wt = time.time()

        # optSoln = solverType.solve(self.optModel,tee=False)
        optSoln = self.optModel.optimize()
        solverFlag = 'normal'
    
        # Elie
        elapsed_solver_pt = \
            str(timedelta(seconds = time.process_time() - start_solver_pt))
        elapsed_solver_wt = \
            str(timedelta(seconds = time.time() - start_solver_wt))
    
            
        # Set of the Nash equilibria
        self.Nash_equilibria = []
        objValue = self.optModel.objective.value
        # TODO: Check if this is the correct way to check if the objective is zero
        #       this assumes .variables are all binary
        if objValue >= 1:
            self.exit_flag = 'objGreaterThanZero'
            for i in self.optModel.indices: 
                if self.optModel.variables[tuple_to_var_name(i)].primal == 1:
                    self.Nash_equilibria.append(list(self.convert_to_payoffMatrix_key(i)))
        elif objValue == 0:
            done = 1
            self.exit_flag = 'objIsZero'
    
        # Time required to run 
        elapsed_run_pt = str(timedelta(seconds = time.process_time() - start_run_pt))
        elapsed_run_wt = str(timedelta(seconds = time.time() - start_run_wt))
    
        if self.stdout_msgs:
            print('NashEqFinder took (hh:mm:ss) (processing/wall) time: pyomo\
                   = {}/{}  ,  solver = {}/{}  ,  run = {}/{} \n'.format(elapsed_optlang_pt,\
                   elapsed_optlang_wt, elapsed_solver_pt,elapsed_solver_wt, \
                   elapsed_run_pt,elapsed_run_wt))
    
    def optlangRun(self):
        """
        Runs the Nash equilibrium finder
        """
        if self.NashEq_type.lower() == 'pure':
            self.optlangFindPure()
        elif self.NashEq_type.lower() == 'mixed':
            pass # To be completed

        return [self.Nash_equilibria, self.exit_flag, self.game.payoff_matrix]


    def show_matrix_2c(self, original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells):  
        fig, ax = mpl.subplots()
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        def payoffs_to_table(payoff_matrix): 
            col_text =  [''] + strategies
            table_text = []
            for i in strategies:
                row = [i]
                for j in strategies:
                    og1 = original_payoff_matrix[(('row', i),('column', j))]['row']
                    og2 = original_payoff_matrix[(('row', i),('column', j))]['column']
                    new1 = payoff_matrix[(('row', i),('column', j))]['row']
                    new2 = payoff_matrix[(('row', i),('column', j))]['column']
                    def check_0(s):
                        if s != '0.0' and s != '-0.0' and s != '0':
                            return f" + {s}"
                        else:
                            return ""

                    def remove_zeroes(s):
                        if s[-2:] == '.0':
                            return s[:-2]
                        else:
                            return s

                    row.append(
                        f"({remove_zeroes(str(round(float(og1), 4)))}" + check_0(str(round(new1 - float(og1), 4))) + f", {remove_zeroes(str(round(float(og2), 4)))}" + check_0(str(round(new2 - float(og2), 4))) + ")"
                    )
                table_text.append(row)
            return col_text, table_text

        if self.stdout_timing:
            print(f"Time before defining the table: {datetime.datetime.now()}")
        table = payoffs_to_table(payoff_matrix)
        the_table = ax.table(#colWidths=[0.3] * (SIZE + 1),
                            cellText=table[1], colLabels=table[0],
                            loc='center', bbox=[-0.05, 0, 1.1, 1.0],
                            rowLoc='center', colLoc='center', cellLoc='center')
        the_table.scale(1, 7)

        # changed cell is of format (
        #   (
        #       (('row', 'S15'), ('column', 'S10')), 'row', 'plus'
        #   ), 7.009999999999547
        # )
        def changed_cell_to_position(changed_cell):
            return(
                strategies.index(changed_cell[0][0][0][1])+1,
                strategies.index(changed_cell[0][0][1][1])+1
            )
        
        if self.stdout_timing:
            print(f"Time before setting the cells colors: {datetime.datetime.now()}")

        changed_cells_positions = [changed_cell_to_position(cell) for cell in changed_cells]
        for (row, col) in changed_cells_positions:
            the_table[row, col].set_facecolor('#90ac44')

        def letter_to_position(letter):
            return strategies.index(letter) + 1
        
        if self.stdout_timing:
            print(f"Time before setting the cells colors again: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].set_linewidth(3)

        if self.stdout_timing:
            print(f"Time before closed grid edges: {datetime.datetime.now()}")
        for (row, col), cell in the_table.get_celld().items():
            cell.visible_edges = ''
            if col != 0 and row != 0:
                cell.visible_edges = 'closed'

        nasheq_positions = []

        if self.stdout_timing:
            print(f"Time before adding nasheq_positions: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            nasheq_positions.append(
                (letter_to_position(nash_equilibria[eq][0][1]), \
                letter_to_position(nash_equilibria[eq][1][1]))
            )

        size = len(strategies)
        the_table.auto_set_font_size(False)
        if size == 2:
            for row in range(size+1):
                for col in range(size+1):
                    if (row, col) not in nasheq_positions:
                        the_table[row, col].set_fontsize(14)
                    else:
                        the_table[row, col].set_fontsize(12)

        elif size == 5:
            for row in range(size+1):
                for col in range(size+1):
                    if (row, col) not in nasheq_positions:
                        the_table[row, col].set_fontsize(8)
                    else:
                        the_table[row, col].set_fontsize(8)
                        
        elif size == 10:
            for row in range(size+1):
                for col in range(size+1):
                    if (row, col) not in nasheq_positions:
                        the_table[row, col].set_fontsize(4)
                    else:
                        the_table[row, col].set_fontsize(4)

        if self.stdout_timing:
            print(f"Time before saving the PNG: {datetime.datetime.now()}")
        if self.stdout_msgs:
            print(f"saving the PNG file for size {len(strategies)}")
        # set title based on iteration
        # if iteration is in method.split(' ') find the word after it
        # if iteration is not in method.split(' ') then just use the method
        if 'iteration' in method.split(' '):
            title = f"Alternative {method.split(' ')[method.split(' ').index('iteration')+1][:-1]}"
        elif 'First' in method:
            title = 'First Solution (Before Finding Alternatives)'
        else:
            title = 'First Solution (Before Finding Alternatives)'
        
        # set title with big font
        ax.set_title(f"{title}", fontsize=20)
        mpl.show()

    
    # a function that checks if the desired cells are nash equilibria
    def check_if_nash_equilibria(self, desired_cells, payoff_matrix, strategies):
        """
        this function returns True if the desired cell is a nash equilibrium
        and in the payoff matrix, and False otherwise

        desired_cells: a list of elements of format (('row', 'S15'), ('column', 'S10'))
        """
        for cell in desired_cells:
            # cell is of format [(('row',f'S{i}'), ('column',f'S{j}'))]
            # check that the first entry is greater than all first entries that are in the same column
            # check that the second entry is greater than all second entries that are in the same row
            # if both are true, then the cell is a Nash Equilibrium
            # for all payoff_matrix cells that start with ('row',f'S{i}'), check that they are greater than cell[0][1]
            row_values = [
                self.game.payoff_matrix[('row', strategy), cell[1]]['row'] <= 
                self.game.payoff_matrix[cell]['row'] for strategy in strategies 
                if ('row', strategy) != cell[0]
                ]
            column_values = [
                self.game.payoff_matrix[cell[0], ('column', strategy)]['column'] <= 
                self.game.payoff_matrix[cell]['column'] for strategy in strategies 
                if ('column', strategy) != cell[1]
                ]
            if all(row_values) and all(column_values):
                return True
            else:
                return False


    # Elie on 9/5/2022
    def validate_2c(self, nasheq_cells, method, iteration, time_spent, size):
        # Validation: needs removal of hard coded methods
        # Helper function
        def string_to_index(string):
            lst = string.split('_')
            player = lst[-2]
            sign = lst[-1]
            return (((lst[0],lst[1]), (lst[2],lst[3]))), player, sign

        original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        
        if self.stdout_timing:
            print(f"Time before creating the matrix of the perterbations: {datetime.datetime.now()}")
        # Creating new payoff matrix (i.e. the result of perterbations)
        # TODO: check if this is the right way to do it
        for var_name, var_primal in self.current_primals.items():
            matrix_key, player, sign = string_to_index(var_name)
            if sign == 'plus':
                self.game.payoff_matrix[matrix_key][player] += var_primal
            if sign == 'minus':
                self.game.payoff_matrix[matrix_key][player] -= var_primal

        # Get the cells that changed
        if self.stdout_timing:
            print(f"Time before getting the changed cells: {datetime.datetime.now()}")
        changed_cells = []
        epsilon = 0.01
        for var_name, var_primal in self.current_primals.items():
            if var_primal >= epsilon and var_name not in (self.current_binary_variables + self.current_y_binary_variables + self.b_i_variables + self.r_c_max_variables):
                changed_cells.append((string_to_index(var_name), var_primal))

        # Write changed cells into a file
        if self.stdout_timing:
            print(f"Time before writing the changed cells: {datetime.datetime.now()}")

        # print changed cells
        if iteration != "NONE":
            print("Changed cells for iteration " + str(int(iteration)+1) + ":")
            pprint(changed_cells)

        
        # Check that the desired cells are Nash Equilibria
        if self.stdout_timing:
            print(f"Time before checking that the the desired cells are Nash Equilibria: {datetime.datetime.now()}")
        
        success = self.check_if_nash_equilibria(nasheq_cells, self.game.payoff_matrix, self.game.players_strategies['row'])
        # if self.stdout_msgs:
        if success:
            print(f"The desired cell(s) {nasheq_cells} are Nash Equilibria\n")
        else:
            print(f"The desired cell(s) {nasheq_cells} are NOT Nash Equilibria\n")
                
        # add to the file
        if self.stdout_timing:
            print(f"Time before writing to Results Summer 2023: {datetime.datetime.now()}")

        self.show_matrix_2c(original_payoff_matrix, self.game.payoff_matrix, nasheq_cells, self.game.players_strategies['row'], method, changed_cells)
  

#     # Elie
#     def newEquilibria(self, nasheq_cells, strategies, original_nasheq_cells):
#         """
#         :param nasheq_cells: a list of elements (('row','C'),('column','C'))
#         Consider a payoff value a. We would like to perturb it such that it can
#         either increase or decrease. To do this, we define two non-negative 
#         variables aa^+ and aa^-. Then we change the payoff as follows:
#                                     β=a+aa^+-aa^-
#         As the objective function then, you minimize sum of all aa^+'s and 
#         aa^-'s for all payoff. 
#         return: a model with the solutions' final model, which contains the 
#                 optimal value of the variables and more
#         """   
        
#         def tuple_to_var_name(t):
#             return str(t).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',','_')
        

#         # set tolerance for optlang
#         print("optlang.interface is: ", optlang.interface)
#         tol = 1e-9
#         self.optModel.configuration.tolerances.feasibility = tol
#         self.optModel.configuration.tolerances.optimality = tol
#         self.optModel.configuration.tolerances.integrality = tol
      
#         # Adding new variables 
#         model = Model(name=f'Original Model')
#         model.players = self.game.players_names
#         if self.stdout_timing:
#             print(f"Time before creating indices for original model: {datetime.datetime.now()}")
#         indices = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()]
#         new_indices = []
#         for index in indices:
#             for player in model.players:
#                 new_indices.append(index + (player, 'plus'))
#                 new_indices.append(index + (player, 'minus'))

#         # Format of indices: e.g. ('row','C','column','C','row','plus')
#         model.indices = new_indices
#         variables_names = []

#         # Add the variables
#         if self.stdout_timing:
#             print(f"Time before adding variables to original model: {datetime.datetime.now()}")
#         for index in model.indices:
#             var_name = tuple_to_var_name(index)
#             var = Variable(var_name, lb=0, type='continuous', problem=model)
#             variables_names.append(var_name)
#             model.add(var)

#         # Add the objective function SLOW
#         model.objective = Objective(expression=sympy.Add(*(sympy.symbols(variables_names))), direction='min')

#         # Add the constraints
#         constraints = []
#         # Each `cell` is of the format (('row','C'),('column','C'))
#         if self.stdout_timing:
#             print(f"Time before adding constraints to original model: {datetime.datetime.now()}")

#         self.b_i_variables = []
#         self.b_i_variables_optlang = []
#         self.r_c_max_variables = []
#         self.r_c_max_variables_optlang = []

#         # Add constraints that keep the original Nash equilibrium entries constant
#         for cell in original_nasheq_cells:
#             root_index = tuple_to_var_name(cell)
#             for player in model.players:
#                 c1 = Constraint(-model.variables[root_index+'_'+player+'_plus'], lb=0)
#                 c2 = Constraint(-model.variables[root_index+'_'+player+'_minus'], lb=0)
#                 constraints.append(c1)
#                 constraints.append(c2)
                
#         for cell in nasheq_cells:
#             root_index = tuple_to_var_name(cell)
#             # For first player we loop over first axis (rows=strategies)
#             for strategy in [x for x in strategies if x != cell[0][1]]:
#                 current_cell = ((cell[0][0], strategy), (cell[1]))
#                 current_index = tuple_to_var_name(current_cell)
#                 player = model.players[0]
#                 # This new snippet precludes cells other than `nasheq_cells` 
#                 # to be a Nash equilibrium
#                 epsilon_nash = 0.01
#                 c = Constraint(
#                         self.game.payoff_matrix[cell][player] \
#                         + model.variables[root_index+'_'+player+'_plus'] \
#                         - model.variables[root_index+'_'+player+'_minus'] \
#                         - 
#                         (self.game.payoff_matrix[current_cell][player] \
#                         + model.variables[current_index+'_'+player+'_plus'] \
#                         - model.variables[current_index+'_'+player+'_minus']
#                         )
#                         -
#                         epsilon_nash, lb=0)
#                 constraints.append(c)

#             # For second player we loop over second axis (rows=strategies)
#             for strategy in [x for x in strategies if x != cell[1][1]]:
#                 current_cell = ((cell[0]), (cell[1][0], strategy))
#                 current_index = tuple_to_var_name(current_cell)
#                 player = model.players[1]
#                 c = Constraint(
#                         self.game.payoff_matrix[cell][player] \
#                         + model.variables[root_index+'_'+player+'_plus'] \
#                         - model.variables[root_index+'_'+player+'_minus'] \
#                         - 
#                         (self.game.payoff_matrix[current_cell][player] \
#                         + model.variables[current_index+'_'+player+'_plus'] \
#                         - model.variables[current_index+'_'+player+'_minus']
#                         )
#                         -
#                         epsilon_nash, lb=0)
#                 constraints.append(c)

#         # b_ij^k=a_ij^k+〖α_ij^k〗^+-〖α_ij^k〗^-,
#         # rmax_rc≥b_ic^1,	∀(r,c)∈U,∀i∈S_1-{i|(i,c)∈U},	(4)
#         # rmax_rc≥b_rc^1,+ϵ	∀(r,c)∈U	(5)
#         # cmax_rc≥b_ri^2,	∀(r,c)∈U,∀i∈S_2-{i|(r,i)∈U},	(6)
#         # cmax_rc≥b_rc^2,+ϵ	∀(r,c)∈U	(7)
#         # rmax and cmax representing the maximum payoff of the Player 1 when 
#         # Player 2 takes strategy corresponding to c and cmax is the maximim 
#         # payoff of Player 2 when Player 1 takes the strategy corresponding 
#         # to r. Constraint (4) requires rmax_rc to be greater than all payoffs 
#         # in (i,c) except for the row r and Constraint (5) requires rmax_rc 
#         # to be strictly greater than the payoff of Player 1 in cell (r,c). 
#         # Constraints (6) and (8) impose similar constriants for cmax_rc 
#         # and the payoff of Player 2in cell (r,c).
#         # U = original_nasheq_cells
#         # D = desired_nasheq_cells
#         # S_1 = strategies of player 1
#         # S_2 = strategies of player 2
#         for cell in original_nasheq_cells:
#             root_index = tuple_to_var_name(cell)
#             # constraint (5): # rmax_rc≥b_rc^1,+ϵ
#             # get rmax_rc
#             # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] for r in strategies if r != cell[0][1]])
#             # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] for r in strategies])
#             # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] + model.variables[tuple_to_var_name((('row', r), ('column', cell[1][1])))+'_row_plus'] - model.variables[tuple_to_var_name((('row', r), ('column', cell[1][1])))+'_row_minus'] for r in strategies])

#             # https://math.stackexchange.com/questions/2446606/linear-programming-set-a-variable-the-max-between-two-another-variables
#             # \begin{align}
#             # U &\ge a_i   &\forall i \in N \\ 
#             # U &\le a_i + (1-b_i)*M  & \forall i \in N \\
#             # \sum_{i \in N} b_i &= 1
#             # \end{align}
#             # where the b_i\in{0,1} is a binary variable that indicates the maximum a_i
#             # (i.e. b_i=1 when a_i is the max value), and M it's a "big number".

#             # create the b_i variables (bi1 and bi2)
#             b_i1_variables = []
#             b_i1_variables_optlang = []
#             b_i2_variables = []
#             b_i2_variables_optlang = []
#             for r in strategies:
#                 # create the b_i1 variables
#                 str_index = root_index + "_binary1_preculde_original_nash_equilibria_" + r
#                 b_i1_variables.append(str_index)
#                 var = Variable(str_index, lb=0, ub=1, type='binary', problem=model)
#                 b_i1_variables_optlang.append(var)
#                 model.add(var)
#                 # create the b_i2 variables
#                 str_index = root_index + "_binary2_preculde_original_nash_equilibria_" + r
#                 b_i2_variables.append(str_index)
#                 var = Variable(str_index, lb=0, ub=1, type='binary', problem=model)
#                 b_i2_variables_optlang.append(var)
#                 model.add(var)



#             for binary_var in b_i1_variables:
#                 self.b_i_variables.append(binary_var)
#             for binary_var in b_i1_variables_optlang:
#                 self.b_i_variables_optlang.append(binary_var)
#             for binary_var in b_i2_variables:
#                 self.b_i_variables.append(binary_var)
#             for binary_var in b_i2_variables_optlang:
#                 self.b_i_variables_optlang.append(binary_var)

#             # define U, and call it root_index + "_rmac_rc"
#             str_rmax_rc = root_index + "_rmax_rc"
#             rmax_rc = Variable(str_rmax_rc, type='continuous', problem=model)
#             model.add(rmax_rc)
#             self.r_c_max_variables.append(str_rmax_rc)
#             self.r_c_max_variables_optlang.append(rmax_rc)
#             # create the constraints for rmax_rc
#             # U &\ge a_i   &\forall i \in N \\
#             # U &\le a_i + (1-b_i)*M  & \forall i \in N \\
#             for r in strategies:
#                 ic_cell = ((('row', r), cell[1]))
#                 ic_index = tuple_to_var_name(ic_cell)
#                 c = Constraint(
#                         rmax_rc \
#                         - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
#                         - model.variables[ic_index + '_row_plus'] \
#                         + model.variables[ic_index + '_row_minus'] \
#                         , lb=0)
#                 constraints.append(c)
                
#                 c = Constraint(
#                         rmax_rc \
#                         - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
#                         - model.variables[ic_index + '_row_plus'] \
#                         + model.variables[ic_index + '_row_minus'] \
#                         - (1 - model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r]) * 1000 \
#                         , ub=0)
#                 constraints.append(c)
#             # \sum_{i \in N} b_i1 <= 1
#             c = Constraint(
#                     sum([model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r] for r in strategies]) \
#                     , ub=1)
#             constraints.append(c)

#             # constraint (4): # rmax_rc≥b_ic^1,	∀(r,c)∈U,∀i∈S_1-{i|(i,c)∈U}
#             for r in strategies:
#                 if r != cell[0][1]:
#                     ic_cell = ((('row', r), cell[1]))
#                     ic_index = tuple_to_var_name(ic_cell)
#                     c = Constraint(
#                             rmax_rc \
#                             - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
#                             - model.variables[ic_index + '_row_plus'] \
#                             + model.variables[ic_index + '_row_minus'] \
#                             , lb=0)
#                     constraints.append(c)

#             # constraint (5): # rmax_rc≥b_rc^1,+ϵ
#             c = Constraint(
#                     rmax_rc \
#                     - self.game.payoff_matrix[cell]['row'] \
#                     - model.variables[root_index+'_row_plus'] \
#                     + model.variables[root_index+'_row_minus'] \
#                     - epsilon_nash, lb=0)
#             constraints.append(c)
          
#             # define U, call it root_index + "_cmax_rc"
#             str_cmax_rc = root_index + "_cmax_rc"
#             cmax_rc = Variable(str_cmax_rc, type='continuous', problem=model)
#             model.add(cmax_rc)
#             self.r_c_max_variables.append(str_cmax_rc)
#             self.r_c_max_variables_optlang.append(cmax_rc)
#             # create the constraints for cmax_rc
#             for strat in strategies:
#                 ri_cell = ((cell[0]), (('column', strat)))
#                 ri_index = tuple_to_var_name(ri_cell)
#                 c = Constraint(
#                         cmax_rc \
#                         - self.game.payoff_matrix[(('row', cell[0][1]), ('column', strat))]['column'] \
#                         - model.variables[ri_index + '_column_plus'] \
#                         + model.variables[ri_index + '_column_minus'] \
#                         , lb=0)
#                 constraints.append(c)

#                 c = Constraint(
#                         cmax_rc \
#                         - self.game.payoff_matrix[(('row', cell[0][1]), ('column', strat))]['column'] \
#                         - model.variables[ri_index + '_column_plus'] \
#                         + model.variables[ri_index + '_column_minus'] \
#                         - (1 - model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat]) * 1000 \
#                         , ub=0)
#                 constraints.append(c)
#             # \sum_{i \in N} b_i2 <= 1
#             c = Constraint(
#                     sum([model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat] for strat in strategies]) \
#                     , ub=1)
#             constraints.append(c)

#             # constraint (6): # cmax_rc≥b_ri^2,	∀(r,c)∈U,∀i∈S_2-{i|(r,i)∈U}
#             for c in strategies:
#                 if c != cell[1][1]:
#                     ri_cell = ((cell[0]), (('column', c)))
#                     ri_index = tuple_to_var_name(ri_cell)
#                     c = Constraint(
#                             cmax_rc \
#                             - self.game.payoff_matrix[(('row', cell[0][1]), ('column', c))]['column'] \
#                             - model.variables[ri_index + '_column_plus'] \
#                             + model.variables[ri_index + '_column_minus'] \
#                             , lb=0)
#                     constraints.append(c)

#             # constraint (7): # cmax_rc≥b_rc^2,+ϵ
#             c = Constraint(
#                     cmax_rc \
#                     - self.game.payoff_matrix[cell]['column'] \
#                     - model.variables[root_index+'_column_plus'] \
#                     + model.variables[root_index+'_column_minus'] \
#                     - epsilon_nash, lb=0)
#             constraints.append(c)

#             # constraint (10 from 12/29/2023): \sum_{i \in N} b_i1 + \sum_{i \in N} b_i2 >= 1
#             c = Constraint(
#                     sum([model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r] for r in strategies]) \
#                     + sum([model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat] for strat in strategies]) \
#                     , lb=1)
#             constraints.append(c)
                    


        
#         # print(constraints)
#         model.add(constraints)  
       
#         self.optModel = model  
#         start_time = datetime.datetime.now()
#         self.optModel.optimize()
#         end_time = datetime.datetime.now() 

#         with open(f"log Summer 2023.txt", "a") as f:
#             f.write(f"\nsize = {len(strategies)}, date: {datetime.datetime.now()}, Original Solution with desired nash equilibria: {nasheq_cells} || Time: {end_time - start_time}")

#         # Neccessary for validate_2c
#         self.current_variables = copy.deepcopy(variables_names)
#         self.current_primals = {}
#         if self.stdout_timing:
#             print(f"Time before self.current_primals[var_name] = var.primal: {datetime.datetime.now()}")
#         for var_name, var in self.optModel.variables.items():
#             # print(var_name, "=", var.primal)
#             self.current_variables.append(var_name)
#             self.current_primals[var_name] = var.primal

#         original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
#         # for validate_2c to work
#         self.current_binary_variables = []
#         self.current_binary_constraints = []
#         self.current_y_binary_variables = []
#         self.current_y_binary_variables_optlang = []
    

#         if self.stdout_timing:
#             print(f"Time before validate_2c on original model: {datetime.datetime.now()}")
#         self.validate_2c(nasheq_cells, method=f'First solution with new equilibria={nasheq_cells} precluded', iteration="NONE", time_spent=end_time-start_time, size=len(strategies))
#         # [NE, exit_flag, g_payoff_matrix] = self.optlangRun()
#         # print("matrix for original solution TOUT is: ", g_payoff_matrix)
#         # print("Nash_equilibria for original solution TOUT is: ", NE)
#         self.game.payoff_matrix = original_payoff_matrix

#         # For each non-zero α whose optimal value is α^opt, add the following
#         # constraints:
#         #                      α≤α^opt-ϵ  &  α≥ α^opt+ϵ

#         #                     α≤〖(1-y)(α〗^opt-ϵ)+y〖UB〗_α    
#         #                     〖α≥y(α〗^opt+ϵ)+(1-y)LB_α 

#         # Note that if y=0, the first constraint is reduced to α≤α^opt-ϵ and 
#         # the second constraint is reduced to α≥〖LB〗_α. On the other hand, 
#         # if y=1, the first constraint is reduced to α≤UB_α and the second 
#         # constraint is reducd to 〖α≥α〗^opt+ϵ. Instead of LB_α and UB_α, 
#         # you can use the big-M approach too, i.e., replaced them with -M and 
#         # M, respectively. 


#         # Adding new attributes to self
#         if self.stdout_timing:
#             print(f"Time before adding new attributes: {datetime.datetime.now()}")
#         self.current_variables = copy.deepcopy(variables_names)
#         self.current_binary_variables = []
#         self.current_binary_constraints = []
#         self.current_y_binary_variables = []
#         self.current_y_binary_variables_optlang = []
#         self.current_primals = {}
#         for var_name, var in self.optModel.variables.items():
#             # print(var_name, "=", var.primal)
#             # TODO: Check if this is needed
#             self.current_variables.append(var_name)
#             # TODO: Check if this is correct because validate_2c is using this.
#             #       Do we need primals for variables that are not in the original
#             #       payoff matrix?
#             self.current_primals[var_name] = var.primal

#         # adding a binary variable for each α
#         if self.stdout_timing:
#             print(f"Time before adding binary variables for each α: {datetime.datetime.now()}") # SLOW
#         for var_name, var_primal in self.current_primals.items():
#             str_index_y_a = var_name + "_binary" + f"_ya" # no iteration needed in this method 
#             y_a = Variable(str_index_y_a, lb=0, type='binary', problem=model)
#             model.add(y_a)
#             self.current_binary_variables.append(str_index_y_a)
#             # Add the binary y variables to the list of binary variables
#             self.current_y_binary_variables.append(str_index_y_a)
#             self.current_y_binary_variables_optlang.append(y_a)


       
#         objective_values = []

#         # For iteration = 1 to 5
#         for iteration in range(50):
#             # print the payoff matrix
#             # Elie on 1/27/2023
#             print(f"\n\n----iteration {iteration + 1}----\n\n")
#             # print the payoff matrix
#             if self.stdout_msgs:
#                 print("Payoff matrix at the beginning of iteration", iteration+1, ":")
#                 pprint(self.game.payoff_matrix)

#             set_binary_vars = set(self.current_binary_variables + self.current_y_binary_variables + self.b_i_variables + self.r_c_max_variables)

#             # print(f"----removing binary variables and constrains to start anew----\n")
#             if self.stdout_timing:
#                 print(f"Time before setting self.current_variables: {datetime.datetime.now()}")

#             self.current_variables = [e for e in self.current_variables if e not in set_binary_vars]
#             if self.stdout_msgs:
#                 print("self.current_variables for iteration", iteration+1, ":", self.current_variables)

#             if self.stdout_timing:
#                 print(f"Time before removing binary variables from self.current_primals: {datetime.datetime.now()}")

#             set_current_primals = set(self.current_primals.keys())
#             for key_to_remove in self.current_binary_variables:
#                 if key_to_remove in set_current_primals:
#                     del self.current_primals[key_to_remove]
#             # On june 20
#             epsilon = 0.01
#             ub = 1000
#             lb = -1000 

#             if self.stdout_timing:
#                 print(f"Time before adding for var_name, var_primal in self.current_primals.items(): {datetime.datetime.now()}")
#             for var_name, var_primal in self.current_primals.items():
#                 # Didn't work without it, which is weird
#                 # Check if the variable optimal is positive and if it is not
#                 # a binary variable
#                 if var_primal >= epsilon and var_name not in set_binary_vars:
#                     # I will introduce two non-linear constraints:
#                     # Define a binary variable y_α such that if y_α=1, then α≠0
#                     # and if y=0, then α=0 in the previous solution (note that 
#                     # here α could be either α^+ or α^-). This can be imposed 
#                     # using the following cosntrinat
#                     #               ϵy_α+(-M)(1-y_α)≤α≤(+M)y_α
#                     # Note that according to this constaint: 
#                     # y_α=1→α≠0
#                     # y_α=0→α=0

#                     # Keep in mind that:
#                     #     Here, ϵ is a small positive value (e.g., 10e-5)
#                     #     Here we define binary variables for each α^+ or α^- 
#                     #     only once, i..e, we do not define new binary variales 
#                     #     in each iteration. 
#                     # For the next solution, we want at least one perturbed 
#                     # payoff value to be different from the solution that we 
#                     # had before. We can easily impose this by using the 
#                     # following constraint:
#                     #              ∑_(α∈NZ)▒y_α ≤(card(NZ_α)-1)_+

#                     # if self.stdout_msgs:
#                     #     print("PIPTH in iteration", iteration+1, ":", var_name, "=", var_primal)
#                     #     print("var_primal", var_primal)

#                     # get the binary variable y_α
#                     y_a = model.variables[var_name + "_binary" + f"_ya"]
#                     if self.stdout_msgs:
#                         print("In iteration", iteration+1, ":", "y_a", "=", y_a)

#                     # get the variable α
#                     # FIXED
#                     a = model.variables[var_name] 

#                     # ϵy_α+(-M)(1-y_α)≤α≤(+M)y_α
#                     # ϵy_α+(-M)(1-y_α)≤α
#                     c_1 = Constraint(
#                             epsilon * y_a + (-ub) * (1 - y_a) - a,
#                             ub=0
#                         )
#                     # α≤(+M)y_α
#                     c_2 = Constraint(
#                             a - ub * y_a,
#                             ub=0
#                         )
                    
#                     # Add the constraints to the model
#                     self.current_binary_constraints.append(c_1)
#                     self.current_binary_constraints.append(c_2)
#                     model.add(c_1)
#                     model.add(c_2)

#                     if self.stdout_msgs:
#                         print("c_1", c_1)
#                         print("c_2", c_2)

#             # let NZ_a be the cardinality of the set of non-zero variables
#             # ∑_(α∈NZ)▒y_α ≤(card(NZ_α)-1)_
#             NZ_a = 0
#             for var_name, var_primal in self.current_primals.items():
#                 if (var_primal >= epsilon) and (var_name not in set_binary_vars) and (var_name not in self.current_y_binary_variables):
#                     NZ_a += 1
#             # print("self.current_y_binary_variables", self.current_y_binary_variables)
#             # print("self.current_binary_variables", self.current_binary_variables)
#             # print("self.current_variables", self.current_variables)
#             # create sypmy expression of NZ_a - 1 - sum(self.current_y_binary_variables)
#             # expression = sympy.Add(NZ_a - 1, - 
#             #     sympy.Add(*sympy.symbols(self.current_y_binary_variables)))
#             # expression = NZ_a - 1 - sympy.Add(*sympy.symbols(self.current_y_binary_variables))

#             expression = NZ_a - 1

#             if self.stdout_timing:
#                 print(f"Time before accumulating expression using -= y_a: {datetime.datetime.now()}")
#             for var_name, var_primal in self.current_primals.items():
#                 # Didn't work without it, which is weird
#                 # Check if the variable optimal is positive and if it is not
#                 # a binary variable
#                 if var_primal > 0 and var_name not in set_binary_vars:
#                 # if (var_primal >= epsilon) and (var_name not in set_binary_vars):
#                     # get the binary variable y_α
#                     y_a = model.variables[var_name + "_binary" + f"_ya"]
#                     expression -= y_a

#             c4 = Constraint(
#                     expression,
#                     lb=0
#                 )
#             if self.stdout_msgs:
#                 print("c4", c4)
#                 print()
#             model.add(c4)
#             self.current_binary_constraints.append(c4)
           
#             model.objective = \
#                 Objective(
#                     expression=sympy.Add(*sympy.symbols(self.current_variables)),
#                     direction='min'
#                 )
            
#             if self.stdout_timing:
#                 print(f"Time before 'add every binary variable to current_variables': {datetime.datetime.now()}")

#             # add every binary variable to current_variables
#             for var in self.current_binary_variables:
#                 self.current_variables.append(var)

#             for var in self.current_y_binary_variables:
#                 self.current_variables.append(var)
            
#             self.optModel = model  
#             # printing model before optimization in iteration + 1
#             if self.stdout_msgs:
#                 print(f"Model before optimization in iteration {iteration + 1} is:", self.optModel)
#             # print the model before optimization
#             start_time = datetime.datetime.now()
#             if self.stdout_msgs:
#                 print("Optimizing the model in iteration", iteration + 1, "started at:", start_time)
#             status = self.optModel.optimize()
#             end_time = datetime.datetime.now()
#             if self.stdout_msgs:
#                 print(f"Status of the optimization in iteration {iteration + 1} is:", status)
#             with open(f"log Summer 2023.txt", "a") as f:
#                 f.write(f"\nsize = {len(strategies)}, date: {datetime.datetime.now()}, iteration {iteration} with desired nash equilibria: {nasheq_cells} || Time: {end_time - start_time}")

#             # add the objective value to the list of objective values
#             objective_values.append(self.optModel.objective.value)
            

#             # break if the model is infeasible
#             if status == "infeasible":
#                 print("The model became infeasible at this point.")
#                 break
            
#             if self.stdout_msgs:
#                 print("self.current_primals for iteration", iteration + 1, "is:")
#                 pprint(self.current_primals)
#             if self.stdout_msgs:
#                 print("self.optModel.variables for iteration", iteration + 1, "is:")
#                 for var_name, var in self.optModel.variables.items():
#                     print(var_name, "=", var)

#             # set the current_primals to the optimal values
#             if self.stdout_msgs:
#                 print("The current primal values for iteration", iteration + 1, "is:")
#             for var_name, var in self.optModel.variables.items():
#                 if self.stdout_msgs:
#                     print(var_name, "=", var.primal)
#                 self.current_primals[var_name] = var.primal

#             # print the binary variables primal values
#             if self.stdout_msgs:
#                 print("The current binary variables primal values for iteration", iteration + 1, "is:")
#                 for var_name, var in self.optModel.variables.items():
#                     if var_name in self.current_binary_variables or var_name in self.current_y_binary_variables:
#                         print(var_name, "=", var.primal)

#             original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
#             self.validate_2c(nasheq_cells, method=f"New Method 2c iteration {iteration + 1}, epsilon={epsilon}, with new equilibria={nasheq_cells} precluded", iteration=str(iteration), time_spent=end_time - start_time, size=len(strategies))
#             # TODO: is this accurate? maybe not because validate_2c writes
#             # the new payoff matrix that we should use into self.game.payoff_matrix
#             # or yes because we start from the original payoff matrix every time
#             #  find nasheq_cells in payoff_matrix
#             [Nash_equilibria, exit_flag, game_payoff_matrix] = self.optlangRun()
#             print("NE is: ", Nash_equilibria)
#             self.game.payoff_matrix = original_payoff_matrix


# def show_matrix_original(original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells):  
#     fig, ax = mpl.subplots()
#     fig.patch.set_visible(False)
#     ax.axis('off')
#     ax.axis('tight')

#     def payoffs_to_table(payoff_matrix): 
#         col_text =  [''] + strategies
#         table_text = []
#         for i in strategies:
#             row = [i]
#             for j in strategies:
#                 og1 = original_payoff_matrix[(('row', i),('column', j))]['row']
#                 og2 = original_payoff_matrix[(('row', i),('column', j))]['column']
#                 new1 = payoff_matrix[(('row', i),('column', j))]['row']
#                 new2 = payoff_matrix[(('row', i),('column', j))]['column']
#                 def check_0(s):
#                     if s != '0.0':
#                         return f" + {s}"
#                     else:
#                         return ""

#                 def remove_zeroes(s):
#                     if s[-2:] == '.0':
#                         return s[:-2]
#                     else:
#                         return s

#                 row.append(
#                     f"({remove_zeroes(str(round(float(og1), 4)))}" + check_0(str(round(new1 - float(og1), 4))) +
#                     f", {remove_zeroes(str(round(float(og2), 4)))}" + check_0(str(round(new2 - float(og2), 4))) + ")"
#                 )
#             table_text.append(row)
#             # print("dis row", row)

#         return col_text, table_text


#     # print(f"Time before defining the table: {datetime.datetime.now()}")
#     table = payoffs_to_table(payoff_matrix)
#     the_table = ax.table(#colWidths=[0.3] * (SIZE + 1),
#                         cellText=table[1], colLabels=table[0],
#                         loc='center', bbox=[-0.05, 0, 1.1, 1.0],
#                         rowLoc='center', colLoc='center', cellLoc='center')
#     the_table.scale(1, 7)

#     # changed cell is of format (
#     #   (
#     #       (('row', 'S15'), ('column', 'S10')), 'row', 'plus'
#     #   ), 7.009999999999547
#     # )
#     def changed_cell_to_position(changed_cell):
#         return (
#             strategies.index(changed_cell[0][0][0][1])+1,
#             strategies.index(changed_cell[0][0][1][1])+1
#         )
#     # print(f"Time before setting the cells colors: {datetime.datetime.now()}")
#     changed_cells_positions = [changed_cell_to_position(cell) for cell in changed_cells]
#     for (row, col) in changed_cells_positions:
#         the_table[row, col].set_facecolor('#d2905d')


#     def letter_to_position(letter):
#         return strategies.index(letter) + 1
    
#     for eq in range(len(nash_equilibria)):
#         the_table[
#             letter_to_position(nash_equilibria[eq][0][1]), 
#             letter_to_position(nash_equilibria[eq][1][1])
#             ].set_linewidth(3)

#     for (row, col), cell in the_table.get_celld().items():
#         cell.visible_edges = ''
#         if col != 0 and row != 0:
#             cell.visible_edges = 'closed'

#     nasheq_positions = []
#     for eq in range(len(nash_equilibria)):
#         nasheq_positions.append(
#             (letter_to_position(nash_equilibria[eq][0][1]), \
#             letter_to_position(nash_equilibria[eq][1][1]))
#         )

#     size = len(strategies)
#     the_table.auto_set_font_size(False)
#     if size == 2:
#         for row in range(size+1):
#             for col in range(size+1):
#                 if (row, col) not in nasheq_positions:
#                     the_table[row, col].set_fontsize(14)
#                 else:
#                     the_table[row, col].set_fontsize(14)
#     elif size == 10:
#         for row in range(size+1):
#             for col in range(size+1):
#                 if (row, col) not in nasheq_positions:
#                     the_table[row, col].set_fontsize(4)
#                 else:
#                     the_table[row, col].set_fontsize(4)

#     # set title based on iteration
#     # if iteration is in method.split(' ') find the word after it
#     # if iteration is not in method.split(' ') then just use the method

#     if 'iteration' in method.split(' '):
#         title = f"Iteration {method.split(' ')[method.split(' ').index('iteration')+1][:-1]}"
#     else:
#         title = 'Original Game'
    
#     # set title with big font
#     ax.set_title(f"{title}", fontsize=20)
#     mpl.show()
        
    # Elie
    def newEquilibria(self, nasheq_cells, strategies, original_nasheq_cells):
        """
        :param nasheq_cells: a list of elements (('row','C'),('column','C'))
        Consider a payoff value a. We would like to perturb it such that it can
        either increase or decrease. To do this, we define two non-negative 
        variables aa^+ and aa^-. Then we change the payoff as follows:
                                    β=a+aa^+-aa^-
        As the objective function then, you minimize sum of all aa^+'s and 
        aa^-'s for all payoff. 
        return: a model with the solutions' final model, which contains the 
                optimal value of the variables and more
        """   
        
        def tuple_to_var_name(t):
            return str(t).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',','_')
        

        # set tolerance for optlang
        print("optlang.interface is: ", optlang.interface)
        tol = 1e-9
        # # optlang.interface.GUROBI_PARAM_FEASIBILITYTOL = tol

        # self.optModel.configuration.tolerances.feasibility = tol
        # self.optModel.configuration.tolerances.optimality = tol
        # self.optModel.configuration.tolerances.integrality = tol

        # self.optModel.configuration._set_feasibility(tol)
        # self.optModel.configuration._set_optimality(tol)
        # self.optModel.configuration._set_integrality(tol)
        # print("dir(self.optModel): ", dir(self.optModel))
        # print("\nself.optModel.__getstate__()", self.optModel.__getstate__())
        # print("\ndir(self.optModel.configuration): ", dir(self.optModel.configuration))

        # print all callable methods of self.optModel



        # Adding new variables 
        model = Model(name=f'Original Model')
        model.players = self.game.players_names
        indices = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()]
        new_indices = []
        for index in indices:
            for player in model.players:
                new_indices.append(index + (player, 'plus'))
                new_indices.append(index + (player, 'minus'))

        # Format of indices: e.g. ('row','C','column','C','row','plus')
        model.indices = new_indices

        variables_names = []

        # Add the variables
        for index in (model.indices):
            var_name = tuple_to_var_name(index)
            var = Variable(var_name, lb=0, type='continuous', problem=model)
            variables_names.append(var_name)
            model.add(var)

        # Add the objective function SLOW
        model.objective = Objective(expression=sympy.Add(*(sympy.symbols(variables_names))), direction='min')
        # make it faster'
        

        # Add the constraints
        constraints = []
        # Each `cell` is of the format (('row','C'),('column','C'))

        self.b_i_variables = []
        self.b_i_variables_optlang = []
        self.r_c_max_variables = []
        self.r_c_max_variables_optlang = []

        # Add constraints that keep the original Nash equilibrium entries constant
        for cell in (original_nasheq_cells):
            root_index = tuple_to_var_name(cell)
            for player in model.players:
                c1 = Constraint(-model.variables[root_index+'_'+player+'_plus'], lb=0)
                c2 = Constraint(-model.variables[root_index+'_'+player+'_minus'], lb=0)
                constraints.append(c1)
                constraints.append(c2)
                # print("ADDED original_nasheq_cells constraint 1: ", c1)
                # print("ADDED original_nasheq_cells constraint 2: ", c2)
            # c1 = Constraint(-model.variables[root_index+'_'+player+'_plus'], lb=0)
            # c2 = Constraint(-model.variables[root_index+'_'+player+'_minus'], lb=0)
            # constraints.append(c1)
            # constraints.append(c2)
            # print("ADDED original_nasheq_cells constraint 1: ", c1)
            # print("ADDED original_nasheq_cells constraint 2: ", c2)


        # exit()
           
        for cell in (nasheq_cells):
            root_index = tuple_to_var_name(cell)
            # For first player we loop over first axis (rows=strategies)
            for strategy in [x for x in strategies if x != cell[0][1]]:
                current_cell = ((cell[0][0], strategy), (cell[1]))
                current_index = tuple_to_var_name(current_cell)
                player = model.players[0]
                # This new snippet precludes cells other than `nasheq_cells` 
                # to be a Nash equilibrium
                epsilon_nash = 0.01
                c = Constraint(
                        self.game.payoff_matrix[cell][player] \
                        + model.variables[root_index+'_'+player+'_plus'] \
                        - model.variables[root_index+'_'+player+'_minus'] \
                        - 
                        (self.game.payoff_matrix[current_cell][player] \
                        + model.variables[current_index+'_'+player+'_plus'] \
                        - model.variables[current_index+'_'+player+'_minus']
                        )
                        -
                        epsilon_nash, lb=0)
                


                constraints.append(c)

            # For second player we loop over second axis (rows=strategies)
            for strategy in [x for x in strategies if x != cell[1][1]]:
                current_cell = ((cell[0]), (cell[1][0], strategy))
                current_index = tuple_to_var_name(current_cell)
                player = model.players[1]
                c = Constraint(
                        self.game.payoff_matrix[cell][player] \
                        + model.variables[root_index+'_'+player+'_plus'] \
                        - model.variables[root_index+'_'+player+'_minus'] \
                        - 
                        (self.game.payoff_matrix[current_cell][player] \
                        + model.variables[current_index+'_'+player+'_plus'] \
                        - model.variables[current_index+'_'+player+'_minus']
                        )
                        -
                        epsilon_nash, lb=0)
                
                constraints.append(c)
            # print("ADDED columns loop")

        # b_ij^k=a_ij^k+〖α_ij^k〗^+-〖α_ij^k〗^-,
        # rmax_rc≥b_ic^1,	∀(r,c)∈U,∀i∈S_1-{i|(i,c)∈U},	(4)
        # rmax_rc≥b_rc^1,+ϵ	∀(r,c)∈U	(5)
        # cmax_rc≥b_ri^2,	∀(r,c)∈U,∀i∈S_2-{i|(r,i)∈U},	(6)
        # cmax_rc≥b_rc^2,+ϵ	∀(r,c)∈U	(7)
        # rmax and cmax representing the maximum payoff of the Player 1 when 
        # Player 2 takes strategy corresponding to c and cmax is the maximim 
        # payoff of Player 2 when Player 1 takes the strategy corresponding 
        # to r. Constraint (4) requires rmax_rc to be greater than all payoffs 
        # in (i,c) except for the row r and Constraint (5) requires rmax_rc 
        # to be strictly greater than the payoff of Player 1 in cell (r,c). 
        # Constraints (6) and (8) impose similar constriants for cmax_rc 
        # and the payoff of Player 2in cell (r,c).
        # U = original_nasheq_cells
        # D = desired_nasheq_cells
        # S_1 = strategies of player 1
        # S_2 = strategies of player 2
        for cell in original_nasheq_cells:
            root_index = tuple_to_var_name(cell)
            # constraint (5): # rmax_rc≥b_rc^1,+ϵ
            # get rmax_rc
            # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] for r in strategies if r != cell[0][1]])
            # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] for r in strategies])
            # rmax_rc = max([self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] + model.variables[tuple_to_var_name((('row', r), ('column', cell[1][1])))+'_row_plus'] - model.variables[tuple_to_var_name((('row', r), ('column', cell[1][1])))+'_row_minus'] for r in strategies])

            # https://math.stackexchange.com/questions/2446606/linear-programming-set-a-variable-the-max-between-two-another-variables
            # \begin{align}
            # U &\ge a_i   &\forall i \in N \\ 
            # U &\le a_i + (1-b_i)*M  & \forall i \in N \\
            # \sum_{i \in N} b_i &= 1
            # \end{align}
            # where the b_i\in{0,1} is a binary variable that indicates the maximum a_i
            # (i.e. b_i=1 when a_i is the max value), and M it's a "big number".

            # create the b_i variables (bi1 and bi2)
            b_i1_variables = []
            b_i1_variables_optlang = []
            b_i2_variables = []
            b_i2_variables_optlang = []
            for r in strategies:
                # create the b_i1 variables
                str_index = root_index + "_binary1_preculde_original_nash_equilibria_" + r
                b_i1_variables.append(str_index)
                var = Variable(str_index, lb=0, ub=1, type='binary', problem=model)
                b_i1_variables_optlang.append(var)
                model.add(var)
                # create the b_i2 variables
                str_index = root_index + "_binary2_preculde_original_nash_equilibria_" + r
                b_i2_variables.append(str_index)
                var = Variable(str_index, lb=0, ub=1, type='binary', problem=model)
                b_i2_variables_optlang.append(var)
                model.add(var)



            for binary_var in b_i1_variables:
                self.b_i_variables.append(binary_var)
            for binary_var in b_i1_variables_optlang:
                self.b_i_variables_optlang.append(binary_var)
            for binary_var in b_i2_variables:
                self.b_i_variables.append(binary_var)
            for binary_var in b_i2_variables_optlang:
                self.b_i_variables_optlang.append(binary_var)

            # define U, and call it root_index + "_rmac_rc"
            str_rmax_rc = root_index + "_rmax_rc"
            rmax_rc = Variable(str_rmax_rc, type='continuous', problem=model)
            model.add(rmax_rc)
            self.r_c_max_variables.append(str_rmax_rc)
            self.r_c_max_variables_optlang.append(rmax_rc)
            # create the constraints for rmax_rc
            # U &\ge a_i   &\forall i \in N \\
            # U &\le a_i + (1-b_i)*M  & \forall i \in N \\
            for r in strategies:
                ic_cell = ((('row', r), cell[1]))
                ic_index = tuple_to_var_name(ic_cell)
                c = Constraint(
                        rmax_rc \
                        - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
                        - model.variables[ic_index + '_row_plus'] \
                        + model.variables[ic_index + '_row_minus'] \
                        , lb=0)
                constraints.append(c)
                
                c = Constraint(
                        rmax_rc \
                        - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
                        - model.variables[ic_index + '_row_plus'] \
                        + model.variables[ic_index + '_row_minus'] \
                        - (1 - model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r]) * 1000 \
                        , ub=0)
                constraints.append(c)
            # \sum_{i \in N} b_i1 <= 1
            c = Constraint(
                    sum([model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r] for r in strategies]) \
                    , ub=1)
            constraints.append(c)

            # constraint (4): # rmax_rc≥b_ic^1,	∀(r,c)∈U,∀i∈S_1-{i|(i,c)∈U}
            for r in strategies:
                if r != cell[0][1]:
                    ic_cell = ((('row', r), cell[1]))
                    ic_index = tuple_to_var_name(ic_cell)
                    c = Constraint(
                            rmax_rc \
                            - self.game.payoff_matrix[(('row', r), ('column', cell[1][1]))]['row'] \
                            - model.variables[ic_index + '_row_plus'] \
                            + model.variables[ic_index + '_row_minus'] \
                            , lb=0)
                    constraints.append(c)


            # constraint (5): # rmax_rc≥b_rc^1,+ϵ
            c = Constraint(
                    rmax_rc \
                    - self.game.payoff_matrix[cell]['row'] \
                    - model.variables[root_index+'_row_plus'] \
                    + model.variables[root_index+'_row_minus'] \
                    - epsilon_nash, lb=0)
            constraints.append(c)


            

            # define U, call it root_index + "_cmax_rc"
            str_cmax_rc = root_index + "_cmax_rc"
            cmax_rc = Variable(str_cmax_rc, type='continuous', problem=model)
            model.add(cmax_rc)
            self.r_c_max_variables.append(str_cmax_rc)
            self.r_c_max_variables_optlang.append(cmax_rc)
            # create the constraints for cmax_rc
            for strat in strategies:
                ri_cell = ((cell[0]), (('column', strat)))
                ri_index = tuple_to_var_name(ri_cell)
                c = Constraint(
                        cmax_rc \
                        - self.game.payoff_matrix[(('row', cell[0][1]), ('column', strat))]['column'] \
                        - model.variables[ri_index + '_column_plus'] \
                        + model.variables[ri_index + '_column_minus'] \
                        , lb=0)
                constraints.append(c)

                c = Constraint(
                        cmax_rc \
                        - self.game.payoff_matrix[(('row', cell[0][1]), ('column', strat))]['column'] \
                        - model.variables[ri_index + '_column_plus'] \
                        + model.variables[ri_index + '_column_minus'] \
                        - (1 - model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat]) * 1000 \
                        , ub=0)
                constraints.append(c)
            # \sum_{i \in N} b_i2 <= 1
            c = Constraint(
                    sum([model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat] for strat in strategies]) \
                    , ub=1)
            constraints.append(c)

            # constraint (6): # cmax_rc≥b_ri^2,	∀(r,c)∈U,∀i∈S_2-{i|(r,i)∈U}
            for c in strategies:
                if c != cell[1][1]:
                    ri_cell = ((cell[0]), (('column', c)))
                    ri_index = tuple_to_var_name(ri_cell)
                    c = Constraint(
                            cmax_rc \
                            - self.game.payoff_matrix[(('row', cell[0][1]), ('column', c))]['column'] \
                            - model.variables[ri_index + '_column_plus'] \
                            + model.variables[ri_index + '_column_minus'] \
                            , lb=0)
                    constraints.append(c)


            # constraint (7): # cmax_rc≥b_rc^2,+ϵ
            c = Constraint(
                    cmax_rc \
                    - self.game.payoff_matrix[cell]['column'] \
                    - model.variables[root_index+'_column_plus'] \
                    + model.variables[root_index+'_column_minus'] \
                    - epsilon_nash, lb=0)
            constraints.append(c)



            # constraint (10 from 12/29/2023): \sum_{i \in N} b_i1 + \sum_{i \in N} b_i2 >= 1
            c = Constraint(
                    sum([model.variables[root_index + "_binary1_preculde_original_nash_equilibria_" + r] for r in strategies]) \
                    + sum([model.variables[root_index + "_binary2_preculde_original_nash_equilibria_" + strat] for strat in strategies]) \
                    , lb=1)
            constraints.append(c)
                    


        
        # print(constraints)
        model.add(constraints)  
       
        self.optModel = model  
        start_time = datetime.datetime.now()
        self.optModel.optimize()
        end_time = datetime.datetime.now() 

        # print the final value of cells (5,2) and (5, 3)
        # print("final row value of cell (5,2) is: ", self.optModel.variables['row_S5_column_S2_row_plus'].primal - self.optModel.variables['row_S5_column_S2_row_minus'].primal + self.game.payoff_matrix[(('row', 'S5'), ('column', 'S2'))]['row'])
        # print("final col value of cell (5,3) is: ", self.optModel.variables['row_S5_column_S3_column_plus'].primal - self.optModel.variables['row_S5_column_S3_column_minus'].primal + self.game.payoff_matrix[(('row', 'S5'), ('column', 'S3'))]['column'])
        # print("final row value of cell (5,3) is: ", self.optModel.variables['row_S5_column_S3_row_plus'].primal - self.optModel.variables['row_S5_column_S3_row_minus'].primal + self.game.payoff_matrix[(('row', 'S5'), ('column', 'S3'))]['row'])
        # print("final col value of cell (5,2) is: ", self.optModel.variables['row_S5_column_S2_column_plus'].primal - self.optModel.variables['row_S5_column_S2_column_minus'].primal + self.game.payoff_matrix[(('row', 'S5'), ('column', 'S2'))]['column'])

        # print the final value of cells (1,1), (1,2), (2,1) and (2,2)
        # print("final row value of cell (1,1) is: ", self.optModel.variables['row_S1_column_S1_row_plus'].primal - self.optModel.variables['row_S1_column_S1_row_minus'].primal + self.game.payoff_matrix[(('row', 'S1'), ('column', 'S1'))]['row'])
        # print("final col value of cell (1,1) is: ", self.optModel.variables['row_S1_column_S1_column_plus'].primal - self.optModel.variables['row_S1_column_S1_column_minus'].primal + self.game.payoff_matrix[(('row', 'S1'), ('column', 'S1'))]['column'])
        # print("final row value of cell (1,2) is: ", self.optModel.variables['row_S1_column_S2_row_plus'].primal - self.optModel.variables['row_S1_column_S2_row_minus'].primal + self.game.payoff_matrix[(('row', 'S1'), ('column', 'S2'))]['row'])
        # print("final col value of cell (1,2) is: ", self.optModel.variables['row_S1_column_S2_column_plus'].primal - self.optModel.variables['row_S1_column_S2_column_minus'].primal + self.game.payoff_matrix[(('row', 'S1'), ('column', 'S2'))]['column'])
        # print("final row value of cell (2,1) is: ", self.optModel.variables['row_S2_column_S1_row_plus'].primal - self.optModel.variables['row_S2_column_S1_row_minus'].primal + self.game.payoff_matrix[(('row', 'S2'), ('column', 'S1'))]['row'])
        # print("final col value of cell (2,1) is: ", self.optModel.variables['row_S2_column_S1_column_plus'].primal - self.optModel.variables['row_S2_column_S1_column_minus'].primal + self.game.payoff_matrix[(('row', 'S2'), ('column', 'S1'))]['column'])
        # print("final row value of cell (2,2) is: ", self.optModel.variables['row_S2_column_S2_row_plus'].primal - self.optModel.variables['row_S2_column_S2_row_minus'].primal + self.game.payoff_matrix[(('row', 'S2'), ('column', 'S2'))]['row'])
        # print("final col value of cell (2,2) is: ", self.optModel.variables['row_S2_column_S2_column_plus'].primal - self.optModel.variables['row_S2_column_S2_column_minus'].primal + self.game.payoff_matrix[(('row', 'S2'), ('column', 'S2'))]['column'])
        # for var_name, var in self.optModel.variables.items():
        #     # if var.primal != 0:
        #     print(var_name, "=", var.primal)





        # Neccessary for validate_2c
        self.current_variables = copy.deepcopy(variables_names)
        self.current_primals = {}
        
        for var_name, var in (self.optModel.variables.items()):
            # print(var_name, "=", var.primal)
            self.current_variables.append(var_name)
            self.current_primals[var_name] = var.primal

        original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        # for validate_2c to work
        self.current_binary_variables = []
        self.current_binary_constraints = []
        self.current_y_binary_variables = []
        self.current_y_binary_variables_optlang = []
    

        self.validate_2c(nasheq_cells, method=f'First solution with new equilibria={nasheq_cells} precluded', iteration="NONE", time_spent=end_time-start_time, size=len(strategies))
        # [NE, exit_flag, g_payoff_matrix] = self.optlangRun()
        # print("matrix for original solution TOUT is: ", g_payoff_matrix)
        # print("Nash_equilibria for original solution TOUT is: ", NE)
        self.game.payoff_matrix = original_payoff_matrix

        # New on May 19
        # Iterative Method 2c

        # For each non-zero α whose optimal value is α^opt, add the following
        # constraints:
        #                      α≤α^opt-ϵ  &  α≥ α^opt+ϵ

        #                     α≤〖(1-y)(α〗^opt-ϵ)+y〖UB〗_α    
        #                     〖α≥y(α〗^opt+ϵ)+(1-y)LB_α 

        # Note that if y=0, the first constraint is reduced to α≤α^opt-ϵ and 
        # the second constraint is reduced to α≥〖LB〗_α. On the other hand, 
        # if y=1, the first constraint is reduced to α≤UB_α and the second 
        # constraint is reducd to 〖α≥α〗^opt+ϵ. Instead of LB_α and UB_α, 
        # you can use the big-M approach too, i.e., replaced them with -M and 
        # M, respectively. 


      
        self.current_variables = copy.deepcopy(variables_names)
        self.current_binary_variables = []
        self.current_binary_constraints = []
        self.current_y_binary_variables = []
        self.current_y_binary_variables_optlang = []
        self.current_primals = {}
        for var_name, var in (self.optModel.variables.items()):
            # print(var_name, "=", var.primal)
            # TODO: Check if this is needed
            self.current_variables.append(var_name)
            # TODO: Check if this is correct because validate_2c is using this.
            #       Do we need primals for variables that are not in the original
            #       payoff matrix?
            self.current_primals[var_name] = var.primal

        # adding a binary variable for each α

        for var_name, var_primal in (self.current_primals.items()):
            str_index_y_a = var_name + "_binary" + f"_ya" # no iteration needed in this method 
            y_a = Variable(str_index_y_a, lb=0, type='binary', problem=model)
            model.add(y_a)
            self.current_binary_variables.append(str_index_y_a)
            # Add the binary y variables to the list of binary variables
            self.current_y_binary_variables.append(str_index_y_a)
            self.current_y_binary_variables_optlang.append(y_a)


       
        objective_values = []

        # For iteration = 1 to 5
        for iteration in range(50):
            # print the payoff matrix
            # Elie on 1/27/2023
            print(f"\n\n\n----iteration {iteration + 1}----\n\n\n")
            # print the payoff matrix


            set_binary_vars = set(self.current_binary_variables + self.current_y_binary_variables + self.b_i_variables + self.r_c_max_variables)

            # print(f"----removing binary variables and constrains to start anew----\n")


            self.current_variables = [e for e in self.current_variables if e not in set_binary_vars]
          

            set_current_primals = set(self.current_primals.keys())
            for key_to_remove in (self.current_binary_variables):
                if key_to_remove in set_current_primals:
                    del self.current_primals[key_to_remove]
                # del self.current_primals[key_to_remove]
            # On june 20
            epsilon = 0.01
            ub = 1000
            lb = -1000 

           
            for var_name, var_primal in (self.current_primals.items()):
                # Didn't work without it, which is weird
                # Check if the variable optimal is positive and if it is not
                # a binary variable
                if var_primal >= epsilon and var_name not in set_binary_vars:
                    # I will introduce two non-linear constraints:
                    # Define a binary variable y_α such that if y_α=1, then α≠0
                    # and if y=0, then α=0 in the previous solution (note that 
                    # here α could be either α^+ or α^-). This can be imposed 
                    # using the following cosntrinat
                    #               ϵy_α+(-M)(1-y_α)≤α≤(+M)y_α
                    # Note that according to this constaint: 
                    # y_α=1→α≠0
                    # y_α=0→α=0

                    # Keep in mind that:
                    #     Here, ϵ is a small positive value (e.g., 10e-5)
                    #     Here we define binary variables for each α^+ or α^- 
                    #     only once, i..e, we do not define new binary variales 
                    #     in each iteration. 
                    # For the next solution, we want at least one perturbed 
                    # payoff value to be different from the solution that we 
                    # had before. We can easily impose this by using the 
                    # following constraint:
                    #              ∑_(α∈NZ)▒y_α ≤(card(NZ_α)-1)_+

                    # if self.stdout_msgs:
                    #     print("PIPTH in iteration", iteration+1, ":", var_name, "=", var_primal)
                    #     print("var_primal", var_primal)

                    # get the binary variable y_α
                    y_a = model.variables[var_name + "_binary" + f"_ya"]
                   
                    # get the variable α
                    # FIXED

                    a = model.variables[var_name] 

                    # ϵy_α+(-M)(1-y_α)≤α≤(+M)y_α
                    # ϵy_α+(-M)(1-y_α)≤α
                    c_1 = Constraint(
                            epsilon * y_a + (-ub) * (1 - y_a) - a,
                            ub=0
                        )
                    # α≤(+M)y_α
                    c_2 = Constraint(
                            a - ub * y_a,
                            ub=0
                        )
                    
                    # Add the constraints to the model
                    self.current_binary_constraints.append(c_1)
                    self.current_binary_constraints.append(c_2)
                    model.add(c_1)
                    model.add(c_2)


            # let NZ_a be the cardinality of the set of non-zero variables
            # ∑_(α∈NZ)▒y_α ≤(card(NZ_α)-1)_
            NZ_a = 0

            for var_name, var_primal in (self.current_primals.items()):
                if (var_primal >= epsilon) and (var_name not in set_binary_vars) and (var_name not in self.current_y_binary_variables):
                    print(var_name, var_primal)
                    NZ_a += 1

            # print("self.current_y_binary_variables", self.current_y_binary_variables)
            # print("self.current_binary_variables", self.current_binary_variables)
            # print("self.current_variables", self.current_variables)
            # create sypmy expression of NZ_a - 1 - sum(self.current_y_binary_variables)
            # expression = sympy.Add(NZ_a - 1, - 
            #     sympy.Add(*sympy.symbols(self.current_y_binary_variables)))
            # expression = NZ_a - 1 - sympy.Add(*sympy.symbols(self.current_y_binary_variables))

            expression = NZ_a - 1
            len_expression = 0

            
            for var_name, var_primal in (self.current_primals.items()):
                # Didn't work without it, which is weird
                # Check if the variable optimal is positive and if it is not
                # a binary variable
                # if var_primal > 0 and var_name not in set_binary_vars:
                if (var_primal >= epsilon) and (var_name not in set_binary_vars):
                    # get the binary variable y_α
                    y_a = model.variables[var_name + "_binary" + f"_ya"]
                    expression -= y_a
                    len_expression += 1
                # check
                    


            c4 = Constraint(
                    expression,
                    lb=0
                )

            model.add(c4)
            self.current_binary_constraints.append(c4)
            # print("current_binary_constraints")
            # for c in self.current_binary_constraints:
            #     print(c)
            # sympy.Add(*sympy.symbols(self.current_variables))

            # print(sympy.Add(*sympy.symbols(self.current_variables)))
            # print(sympy.symbols(self.current_y_binary_variables))


            model.objective = \
                Objective(
                    expression=sympy.Add(*sympy.symbols(self.current_variables)),
                    direction='min'
                )
            

            # add every binary variable to current_variables
            for var in self.current_binary_variables:
                self.current_variables.append(var)

            for var in self.current_y_binary_variables:
                self.current_variables.append(var)
            
            self.optModel = model  
           
            # print the model before optimization
            start_time = datetime.datetime.now()
           
            status = self.optModel.optimize()
            end_time = datetime.datetime.now()
           
            # for i in range(5):
            #     start_time = datetime.datetime.now()
            #     status = self.optModel.optimize()
            #     end_time = datetime.datetime.now()
            #     with open(f"log Summer 2023.txt", "a") as f:
            #         f.write(f"\nsize = {len(strategies)}, repeated iteration {iteration} with desired nash equilibria: {nasheq_cells} || Time: {end_time - start_time}")


            # add the objective value to the list of objective values
            objective_values.append(self.optModel.objective.value)
            

            # break if the model is infeasible
            if status == "infeasible":
                print("The model is infeasible")
                # [Nash_equilibria, exit_flag, game_payoff_matrix] = self.optlangRun()
                # print("matrix TOUT is: ", game_payoff_matrix)
                # print("Nash_equilibria TOUT is: ", Nash_equilibria)
                break
            
          

            # set the current_primals to the optimal values
            
            for var_name, var in self.optModel.variables.items():
                self.current_primals[var_name] = var.primal

            # print the binary variables primal values
          
            original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
            self.validate_2c(nasheq_cells, method=f"New Method 2c iteration {iteration + 1}, epsilon={epsilon}, with new equilibria={nasheq_cells} precluded", iteration=str(iteration), time_spent=end_time - start_time, size=len(strategies))

            # TODO: is this accurate? maybe not because validate_2c writes
            # the new payoff matrix that we should use into self.game.payoff_matrix
            # or yes because we start from the original payoff matrix every time
            #  find nasheq_cells in payoff_matrix
            [Nash_equilibria, exit_flag, game_payoff_matrix] = self.optlangRun()
            # print("matrix TOUT is: ", game_payoff_matrix)

            self.game.payoff_matrix = original_payoff_matrix

        # print the objective values
        # if self.stdout_msgs:

def show_matrix_original(original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells):  
    fig, ax = mpl.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')

    def payoffs_to_table(payoff_matrix): 
        col_text =  [''] + strategies
        table_text = []
        for i in strategies:
            row = [i]
            for j in strategies:
                og1 = original_payoff_matrix[(('row', i),('column', j))]['row']
                og2 = original_payoff_matrix[(('row', i),('column', j))]['column']
                new1 = payoff_matrix[(('row', i),('column', j))]['row']
                new2 = payoff_matrix[(('row', i),('column', j))]['column']
                def check_0(s):
                    if s != '0.0':
                        return f" + {s}"
                    else:
                        return ""

                def remove_zeroes(s):
                    if s[-2:] == '.0':
                        return s[:-2]
                    else:
                        return s

                row.append(
                    f"({remove_zeroes(str(round(float(og1), 4)))}" + check_0(str(round(new1 - float(og1), 4))) +
                    f", {remove_zeroes(str(round(float(og2), 4)))}" + check_0(str(round(new2 - float(og2), 4))) + ")"
                )
            table_text.append(row)
            # print("dis row", row)

        return col_text, table_text


    # print(f"Time before defining the table: {datetime.datetime.now()}")
    table = payoffs_to_table(payoff_matrix)
    the_table = ax.table(#colWidths=[0.3] * (SIZE + 1),
                        cellText=table[1], colLabels=table[0],
                        loc='center', bbox=[-0.05, 0, 1.1, 1.0],
                        rowLoc='center', colLoc='center', cellLoc='center')
    the_table.scale(1, 7)

    # changed cell is of format (
    #   (
    #       (('row', 'S15'), ('column', 'S10')), 'row', 'plus'
    #   ), 7.009999999999547
    # )
    def changed_cell_to_position(changed_cell):
        return (
            strategies.index(changed_cell[0][0][0][1])+1,
            strategies.index(changed_cell[0][0][1][1])+1
        )
    # print(f"Time before setting the cells colors: {datetime.datetime.now()}")
    changed_cells_positions = [changed_cell_to_position(cell) for cell in changed_cells]
    for (row, col) in changed_cells_positions:
        the_table[row, col].set_facecolor('#d2905d')


    def letter_to_position(letter):
        return strategies.index(letter) + 1
    
    for eq in range(len(nash_equilibria)):
        the_table[
            letter_to_position(nash_equilibria[eq][0][1]), 
            letter_to_position(nash_equilibria[eq][1][1])
            ].set_linewidth(3)

    for (row, col), cell in the_table.get_celld().items():
        cell.visible_edges = ''
        if col != 0 and row != 0:
            cell.visible_edges = 'closed'

    nasheq_positions = []
    for eq in range(len(nash_equilibria)):
        nasheq_positions.append(
            (letter_to_position(nash_equilibria[eq][0][1]), \
            letter_to_position(nash_equilibria[eq][1][1]))
        )

    size = len(strategies)
    the_table.auto_set_font_size(False)
    if size == 2:
        for row in range(size+1):
            for col in range(size+1):
                if (row, col) not in nasheq_positions:
                    the_table[row, col].set_fontsize(14)
                else:
                    the_table[row, col].set_fontsize(14)
    elif size == 10:
        for row in range(size+1):
            for col in range(size+1):
                if (row, col) not in nasheq_positions:
                    the_table[row, col].set_fontsize(4)
                else:
                    the_table[row, col].set_fontsize(4)

    # set title based on iteration
    # if iteration is in method.split(' ') find the word after it
    # if iteration is not in method.split(' ') then just use the method

    if 'iteration' in method.split(' '):
        title = f"Iteration {method.split(' ')[method.split(' ').index('iteration')+1][:-1]}"
    else:
        title = 'Original Game'
    
    # set title with big font
    ax.set_title(f"{title}", fontsize=20)
    mpl.show()