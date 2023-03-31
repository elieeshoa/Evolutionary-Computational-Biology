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
from pyomo.environ import *
from pyomo.opt import *
from sympy.sets.sets import FiniteSet
sys.path.append('/Users/elieeshoa/Dropbox/Elie_Eshoa/Ali_codes/')
from pyomoSolverCreator import *
import optlang # 1.5.2
print(optlang.__version__)
print(optlang.__file__)

import sympy # 1.8
print(sympy.__version__)
# exit()
import matplotlib.pyplot as mpl

# The following lines change the temporary directory for pyomo
# from pyutilib.services import TempfileManager
# TempfileManager.tempdir = pyomo_tmp_dir


SIZE = 2



class NashEqFinder(object):
    """
    General class for NashEq Finder. Sample usage is provided at the end 
    """   

    def __init__(self, game, NashEq_type = 'pure', optimization_solver = 'gurobi', warnings = True, stdout_msgs = True, output_file = ''):
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
        of the game, i.e., ('p1','s1','p2','s2') is converted to (('p1','s1'),('p2','s2')) 
        (see optModel.I for details)
        """
        gameState = []
        done = 0
        k1 = list(i)
        while done == 0:
            gameState.append(tuple(k1[0:2]))
            # print('tuple(k1[0:2])', tuple(k1[0:2]))
            del k1[0:2]
            if len(k1) == 0:
                done = 1
        return tuple(gameState)

        
    def createPyomoModel(self):
        """
        This creates a pyomo optimization model 

        Instead of several indicies for binary variables (y), we just define a single set I containing all
        possible labels of the payoff matrix (combinations of players and strategies). 
        """   
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of players
        optModel.P = Set(initialize = self.game.players_names) 

        # Set of players' strategy combinations 
        # Keys of the game.payoff_matrix are in the form of a list of tuples, where each
        # tuple is compased of inner tuple of length two, e.g., 
        # [(('p1','s1'),('p2','s2')),(('p1','s2'),('p2','s1')),...]
        # These keys should serve as the elements of the set I in the optimization model,
        # however, pyomo does not accept list of tuples with nested tuples. Therefore, we 
        # need to convert this to a list of tuples with no inner tuples, i.e.,
        # [('p1','s1','p2','s2'),('p1','s2','p2','s1'),...]
        
        # optlangInit = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()]
        optModel.I = Set(initialize = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()])   

        #--- Define the variables --- 
        optModel.y = Var(optModel.I, domain=Boolean)

        #--- Define the objective function and constraints ----
        # Objective function
        optModel.objective_rule = Objective(rule = lambda optModel: sum(optModel.y[i] for i in optModel.I), sense = maximize)

        # Constraint checking the best strategy of player p given the strategy of 
        # all other players 
        def NashCond_rule(optModel,p,*i):

            # Convert the game state to the format of keys of the payoff matrix
            i = self.convert_to_payoffMatrix_key(i)

            # All possible responses of P to the action all other players
            # have taken in i
            responseP = [k for k in self.game.payoff_matrix.keys() if False not in [dict(k)[pp] == dict(i)[pp] for  pp in dict(i).keys() if pp != p]] 

            # Find the payoff of the best response of player P 
            bestResP = max([self.game.payoff_matrix[k][p] for k in responseP])

            return self.game.payoff_matrix[i][p] >= bestResP*optModel.y[i] + self.payoffLB*(1 - optModel.y[i])

        optModel.NashCond = Constraint(optModel.P,optModel.I, rule=NashCond_rule)

        self.optModel = optModel 
        

    # Elie
    def createOptlangModel(self):
        """
        This creates a optlang optimization model 
        """   
        
        def add_optlang_NashCond_rule(optlangOptModel,p,i):

            # Convert the game state to the format of keys of the payoff matrix
            i = self.convert_to_payoffMatrix_key(i)
            
            responseP = [k for k in self.game.payoff_matrix.keys() if False not in [dict(k)[pp] == dict(i)[pp] for  pp in dict(i).keys() if pp != p]]
            # print('optlang responseP', responseP)

            # Find the payoff of the best response of player P 
            bestResP = max([self.game.payoff_matrix[k][p] for k in responseP])

            model.add(
                optlang.Constraint(
                    self.game.payoff_matrix[i][p] - \
                      bestResP*optlangOptModel.variables[str(i).replace(" ", "")\
                      .replace('(', "").replace(')', "").replace("'","")\
                      .replace(',','_')] - \
                        self.payoffLB * (1 - \
                        optlangOptModel.variables[str(i).replace(" ", "").\
                        replace('(', "").replace(')', "").replace("'","").\
                        replace(',','_')]),
                    lb=0))

      
        model = optlang.Model(name='Original Model')
        # exit()
        model.players = self.game.players_names
        model.indices = [tuple([k3 for k2 in k1 for k3 in k2]) for k1 in self.game.payoff_matrix.keys()]

        variables_names = []

        # Add the variables
        for index in model.indices:
            str_index = str(index).replace(" ", "").replace('(', "").replace(')', "").replace("'","")
            var = optlang.Variable(str_index.replace(',', '_'), type='binary', problem=model)
            model.add(var)
            variables_names.append(str_index.replace(',', '_'))

        # Add the objective function
        model.objective = \
            optlang.Objective(
                expression=sympy.Add(*sympy.symbols(variables_names)), \
                direction='max')
   
        # Add the constraints
        for player in model.players:
            for index in model.indices:
                add_optlang_NashCond_rule(model, player, index)

        # print('FINAL optlang model', model)
        
        self.optModel = model
    
    
    def findPure(self):
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
        # Processing and wall time
        # Elie
        # start_run_pt = time.clock()
        start_run_pt = time.process_time()
        start_run_wt = time.time()

        #---- Creating and instantiating the optModel ----
        # Elie
        # start_pyomo_pt = time.clock()
        start_pyomo_pt = time.process_time()
        start_pyomo_wt = time.time()

        # Create the optModel model        
        self.createPyomoModel()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optimization_solver)

        # Elie
        # elapsed_pyomo_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
        elapsed_pyomo_pt = str(timedelta(seconds = time.process_time() - start_pyomo_pt))
        elapsed_pyomo_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))

        #-- Some initializations --
        # Instantiate the optModel with new fixed variables
        self.optModel.preprocess()

        #- Solve the optModel (tee=True shows the solver output) -
        try:
            # Elie
            # start_solver_pt = time.clock()
            start_solver_pt = time.process_time()
            start_solver_wt = time.time()

            optSoln = solverType.solve(self.optModel,tee=False)
            solverFlag = 'normal'
    
        # In the case of an error switch the solver
        except:
            if self.warnings:
                print ("WARNING! ",self.optimization_solver," failed. An alternative solver is tried")  
    
            if self.optimization_solver.lower() == 'gurobi':
                self.optimization_solver = 'cplex'
            elif self.optimization_solver.lower() == 'cplex':
                self.optimization_solver = 'gurobi'
    
            # Try solving with the alternative solver
            solverType = pyomoSolverCreator(self.optimization_solver)
            try:
                # Elie
                # start_solver_pt = time.clock()
                start_solver_pt = time.process_time()
                start_solver_wt = time.time()

                optSoln = solverType.solve(self.optModel,tee=False)
                solverFlag = 'normal'
            except:
                solverFlag = 'solverError'
                if self.warnings:
                    print ('\nWARNING! The alternative solver failed. No solution was returned')
        # Elie
        # elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_pt = str(timedelta(seconds = time.process_time() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))
    
        #----- Print the results in the output (screen, file and/or variable) ------
        # Load the results (model.load() is dprecated)
        #self.optModel.load(optSoln)
            
        # Set of the Nash equilibria
        self.Nash_equilibria = []
        
        if solverFlag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':
            
            optimExitflag = 'globallyOptimal'
    
            # Value of the objective function
            objValue = self.optModel.objective_rule()
    
            # Print the results on the screen 
            if self.stdout_msgs:
                print ("\nsolver.status = ",optSoln.solver.termination_condition,"\n")
                print ("objective value = ",objValue)

            if objValue >= 1:
                self.exit_flag = 'objGreaterThanZero'
                for i in self.optModel.I.value: 
                    if self.optModel.y[i].value == 1:
                        self.Nash_equilibria.append(list(self.convert_to_payoffMatrix_key(i)))
            elif objValue == 0:
                done = 1
                self.exit_flag = 'objIsZero'
                      
            # Write the results into the output file 
            if self.output_file != '': 
                pass   # To be added 

        # If the optimization problem was not solved successfully
        else:

            if solverFlag == 'solverError':
                optimExitflag = solverFlag
                self.exit_flag = solverFlag
            else:
                optimExitflag = str(optSoln.solver.termination_condition)
                self.exit_flag = str(optSoln.solver.termination_condition)
 
            objValue = None 
    
            # Write on the screen
            if self.warnings:
                print ("\nWARNING! No optimal solutions found (solution.solver.status = ",optSoln.Solution.status,", solver.status =",optSoln.solver.status,", solver.termination_condition = ",optSoln.solver.termination_condition,")\n")
    
            # Write the results into the output file
            if self.output_file != None: 
                pass    # *** To be completed ***
            else:
                pass
    
        # Time required to run 
        # Elie
        # elapsed_run_pt = str(timedelta(seconds = time.clock() - start_run_pt))
        elapsed_run_pt = str(timedelta(seconds = time.process_time() - start_run_pt))
        elapsed_run_wt = str(timedelta(seconds = time.time() - start_run_wt))
    
        if self.stdout_msgs:
           print ('NashEqFinder took (hh:mm:ss) (processing/wall) time: pyomo = {}/{}  ,  solver = {}/{}  ,  run = {}/{} for a game with {} cells in its payoff matrix\n'.format(elapsed_pyomo_pt,elapsed_pyomo_wt,elapsed_solver_pt,elapsed_solver_wt,elapsed_run_pt,elapsed_run_wt, len(self.game.payoff_matrix)) )

    

    # Elie
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
        # Processing and wall time
        start_run_pt = time.process_time()
        start_run_wt = time.time()

        #---- Creating and instantiating the optModel ----
        start_optlang_pt = time.process_time()
        start_optlang_wt = time.time()

        # Create the optModel model    
        self.createOptlangModel()
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

        # Print the results on the screen 
        # print("status:", self.optModel.status)
        # print("objective value:", self.optModel.objective.value)
        # print("--------------------------------------------------")
        # for var_name, var in self.optModel.variables.items():
        #     print(var_name, "=", var.primal)
        solverFlag = 'normal'
    
        # Elie
        elapsed_solver_pt = \
            str(timedelta(seconds = time.process_time() - start_solver_pt))
        elapsed_solver_wt = \
            str(timedelta(seconds = time.time() - start_solver_wt))
    
            
        # Set of the Nash equilibria
        self.Nash_equilibria = []
        objValue = self.optModel.objective.value
        if objValue >= 1:
            self.exit_flag = 'objGreaterThanZero'
            for i in self.optModel.indices: 
                if self.optModel.variables[
                    str(i).replace(" ", "").replace('(', "").replace(')', "").\
                    replace("'","").replace(',', '_')]\
                    .primal == 1:
                    self.Nash_equilibria.append(list(self.convert_to_payoffMatrix_key(i)))
        elif objValue == 0:
            done = 1
            self.exit_flag = 'objIsZero'
    
        # Time required to run 
        elapsed_run_pt = str(timedelta(seconds = time.process_time() - start_run_pt))
        elapsed_run_wt = str(timedelta(seconds = time.time() - start_run_wt))
    
        if self.stdout_msgs:
        #    print ('NashEqFinder took (hh:mm:ss) (processing/wall) time: pyomo\
        #            = {}/{}  ,  solver = {}/{}  ,  run = {}/{} for a game with\
        #            {} cells in its payoff matrix\n'.format(elapsed_optlang_pt,\
        #            elapsed_optlang_wt, elapsed_solver_pt,elapsed_solver_wt, \
        #            elapsed_run_pt,elapsed_run_wt, \
        #            len(self.game.payoff_matrix)) )
            print('NashEqFinder took (hh:mm:ss) (processing/wall) time: pyomo\
                   = {}/{}  ,  solver = {}/{}  ,  run = {}/{} \n'.format(elapsed_optlang_pt,\
                   elapsed_optlang_wt, elapsed_solver_pt,elapsed_solver_wt, \
                   elapsed_run_pt,elapsed_run_wt) )
        
        # return show_matrix(self, self.game.payoff_matri, self.Nash_equilibria, self.strategies, method + "called by optlangFindPure")       
    
    def run(self):
        """
        Runs the Nash equilibrium finder
        """
        if self.NashEq_type.lower() == 'pure':
            self.findPure()
        elif self.NashEq_type.lower() == 'mixed':
            pass # To be completed

        return [self.Nash_equilibria,self.exit_flag]

    
    # Elie
    def optlangRun(self):
        """
        Runs the Nash equilibrium finder
        """
        if self.NashEq_type.lower() == 'pure':
            self.optlangFindPure()
        elif self.NashEq_type.lower() == 'mixed':
            pass # To be completed

        return [self.Nash_equilibria, self.exit_flag, self.game.payoff_matrix]


    # Elie
    def show_matrix(self, payoff_matrix, nash_equilibria, strategies, method): 
        print(f"Time before writing to limited_show_matrix.txt: {datetime.datetime.now()}")
        with open(f"limited_show_matrix.txt", "a") as f:
            f.write(f"\n\nSize = {len(strategies)}\n Method: {method} \nNash Equilibria: {nash_equilibria} \nTime: {datetime.datetime.now()}")

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
                    row.append(
                        (round(payoff_matrix[(('row', i),('column', j))]['row'], 4), \
                        round(payoff_matrix[(('row', i),('column', j))]['column'], 4))
                    )
                table_text.append(row)
                # print("dis row", row)

            return col_text, table_text


        table = payoffs_to_table(payoff_matrix)
        the_table = ax.table(#colWidths=[.3] * (SIZE + 1),
                            cellText=table[1], colLabels=table[0],
                            loc='center', bbox=[-0.05, 0, 1.1, 1.0],
                            rowLoc='center', colLoc='center', cellLoc='center')
        the_table.scale(1, 7)

        def letter_to_position(letter):
            return strategies.index(letter) + 1
            
        for eq in range(len(nash_equilibria)):
            the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].set_facecolor('#5dbcd2')

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

        the_table.auto_set_font_size(False)
        for row in range(SIZE+1):
            for col in range(SIZE+1):
                if (row, col) not in nasheq_positions:
                    the_table[row, col].set_fontsize(3.5)
                else:
                    the_table[row, col].set_fontsize(1.4)

        # for cell in the_table._cells:
            if the_table._cells[cell].xy not in nasheq_positions:
                text = the_table._cells[cell].get_text()
                text.set_fontsize(3.5)

        mpl.savefig(f'PNGs3/size = {len(strategies)}, {method}.png', dpi=1000)
        mpl.show()
        for eq in range(len(nash_equilibria)):
            fontsize = the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].get_fontsize()
            print("Table fontsize", fontsize)
            exit()


    # Elie
    def show_matrix_2c(self, original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells):  
        print(f"Time before writing to limited_show_matrix_2c.txt: {datetime.datetime.now()}")
        # commented it out on 10/20/2022
        # with open(f"limited_show_matrix_2c.txt", "a") as f:
        #     f.write(f"\n\nSize = {len(strategies)}\n Method: {method} \nnash_equilibria: {nash_equilibria} \nchanged_cells: {changed_cells} \nTime: {datetime.datetime.now()}")
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
                        # ("({:s} " + check_0(round(new1 - float(og1), 4)) + ", {:s}" + check_0(round(new2 - float(og2), 4)) + ")").format(str(round(float(og1), 4)), str(round(new1 - float(og1), 4)), str(round(float(og2), 4)), str(round(new2 - float(og2), 4)))
                        f"({remove_zeroes(str(round(float(og1), 4)))}" + check_0(str(round(new1 - float(og1), 4))) + f", {remove_zeroes(str(round(float(og2), 4)))}" + check_0(str(round(new2 - float(og2), 4))) + ")"
                        # ("{:s} + {:s}".format(str(round(og, 4)), str(round(new1 - og, 4))), \
                        #  "{:s} + {:s}".format(str(round(og, 4)), str(round(new2 - og, 4))))
                    )
                table_text.append(row)
                # print("dis row", row)

            return col_text, table_text


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
            return (
                strategies.index(changed_cell[0][0][0][1])+1,
                strategies.index(changed_cell[0][0][1][1])+1
            )
        print(f"Time before setting the cells colors: {datetime.datetime.now()}")
        changed_cells_positions = [changed_cell_to_position(cell) for cell in changed_cells]
        for (row, col) in changed_cells_positions:
            the_table[row, col].set_facecolor('#d2905d')


        def letter_to_position(letter):
            return strategies.index(letter) + 1
        
        print(f"Time before setting the cells colors again: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].set_facecolor('#5dbcd2')

        
        print(f"Time before closed grid edges: {datetime.datetime.now()}")
        for (row, col), cell in the_table.get_celld().items():
            cell.visible_edges = ''
            if col != 0 and row != 0:
                cell.visible_edges = 'closed'

        nasheq_positions = []
        print(f"Time before adding nasheq_positions: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            nasheq_positions.append(
                (letter_to_position(nash_equilibria[eq][0][1]), \
                letter_to_position(nash_equilibria[eq][1][1]))
            )

        # Old before July 6
        # the_table.auto_set_font_size(False)
        # for row in range(SIZE+1):
        #     for col in range(SIZE+1):
        #         # if (row, col) not in nasheq_positions and (row, col) not in changed_cells_positions:
        #         #     the_table[row, col].set_fontsize(3.5)
        #         # else:
        #         #     the_table[row, col].set_fontsize(1.4)


        size = len(strategies)
        the_table.auto_set_font_size(False)
        if size == 2:
            for row in range(size+1):
                for col in range(size+1):
                    if (row, col) not in nasheq_positions:
                        the_table[row, col].set_fontsize(14)
                    else:
                        the_table[row, col].set_fontsize(12)
        elif size == 10:
            for row in range(size+1):
                for col in range(size+1):
                    if (row, col) not in nasheq_positions:
                        the_table[row, col].set_fontsize(4)
                    else:
                        the_table[row, col].set_fontsize(4)

        print(f"Time before saving the PNG: {datetime.datetime.now()}")
        print(f"saving the PNG file for size {len(strategies)}")
        # set title based on iteration
        # if iteration is in method.split(' ') find the word after it
        # if iteration is not in method.split(' ') then just use the method

        if 'iteration' in method.split(' '):
            title = f"Alternative {method.split(' ')[method.split(' ').index('iteration')+1][:-1]}"
        elif 'First' in method:
            title = 'Optimized Game'
        else:
            title = 'Original Game'
        
        # set title with big font
        ax.set_title(f"{title}", fontsize=20)
        mpl.savefig(f'PNGs3/size = {len(strategies)}, {method}.png', dpi=1000)
        mpl.show()


        

    # Elie
    def validate(self, nasheq_cells, method):
        # Validation: needs removal of hard coded methods
        print(f"\n Validation for method {method}\n")

        # Helper function
        def string_to_index(string):
            lst = string.split('_')
            player = lst[-2]
            sign = lst[-1]
            return (((lst[0],lst[1]), (lst[2],lst[3]))), player, sign

        original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        

        # Creating new payoff matrix (i.e. the result of perterbations)
        for var_name, var in self.optModel.variables.items():
            # print(var_name, "=", var.primal)
            matrix_key, player, sign = string_to_index(var_name)
            if sign == 'plus':
                self.game.payoff_matrix[matrix_key][player] += var.primal
            if sign == 'minus':
                self.game.payoff_matrix[matrix_key][player] -= var.primal

        # print('New payoff matrix for', method, self.game.payoff_matrix)

        
        # Define an instance of the NashEqFinder
        NashEqFinderInst = NashEqFinder(self.game, stdout_msgs = True)
        [Nash_equilibria, exit_flag, _game_payoff_matrix] = NashEqFinderInst.optlangRun()


        self.show_matrix_2c(original_payoff_matrix, self.game.payoff_matrix, Nash_equilibria, self.game.players_strategies['row'], method="Original_solution")


        
        # print(f"Validation for {method} returned these results:")
        print ('exit_flag = ', exit_flag)
        print ('Nash_equilibria = ', Nash_equilibria)
        # for desired_state in nasheq_cells:
        #     if list(desired_state) in Nash_equilibria:
        #         # print('DESIRED STATE', desired_state, "ACHIEVED with these perturbations")
        #         for var_name, var in self.optModel.variables.items():
        #             # print(var_name, "=", var.primal)
        #     else:
        #         # print('DESIRED STATE', desired_state, "FAILED")
        #         for var_name, var in self.optModel.variables.items():
        #             # print(var_name, "=", var.primal)

        # print('self.game.players_strategies', self.game.players_strategies)
        self.show_matrix(self.game.payoff_matrix, Nash_equilibria, self.game.players_strategies['row'], method)
        print(f"DONE Validation of {method}")


    
    # # Elie
    # def validate_2c(self, nasheq_cells, method, iteration, time_spent, size):
    #     # Validation: needs removal of hard coded methods
    #     print(f"\n Validation for method {method}\n")

    #     # Helper function
    #     def string_to_index(string):
    #         lst = string.split('_')
    #         player = lst[-2]
    #         sign = lst[-1]
    #         return (((lst[0],lst[1]), (lst[2],lst[3]))), player, sign

    #     original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        

    #     print(f"Time before creating the matrix of the perterbations: {datetime.datetime.now()}")
    #     # Creating new payoff matrix (i.e. the result of perterbations)
    #     for var_name, var_primal in self.current_primals.items():
    #         matrix_key, player, sign = string_to_index(var_name)
    #         if sign == 'plus':
    #             self.game.payoff_matrix[matrix_key][player] += var_primal
    #         if sign == 'minus':
    #             self.game.payoff_matrix[matrix_key][player] -= var_primal

    #     # Get the cells that changed
    #     print(f"Time before getting the changed cells: {datetime.datetime.now()}")
    #     changed_cells = []
    #     for var_name, var_primal in self.current_primals.items():
    #         if var_primal != 0:
    #             changed_cells.append((string_to_index(var_name), var_primal))

    #     # Write changed cells into a file
    #     print(f"Time before writing the changed cells: {datetime.datetime.now()}")
    #     with open(f"results2/size = {size}, {method}_{iteration}_results.txt", "w") as f:
    #         f.write("Changed cells:\n")
    #         for cell in changed_cells:
    #             f.write(str(cell) + "\n")

        
    #     # Define an instance of the NashEqFinder
    #     print(f"Time before defining the NashEqFinder: {datetime.datetime.now()}")
    #     NashEqFinderInst = NashEqFinder(self.game, stdout_msgs = True)
    #     print(f"Time before calling optlangRun: {datetime.datetime.now()}")
    #     [Nash_equilibria, exit_flag, _game_payoff_matrix] = NashEqFinderInst.optlangRun()


    #     # self.show_matrix(original_payoff_matrix, Nash_equilibria, self.game.players_strategies['row'], method="Original_solution")

    #     # add to the file
    #     print(f"Time before writing to results2: {datetime.datetime.now()}")
    #     with open(f"results2/size = {size}, {method}_{iteration}_results.txt", "a") as f:
    #         f.write(f"\n\nNash_equilibria = {Nash_equilibria}\n")
    #         f.write(f"exit_flag = {exit_flag}\n")
    #         f.write(f"time_spent = {time_spent}\n")
  
    #     print(f"Time before calling show_matrix_2c: {datetime.datetime.now()}")
    #     self.show_matrix_2c(original_payoff_matrix, self.game.payoff_matrix, Nash_equilibria, self.game.players_strategies['row'], method, changed_cells)
    #     print(f"DONE Validation of {method}")

    # Elie in 9/5/2022
    def validate_2c(self, nasheq_cells, method, iteration, time_spent, size):
        # Validation: needs removal of hard coded methods
        print(f"\n Validation for method {method}\n")

        # Helper function
        def string_to_index(string):
            lst = string.split('_')
            player = lst[-2]
            sign = lst[-1]
            return (((lst[0],lst[1]), (lst[2],lst[3]))), player, sign

        original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        

        print(f"Time before creating the matrix of the perterbations: {datetime.datetime.now()}")
        # Creating new payoff matrix (i.e. the result of perterbations)
        for var_name, var_primal in self.current_primals.items():
            matrix_key, player, sign = string_to_index(var_name)
            if sign == 'plus':
                self.game.payoff_matrix[matrix_key][player] += var_primal
            if sign == 'minus':
                self.game.payoff_matrix[matrix_key][player] -= var_primal

        # Get the cells that changed
        print(f"Time before getting the changed cells: {datetime.datetime.now()}")
        changed_cells = []
        for var_name, var_primal in self.current_primals.items():
            if var_primal != 0:
                changed_cells.append((string_to_index(var_name), var_primal))

        # Write changed cells into a file
        print(f"Time before writing the changed cells: {datetime.datetime.now()}")
        with open(f"results4/size = {size}, {method}_{iteration}_results.txt", "w") as f:
            f.write("Changed cells:\n")
            for cell in changed_cells:
                f.write(str(cell) + "\n")

        
        # Define an instance of the NashEqFinder
        print(f"Time before checking that the the desired cells are Nash Equilibria: {datetime.datetime.now()}")
        
        # 10/20/2022
        # for cell in nasheq_cells:
        #     # cell is of format [(('row',f'S{i}'), ('column',f'S{j}'))]
        #     # check that the first entry is greater than all first entries that are in the same column
        #     # check that the second entry is greater than all second entries that are in the same row
        #     # if both are true, then the cell is a Nash Equilibrium
        #     # for all payoff_matrix cells that start with ('row',f'S{i}'), check that they are greater than cell[0][1]
        #     row_values = [self.game.payoff_matrix[('row', f"S{i}"), cell[1]]['row'] <= self.game.payoff_matrix[cell]['row'] for i in range(1, size+1) if ('row', f"S{i}") != cell[0]]
        #     column_values = [self.game.payoff_matrix[cell[0], ('column', f"S{i}")]['column'] <= self.game.payoff_matrix[cell]['column'] for i in range(1, size+1) if ('column', f"S{i}") != cell[1]]
        #     if all(row_values) and all(column_values):
        #         print(f"Success: cell {cell} is a Nash Equilibrium")
        #     else:
        #         print(f"Failure: cell {cell} is not a Nash Equilibrium")
                
        # add to the file
        print(f"Time before writing to results4: {datetime.datetime.now()}")
        with open(f"results4/size = {size}, {method}_{iteration}_results.txt", "a") as f:
            f.write(f"\n\nDesired nash_equilibria = {nasheq_cells}\n")

        self.show_matrix_2c(original_payoff_matrix, self.game.payoff_matrix, nasheq_cells, self.game.players_strategies['row'], method, changed_cells)
  
        print(f"DONE Validation of {method}")

    # Elie
    def newEquilibria(self, nasheq_cells, strategies):
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
        

        # Adding new variables 
        model = optlang.Model(name=f'Original Model')
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
        for index in model.indices:
            str_index = str(index).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',', '_')
            var = optlang.Variable(str_index, lb=0, type='continuous', problem=model)
            model.add(var)
            variables_names.append(str_index)

        # Add the objective function
        model.objective = optlang.Objective(expression=sympy.Add(*(sympy.symbols(variables_names))), direction='min')
   
        def strip_down(index):
            return str(index).replace(" ", "").replace('(', "").replace(')', "").replace("'","").replace(',', '_')

        # Add the constraints
        constraints = []
        # Each `cell` is of the format (('row','C'),('column','C'))
        for cell in nasheq_cells:    
            root_index = strip_down(cell)
            # For first player we loop over first axis (rows=strategies)
            for strategy in [x for x in strategies if x != cell[0][1]]:
                current_cell = ((cell[0][0], strategy), (cell[1]))
                current_index = strip_down(current_cell)
                player = model.players[0]
                # This new snippet precludes cells other than `nasheq_cells` 
                # to be a Nash equilibrium
                epsilon = 0.01
                c = optlang.Constraint(
                        self.game.payoff_matrix[cell][player] \
                        + model.variables[root_index+'_'+player+'_plus'] \
                        - model.variables[root_index+'_'+player+'_minus'] \
                        - 
                        (self.game.payoff_matrix[current_cell][player] \
                        + model.variables[current_index+'_'+player+'_plus'] \
                        - model.variables[current_index+'_'+player+'_minus']
                        )
                        -
                        epsilon, lb=0)
                constraints.append(c)

            # For second player we loop over second axis (rows=strategies)
            for strategy in [x for x in strategies if x != cell[1][1]]:
                current_cell = ((cell[0]), (cell[1][0], strategy))
                current_index = strip_down(current_cell)
                player = model.players[1]
                c = optlang.Constraint(
                        self.game.payoff_matrix[cell][player] \
                        + model.variables[root_index+'_'+player+'_plus'] \
                        - model.variables[root_index+'_'+player+'_minus'] \
                        - 
                        (self.game.payoff_matrix[current_cell][player] \
                        + model.variables[current_index+'_'+player+'_plus'] \
                        - model.variables[current_index+'_'+player+'_minus']
                        )
                        -
                        epsilon, lb=0)
                constraints.append(c)
            # print("ADDED columns loop")
        
        # print(constraints)
        model.add(constraints)  

       

        self.optModel = model  
        start_time = datetime.datetime.now()
        self.optModel.optimize()
        end_time = datetime.datetime.now() 

        # Neccessary for validate_2c
        self.current_variables = copy.deepcopy(variables_names)
        self.current_primals = {}
        for var_name, var in self.optModel.variables.items():
            # print(var_name, "=", var.primal)
            self.current_variables.append(var_name)
            self.current_primals[var_name] = var.primal

        original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
        # print("NONONO")
        self.validate_2c(nasheq_cells, method=f'First solution with new equilibria={nasheq_cells} precluded', iteration="NONE", time_spent=end_time-start_time, size=len(strategies))
        self.game.payoff_matrix = original_payoff_matrix

        # New on May 19
        # Iterative Method 2c
        print("\n Using binary variables \n")
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

        
        # Adding new attributes to self
        self.current_variables = copy.deepcopy(variables_names)
        self.current_binary_variables = []
        self.current_binary_constraints = []
        self.current_y_binary_variables = []
        self.current_y_binary_variables_optlang = []
        self.current_primals = {}
        for var_name, var in self.optModel.variables.items():
            # print(var_name, "=", var.primal)
            self.current_variables.append(var_name)
            self.current_primals[var_name] = var.primal

        # adding a binary variable for each α
        for var_name, var_primal in self.current_primals.items():
            a = model.variables[var_name]
            str_index_y_a = var_name + "_binary" + f"_ya" # no iteration needed in this method 
            y_a = optlang.Variable(str_index_y_a, lb=0, type='binary', problem=model)
            model.add(y_a)
            self.current_binary_variables.append(str_index_y_a)
            # Add the binary y variables to the list of binary variables
            self.current_y_binary_variables.append(str_index_y_a)
            self.current_y_binary_variables_optlang.append(y_a)

        # For iteration = 1 to 10
        for iteration in range(5):
            # print the payoff matrix
            # Elie on 1/27/2023
            print(f"\n\n\n----iteration {iteration + 1}----\n\n\n")
            # print(f"----removing binary variables and constrains to start anew----\n")
            self.current_variables = [e for e in self.current_variables if e not in (self.current_binary_variables + self.current_y_binary_variables)]
            print("self.current_variables for iteration", iteration+1, ":", self.current_variables)
            for key_to_remove in self.current_binary_variables:
                if key_to_remove in self.current_primals:
                    del self.current_primals[key_to_remove]
                # del self.current_primals[key_to_remove]
            # On june 20
            epsilon = 1
            ub = 1000
            lb = -1000 
            for var_name, var_primal in self.current_primals.items():
                # Didn't work without it, which is weird
                # Check if the variable optimal is positive and if it is not
                # a binary variable
                if var_primal > 0 and var_name not in (self.current_binary_variables + self.current_y_binary_variables):
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

                    # get the binary variable y_α
                    y_a = model.variables[var_name + "_binary" + f"_ya"]

                    # ϵy_α+(-M)(1-y_α)≤α≤(+M)y_α
                    # ϵy_α+(-M)(1-y_α)≤α
                    c_1 = optlang.Constraint(
                            epsilon * y_a + (-ub) * (1 - y_a) - a,
                            ub=0
                        )
                    # α≤(+M)y_α
                    c_2 = optlang.Constraint(
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
            print("Calculating NZ_a")
            for var_name, var_primal in self.current_primals.items():
                # print(var_name, var_primal)
                if (var_primal > 0) and (var_name not in (self.current_binary_variables + self.current_y_binary_variables)):
                    NZ_a += 1
            print("NZ_a", NZ_a, "\n")
            # print("self.current_y_binary_variables", self.current_y_binary_variables)
            # print("self.current_binary_variables", self.current_binary_variables)
            # print("self.current_variables", self.current_variables)
            # create sypmy expression of NZ_a - 1 - sum(self.current_y_binary_variables)
            # expression = sympy.Add(NZ_a - 1, - 
            #     sympy.Add(*sympy.symbols(self.current_y_binary_variables)))
            # expression = NZ_a - 1 - sympy.Add(*sympy.symbols(self.current_y_binary_variables))

            expression = NZ_a - 1

            for var_name, var_primal in self.current_primals.items():
                # Didn't work without it, which is weird
                # Check if the variable optimal is positive and if it is not
                # a binary variable
                if var_primal > 0 and var_name not in (self.current_binary_variables + self.current_y_binary_variables):
                    # get the binary variable y_α
                    y_a = model.variables[var_name + "_binary" + f"_ya"]
                    expression -= y_a
                # check

            c4 = optlang.Constraint(
                    expression,
                    lb=0
                )
            print("c4", c4)
            print()
            model.add(c4)
            self.current_binary_constraints.append(c4)
            # print("current_binary_constraints")
            # for c in self.current_binary_constraints:
            #     print(c)
            # sympy.Add(*sympy.symbols(self.current_variables))

            # print(sympy.Add(*sympy.symbols(self.current_variables)))
            # print(sympy.symbols(self.current_y_binary_variables))


            # print('binary_variables_names', self.current_binary_variables)
            model.objective = \
                optlang.Objective(
                    expression=sympy.Add(*sympy.symbols(self.current_variables)),
                    direction='min'
                )

            # add every binary variable to current_variables
            for var in self.current_binary_variables:
                self.current_variables.append(var)

            for var in self.current_y_binary_variables:
                self.current_variables.append(var)
            # print(f"Binary variablessss after iteration {iteration+1}:", self.current_binary_variables)

            # print('self.current_binary_constraints:')
            # for con in self.current_binary_constraints:
            #     print(con) 
            
            self.optModel = model  
            # printing model before optimization in iteration + 1
            print(f"Model before optimization in iteration {iteration + 1} is:", self.optModel)
            # print the model before optimization
            start_time = datetime.datetime.now()
            status = self.optModel.optimize()
            print(f"Status of the optimization in iteration {iteration + 1} is:", status)

            end_time = datetime.datetime.now()
            print("self.current_primals for iteration", iteration + 1, "is:", self.current_primals)
            print()
            print("self.optModel.variables for iteration", iteration + 1, "is:")
            for var_name, var in self.optModel.variables.items():
                print(var_name, "=", var)
            for var_name, var in self.optModel.variables.items():
                print(var_name, "=", var.primal)
                self.current_primals[var_name] = var.primal


            changed_binary_variables = []
            for var_name, var_primal in self.current_primals.items():
                if var_name in self.current_binary_variables:
                    if var_primal > 0:
                        changed_binary_variables.append(var_name)

            print("changed_binary_variables for iteration", iteration + 1, "is:", changed_binary_variables)

            
            # with open(f"log2.txt", "a") as f:
            #     f.write(f"\nsize = {len(strategies)}, 2c iteration {iteration} with desired nach equilibria: {nasheq_cells} || Time: {end_time - start_time}")

            
            original_payoff_matrix = copy.deepcopy(self.game.payoff_matrix)
            self.validate_2c(nasheq_cells, method=f"New Method 2c iteration {iteration + 1}, epsilon={epsilon}, with new equilibria={nasheq_cells} precluded", iteration=str(iteration), time_spent=end_time - start_time, size=len(strategies))
            self.game.payoff_matrix = original_payoff_matrix
        


# Elie on 11/18/2022
def show_matrix_original(original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells):  
        print(f"Time before writing to limited_show_matrix_2c.txt: {datetime.datetime.now()}")
        # commented it out on 10/20/2022
        # with open(f"limited_show_matrix_2c.txt", "a") as f:
        #     f.write(f"\n\nSize = {len(strategies)}\n Method: {method} \nnash_equilibria: {nash_equilibria} \nchanged_cells: {changed_cells} \nTime: {datetime.datetime.now()}")
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
                        # ("({:s} " + check_0(round(new1 - float(og1), 4)) + ", {:s}" + check_0(round(new2 - float(og2), 4)) + ")").format(str(round(float(og1), 4)), str(round(new1 - float(og1), 4)), str(round(float(og2), 4)), str(round(new2 - float(og2), 4)))
                        f"({remove_zeroes(str(round(float(og1), 4)))}" + check_0(str(round(new1 - float(og1), 4))) + f", {remove_zeroes(str(round(float(og2), 4)))}" + check_0(str(round(new2 - float(og2), 4))) + ")"
                        # ("{:s} + {:s}".format(str(round(og, 4)), str(round(new1 - og, 4))), \
                        #  "{:s} + {:s}".format(str(round(og, 4)), str(round(new2 - og, 4))))
                    )
                table_text.append(row)
                # print("dis row", row)

            return col_text, table_text


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
            return (
                strategies.index(changed_cell[0][0][0][1])+1,
                strategies.index(changed_cell[0][0][1][1])+1
            )
        print(f"Time before setting the cells colors: {datetime.datetime.now()}")
        changed_cells_positions = [changed_cell_to_position(cell) for cell in changed_cells]
        for (row, col) in changed_cells_positions:
            the_table[row, col].set_facecolor('#d2905d')


        def letter_to_position(letter):
            return strategies.index(letter) + 1
        
        print(f"Time before setting the cells colors again: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].set_facecolor('#5dbcd2')

        
        print(f"Time before closed grid edges: {datetime.datetime.now()}")
        for (row, col), cell in the_table.get_celld().items():
            cell.visible_edges = ''
            if col != 0 and row != 0:
                cell.visible_edges = 'closed'

        nasheq_positions = []
        print(f"Time before adding nasheq_positions: {datetime.datetime.now()}")
        for eq in range(len(nash_equilibria)):
            nasheq_positions.append(
                (letter_to_position(nash_equilibria[eq][0][1]), \
                letter_to_position(nash_equilibria[eq][1][1]))
            )

        # Old before July 6
        # the_table.auto_set_font_size(False)
        # for row in range(SIZE+1):
        #     for col in range(SIZE+1):
        #         # if (row, col) not in nasheq_positions and (row, col) not in changed_cells_positions:
        #         #     the_table[row, col].set_fontsize(3.5)
        #         # else:
        #         #     the_table[row, col].set_fontsize(1.4)


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

        print(f"Time before saving the PNG: {datetime.datetime.now()}")
        print(f"saving the PNG file for size {len(strategies)}")
        # set title based on iteration
        # if iteration is in method.split(' ') find the word after it
        # if iteration is not in method.split(' ') then just use the method

        if 'iteration' in method.split(' '):
            title = f"Iteration {method.split(' ')[method.split(' ').index('iteration')+1][:-1]}"
        else:
            title = 'Original Game'
        
        # set title with big font
        ax.set_title(f"{title}", fontsize=20)
        mpl.savefig(f'PNGs3/size = {len(strategies)}, {method}.png', dpi=1000)
        mpl.show()



# Elie
def show_matrix(payoff_matrix, nash_equilibria, strategies, method):  
    # # Open and append to the file
    # with open(f"log.txt", "a") as f:
    #     f.write(f"\n\nMethod: {method} || Nash Equilibria: {nash_equilibria} || Time: {datetime.datetime.now()}")
    print(f"Time before writing to limited_show_matrix.txt: {datetime.datetime.now()}")
    with open(f"limited_show_matrix.txt", "a") as f:
            f.write(f"\n\nSize = {len(strategies)}\n Method: {method} \nNash Equilibria: {nash_equilibria} \nTime: {datetime.datetime.now()}")
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
                row.append(
                    (round(payoff_matrix[(('row', i),('column', j))]['row'], 4), \
                     round(payoff_matrix[(('row', i),('column', j))]['column'], 4))
                )
            table_text.append(row)

        return col_text, table_text


    table = payoffs_to_table(payoff_matrix)
    the_table = ax.table(#colWidths=[0.3] * (SIZE + 1),
                        cellText=table[1], colLabels=table[0],
                        loc='center', bbox=[-0.05, 0, 1.1, 1.0],
                        rowLoc='center', colLoc='center', cellLoc='center')
    the_table.scale(1, 7)

    def letter_to_position(letter):
        return strategies.index(letter) + 1
        
    for eq in range(len(nash_equilibria)):
        the_table[
            letter_to_position(nash_equilibria[eq][0][1]), 
            letter_to_position(nash_equilibria[eq][1][1])
            ].set_facecolor('#5dbcd2')

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


    the_table.auto_set_font_size(False)
    # the_table.set_fontsize(3.5)
    for row in range(SIZE+1):
        for col in range(SIZE+1):
            if (row, col) not in nasheq_positions:
                the_table[row, col].set_fontsize(3.5)
            else:
                the_table[row, col].set_fontsize(1.4)

    # for (row, col), cell in the_table.get_celld().items():
    #     if (row, col) not in nasheq_positions:
    for cell in the_table._cells:
        # print("cell", cell)
        # print("the_table._cells[cell]", the_table._cells[cell])
        # if the_table._cells[cell].facecolor != '#5dbcd2':
        if the_table._cells[cell].xy not in nasheq_positions:
            text = the_table._cells[cell].get_text()
            text.set_fontsize(3.5)

    mpl.savefig(f'PNGs3/size = {len(strategies)}, {method}.png', dpi = 1000)
    mpl.show()
    for eq in range(len(nash_equilibria)):
            fontsize = the_table[
                letter_to_position(nash_equilibria[eq][0][1]), 
                letter_to_position(nash_equilibria[eq][1][1])
                ].get_fontsize()
            print("Table fontsize", fontsize)
            exit()

#--------- Sample implementation ------
if __name__ == "__main__":
    print("(NashEqFinder")

    from game import *
    print ("\n\n\n\n\n\n\n\n\n\n")
    
    #---------------------------------- 
    # print ("\n-- Prisoner's Dilemma ---")
    # Pure strategy Nash eq = (D,D)
    
    game_name = "Prisoner's Dilemma"
    numberOfPlayers = 2
    players_names = ['row','column']
    
    players_strategies = {}
    players_strategies['row'] = ['C','D']
    players_strategies['column'] = ['C','D']
    
    payoff_matrix = {}
    payoff_matrix[(('row','C'),('column','C'))] = {'row':-1,'column':-1}
    payoff_matrix[(('row','C'),('column','D'))] = {'row':-4,'column':0}
    payoff_matrix[(('row','D'),('column','C'))] = {'row':0,'column':-4}
    payoff_matrix[(('row','D'),('column','D'))] = {'row':-3,'column':-3}
    
    # Define an instance of the game
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # (self, original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells)
    show_matrix_original(game_payoff_matrix, game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Payoff Matrix", [])       

    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','C'), ('column','C'))], strategies=['C', 'D'])
    
    
    # ----------------------------------
    # print ("\n-- Game of size 10 ---")
    # SIZE = 10
    # game_name = f"{SIZE} strategies"
    # players_names = ['row','column']
    # players_strategies = {}
    # strategies = ['S' + str(i) for i in range(1,SIZE+1)]
    # players_strategies['row'] = strategies
    # players_strategies['column'] = strategies
    # payoff_matrix = {}
    # def read_matrix():
    #     print("Reading the payoff matrix from the file")
    #     with open(f"random/payoff_matrix_{SIZE}.txt", "r") as f:
    #         lines = f.readlines()
    #         print(len(lines))
    #         print(len(lines[0].split()))
    #         for i in range(SIZE):
    #             for j in range(SIZE):
    #                 payoff_matrix[('row', f"S{i+1}"), ('column', f"S{j+1}")] = \
    #                     {'row': int(lines[i].split()[2*j]), 'column': int(lines[i].split()[2*j+1])}

    # # write_matrix()
    # read_matrix()
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # [Nash_equilibria_, exit_flag_] = NashEqFinderInst.run()
    # show_matrix(game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Game called by optlangFindPure")
    # show_matrix(payoff_matrix, Nash_equilibria_, players_strategies['row'], "Original Game called by FindPure")
    # print("Nash_equilibria optlangFindPure = ", Nash_equilibria)
    # print("Nash_equilibria FindPure        = ", Nash_equilibria_)
    # # order the equilibria
    # Nash_equilibria.sort()
    # Nash_equilibria_.sort()
    # print(Nash_equilibria_==Nash_equilibria)
    # # exit()

    # NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','S5'), ('column','S3'))], strategies=strategies)
    
    exit()



    # # Validation:
    # payoff_matrix = {}
    # payoff_matrix[(('row','C'),('column','C'))] = {'row':-1,'column':-1}
    # payoff_matrix[(('row','C'),('column','D'))] = {'row':-4 + 1,'column':0} #validation
    # payoff_matrix[(('row','D'),('column','C'))] = {'row':0,'column':-4}
    # payoff_matrix[(('row','D'),('column','D'))] = {'row':-3,'column':-3}
    
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.optlangRun()

    # print("DONE")
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )


    # ----------------------------------
    
    # create a SIZE strategies names S1, S2, S3, ..., S100
    import random
    game_name = f"{SIZE} strategies"
    players_names = ['row','column']
    
    players_strategies = {}
    strategies = ['S' + str(i) for i in range(1,SIZE+1)]
    players_strategies['row'] = strategies
    players_strategies['column'] = strategies
    payoff_matrix = {}
    # pupulate the payoff matrix with random payoffs upper triangle
    for i in range(1,SIZE+1):
        for j in range(i,SIZE+1):
            payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] = \
                {'row': random.randint(-15,0), 'column': random.randint(-15,0)}
        
    # copy into the lower triangle
    for i in range(1,SIZE+1):
        for j in range(i,SIZE+1):
            payoff_matrix[('row', f"S{j}"), ('column', f"S{i}")] = \
                payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] 

    # write the payoff matrix to a file
    def write_matrix():
        print("Writing matrix to file")
        with open(f"random/payoff_matrix_{SIZE}.txt", "w") as f:
            for i in range(1,SIZE+1):
                for j in range(1,SIZE+1):
                    f.write(f"""{payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")]['row']} """)
                    f.write(f"""{payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")]['column']} """)
                f.write("\n")

    # read the payoff matrix from the file
    def read_matrix():
        print("Reading the payoff matrix from the file")
        with open(f"random/payoff_matrix_{SIZE}.txt", "r") as f:
            lines = f.readlines()
            print(len(lines))
            print(len(lines[0].split()))
            for i in range(SIZE):
                for j in range(SIZE):
                    payoff_matrix[('row', f"S{i+1}"), ('column', f"S{j+1}")] = \
                        {'row': int(lines[i].split()[2*j]), 'column': int(lines[i].split()[2*j+1])}

    # write_matrix()
    read_matrix()



    # # read the payoff matrix from a file
    # with open(f"payoff_matrix_{SIZE}.txt", "r") as f:
    #     lines = f.readlines()
    #     for i in range(1,SIZE+1):
    #         for j in range(1,SIZE+1):
    #             payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] = \
    #                 {'row': int(lines[i-1][j-1]), 'column': int(lines[i-1][j-1])}
    
    # Define an instance of the game
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    [Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    [Nash_equilibria_, exit_flag_] = NashEqFinderInst.run()
    show_matrix(game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Game called by optlangFindPure")
    show_matrix(payoff_matrix, Nash_equilibria_, players_strategies['row'], "Original Game called by FindPure")
    print("Nash_equilibria optlangFindPure = ", Nash_equilibria)
    print("Nash_equilibria FindPure        = ", Nash_equilibria_)
    # order the equilibria
    Nash_equilibria.sort()
    Nash_equilibria_.sort()
    print(Nash_equilibria_==Nash_equilibria)
    # exit()

    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','S15'), ('column','S10'))], strategies=strategies)
    # NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','S5'), ('column','S3'))], strategies=strategies)



    
    # #---------------------------------- 
    # print ("\n-- Game of pure coordination ---")
    # # Pure strategy Nash eq: (Left,Left) and (Right,Right)
    
    # game_name = "Pure coordination"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['Left','Right']
    # players_strategies['column'] = ['Left','Right']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','Left'),('column','Left'))] = {'row':1,'column':1}
    # payoff_matrix[(('row','Left'),('column','Right'))] = {'row':0,'column':0}
    # payoff_matrix[(('row','Right'),('column','Left'))] = {'row':0,'column':0}
    # payoff_matrix[(('row','Right'),('column','Right'))] = {'row':1,'column':1}
    
    # # Define an instance of the game
    # PC = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PC, stdout_msgs = True)
    # # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.optlangRun()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )
    
    # #---------------------------------- 
    # print ("\n-- Game of Battle of the sexes ---")
    # # Pure strategy Nash eq: (B,B) and (F,F)
    
    # game_name = "Battle of the sexes"
    # numberOfPlayers = 2
    # players_names = ['husband','wife']
    
    # players_strategies = {}
    # players_strategies['husband'] = ['B','F']
    # players_strategies['wife'] = ['B','F']
    
    # payoff_matrix = {}
    # payoff_matrix[(('husband','B'),('wife','B'))] = {'husband':2,'wife':1}
    # payoff_matrix[(('husband','B'),('wife','F'))] = {'husband':0,'wife':0}
    # payoff_matrix[(('husband','F'),('wife','B'))] = {'husband':0,'wife':0}
    # payoff_matrix[(('husband','F'),('wife','F'))] = {'husband':1,'wife':2}
    
    # # Define an instance of the game
    # BS = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(BS, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )
    
    # #---------------------------------- 
    # print ("\n-- Game of Matching pennies ---")
    # # Pure strategy Nash eq: None
    
    # game_name = "Matching pennies"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['Heads','Tails']
    # players_strategies['column'] = ['Heads','Tails']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','Heads'),('column','Heads'))] = {'row':1,'column':-1}
    # payoff_matrix[(('row','Heads'),('column','Tails'))] = {'row':-1,'column':1}
    # payoff_matrix[(('row','Tails'),('column','Heads'))] = {'row':-1,'column':1}
    # payoff_matrix[(('row','Tails'),('column','Tails'))] = {'row':1,'column':-1}
    
    # # Define an instance of the game
    # MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    # [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )
    
    # #---------------------------------- 
    # print ("\n-- Problem 4- Homework 1 (game theory I) ---")
    # # This is a game with two players and multiple strategies
    # # Pure strategy Nash eq: (c,y)
    # game_name = "Hw1Prob4"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['a','b','c','d']
    # players_strategies['column'] = ['x','y','z']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','a'),('column','x'))] = {'row':1,'column':2}
    # payoff_matrix[(('row','a'),('column','y'))] = {'row':2,'column':2}
    # payoff_matrix[(('row','a'),('column','z'))] = {'row':5,'column':1}
    
    # payoff_matrix[(('row','b'),('column','x'))] = {'row':4,'column':1}
    # payoff_matrix[(('row','b'),('column','y'))] = {'row':3,'column':5}
    # payoff_matrix[(('row','b'),('column','z'))] = {'row':3,'column':3}
    
    # payoff_matrix[(('row','c'),('column','x'))] = {'row':5,'column':2}
    # payoff_matrix[(('row','c'),('column','y'))] = {'row':4,'column':4}
    # payoff_matrix[(('row','c'),('column','z'))] = {'row':7,'column':0}
    
    # payoff_matrix[(('row','d'),('column','x'))] = {'row':2,'column':3}
    # payoff_matrix[(('row','d'),('column','y'))] = {'row':0,'column':4}
    # payoff_matrix[(('row','d'),('column','z'))] = {'row':3,'column':0}
    
    # # Define an instance of the game
    # Hw1Pb4 = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(Hw1Pb4, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )
    
    
    # #---------------------------------- 
    # print ("\n-- Problem 9- Homework 1 (game theory I) ---")
    # # This is a game with three players and two strategies
    # # Pure strategy Nash eq: (c,y)
    
    # game_name = "Pure coordination"
    # numberOfPlayers = 3
    # players_names = ['voter1','voter2','voter3']
    
    # players_strategies = {}
    # players_strategies['voter1'] = ['candidateA','candidateB']
    # players_strategies['voter2'] = ['candidateA','candidateB']
    # players_strategies['voter3'] = ['candidateA','candidateB']
    
    # payoff_matrix = {}
    # payoff_matrix[(('voter1','candidateA'),('voter2','candidateA'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    # payoff_matrix[(('voter1','candidateA'),('voter2','candidateA'),('voter3','candidateB'))] = {'voter1':1,'voter2':0,'voter3':0}
    # payoff_matrix[(('voter1','candidateA'),('voter2','candidateB'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    # payoff_matrix[(('voter1','candidateA'),('voter2','candidateB'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    # payoff_matrix[(('voter1','candidateB'),('voter2','candidateA'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    # payoff_matrix[(('voter1','candidateB'),('voter2','candidateA'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    # payoff_matrix[(('voter1','candidateB'),('voter2','candidateB'),('voter3','candidateA'))] = {'voter1':0,'voter2':1,'voter3':1}
    # payoff_matrix[(('voter1','candidateB'),('voter2','candidateB'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    
    # # Define an instance of the game
    # Hw1Pb9 = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(Hw1Pb9, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    # #---------------------------------- 
    # print ("\n-- Mutualism 1 ---")
    # # Pure strategy Nash eq: None
    
    # game_name = "Mutualism"
    # numberOfPlayers = 2
    # players_names = ['m1','m2']
    
    # players_strategies = {}
    # players_strategies['m1'] = ['C','D']
    # players_strategies['m2'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('m1','C'),('m2','C'))] = {'m1':5,'m2':6}
    # payoff_matrix[(('m1','C'),('m2','D'))] = {'m1':-1,'m2':8}
    # payoff_matrix[(('m1','D'),('m2','C'))] = {'m1':7,'m2':-2}
    # payoff_matrix[(('m1','D'),('m2','D'))] = {'m1':0,'m2':0}
    
    # # Define an instance of the game
    # MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    # [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    # #---------------------------------- 
    # print ("\n-- Mutualism 2 ---")
    # # Pure strategy Nash eq: None
    
    # game_name = "Mutualism"
    # numberOfPlayers = 2
    # players_names = ['p1','p2']
    
    # players_strategies = {}
    # players_strategies['p1'] = ['C1','C2','D1','D2']
    # players_strategies['p2'] = ['C1','C2','D1','D2']
    
    # payoff_matrix = {}
    # payoff_matrix[(('p1','C1'),('p2','C1'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','C1'),('p2','D1'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','C1'),('p2','C2'))] = {'p1':5,'p2':3}
    # payoff_matrix[(('p1','C1'),('p2','D2'))] = {'p1':-1,'p2':8}

    # payoff_matrix[(('p1','D1'),('p2','C1'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','D1'),('p2','D1'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','D1'),('p2','C2'))] = {'p1':7,'p2':-2}
    # payoff_matrix[(('p1','D1'),('p2','D2'))] = {'p1':0,'p2':0}

    # payoff_matrix[(('p1','C2'),('p2','C1'))] = {'p1':3,'p2':5}
    # payoff_matrix[(('p1','C2'),('p2','D1'))] = {'p1':-2,'p2':7}
    # payoff_matrix[(('p1','C2'),('p2','C2'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','C2'),('p2','D2'))] = {'p1':0,'p2':0}

    # payoff_matrix[(('p1','D2'),('p2','C1'))] = {'p1':8,'p2':-1}
    # payoff_matrix[(('p1','D2'),('p2','D1'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','D2'),('p2','C2'))] = {'p1':0,'p2':0}
    # payoff_matrix[(('p1','D2'),('p2','D2'))] = {'p1':0,'p2':0}
    
    # # Define an instance of the game
    # MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    # [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )


    # #---------------------------------- 
    # print ("\n-- Synergism ---")
    # # Pure strategy Nash eq: None
    
    # game_name = "Synergism"
    # numberOfPlayers = 2
    # players_names = ['m1','m2']
    
    # players_strategies = {}
    # players_strategies['m1'] = ['C','D']
    # players_strategies['m2'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('m1','C'),('m2','C'))] = {'m1':5,'m2':6}
    # payoff_matrix[(('m1','C'),('m2','D'))] = {'m1':1,'m2':6}
    # payoff_matrix[(('m1','D'),('m2','C'))] = {'m1':5,'m2':2}
    # payoff_matrix[(('m1','D'),('m2','D'))] = {'m1':1,'m2':2}
    
    # # Define an instance of the game
    # MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    # [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    # #---------------------------------- 
    # print ("\n-- Commensalism ---")
    # # Pure strategy Nash eq: None
    
    # game_name = "Commensalism"
    # numberOfPlayers = 2
    # players_names = ['m1','m2']
    
    # players_strategies = {}
    # players_strategies['m1'] = ['C','D']
    # players_strategies['m2'] = ['D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('m1','C'),('m2','D'))] = {'m1':5,'m2':6}
    # payoff_matrix[(('m1','D'),('m2','D'))] = {'m1':5,'m2':2}
    
    # # Define an instance of the game
    # MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    # [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    
    # #---------------------------------- 
    # print ("\n-- Maya's game (e = 0.8 sATP = 5) --> Mutually Beneficial ---")
    
    # game_name = "Maya's game: e = 0.8 sATP = 5"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['C','D']
    # players_strategies['column'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','C'),('column','C'))] = {'row':0.044,'column':0.044}
    # payoff_matrix[(('row','C'),('column','D'))] = {'row':0.039,'column':0.008}
    # payoff_matrix[(('row','D'),('column','C'))] = {'row':0.008,'column':0.039}
    # payoff_matrix[(('row','D'),('column','D'))] = {'row':-0.0016,'column':-0.0016}
    
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )
    
    # #---------------------------------- 
    # print ("\n-- Maya's game (e = 0.01 sATP = 5) --> Prisoner's Dilemma ---")
    
    # game_name = "Maya's game: e = 0.01 sATP = 5"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['C','D']
    # players_strategies['column'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','C'),('column','C'))] = {'row':0.28,'column':0.28}
    # payoff_matrix[(('row','C'),('column','D'))] = {'row':0.0011,'column':0.27}
    # payoff_matrix[(('row','D'),('column','C'))] = {'row':0.27,'column':0.0011}
    # payoff_matrix[(('row','D'),('column','D'))] = {'row':-0.0016,'column':-0.0016}
    
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    # #---------------------------------- 
    # print ("\n-- Maya's game (e = 0.4 sATP = 5) - Snowdirft ---")
    
    # game_name = "Maya's game: e = 0.8 sATP = 5"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['C','D']
    # players_strategies['column'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','C'),('column','C'))] = {'row':0.063,'column':0.063}
    # payoff_matrix[(('row','C'),('column','D'))] = {'row':0.022,'column':0.060}
    # payoff_matrix[(('row','D'),('column','C'))] = {'row':0.060,'column':0.022}
    # payoff_matrix[(('row','D'),('column','D'))] = {'row':-0.0016,'column':-0.0016}
    
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )

    # #---------------------------------- 
    # print ("\n-- Elie's game (e = 0.01 sATP = 5) --> Prisoner's Dilemma ---")
    
    # game_name = "Elie's game: e = 0.01 sATP = 5"
    # numberOfPlayers = 2
    # players_names = ['row','column']
    
    # players_strategies = {}
    # players_strategies['row'] = ['C','D']
    # players_strategies['column'] = ['C','D']
    
    # payoff_matrix = {}
    # payoff_matrix[(('row','C'),('column','C'))] = {'row':0.28,'column':0.28}
    # payoff_matrix[(('row','C'),('column','D'))] = {'row':0.0011,'column':0.27}
    # payoff_matrix[(('row','D'),('column','C'))] = {'row':0.27,'column':0.0011}
    # payoff_matrix[(('row','D'),('column','D'))] = {'row':-0.0016,'column':-0.0016}
    
    # # Define an instance of the game
    # PD = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # # Define an instance of the NashEqFinder
    # NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    # print ('exit_flag = ',exit_flag)
    # print ('Nash_equilibria = ',Nash_equilibria )


