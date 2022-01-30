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
import re, sys, math, copy, time, random
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from coopr.pyomo import *
from coopr.opt import *
from pyomoSolverCreator import *

class NashEqFinder(object):
    """
    General class for NashEq Finder. Sample usage is provided at the end 
    """   

    def __init__(self,game, NashEq_type = 'pure', optimization_solver = 'gurobi', warnings = True, stdout_msgs = True, output_file = ''):
        """
        INPUTS 
        ------
        game: 
        An instance of the class game (see game.py for details) 

        NashEq_type: 
        Type of the Nash equilibrium to find (currently only pure strategy Nash equilibria)

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
        file and its name (e.g., 'results/fbaResults.txt'), where
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
        start_run_pt = time.clock()
        start_run_wt = time.time()

        #---- Creating and instantiating the optModel ----
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()

        # Create the optModel model        
        self.createPyomoModel()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optimization_solver)

        elapsed_pyomo_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
        elapsed_pyomo_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))

        #-- Some initializations --
        # Instantiate the optModel with new fixed variables
        self.optModel.preprocess()

        #- Solve the optModel (tee=True shows the solver output) -
        try:
            start_solver_pt = time.clock()
            start_solver_wt = time.time()

            optSoln = solverType.solve(self.optModel,tee=False)
            solverFlag = 'normal'
    
        # In the case of an error switch the solver
        except:
            if self.warnings:
                print "WARNING! ",self.optimization_solver," failed. An alternative solver is tried"        
    
            if self.optimization_solver.lower() == 'gurobi':
                self.optimization_solver = 'cplex'
            elif self.optimization_solver.lower() == 'cplex':
                self.optimization_solver = 'gurobi'
    
            # Try solving with the alternative solver
            solverType = pyomoSolverCreator(self.optimization_solver)
            try:
                start_solver_pt = time.clock()
                start_solver_wt = time.time()

                optSoln = solverType.solve(self.optModel,tee=False)
                solverFlag = 'normal'
            except:
                solverFlag = 'solverError'
                if self.warnings:
                    print '\nWARNING! The alternative solver failed. No solution was returned'

        elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))
    
        #----- Print the results in the output (screen, file and/or variable) ------
        # Load the results
        self.optModel.load(optSoln)
            
        # Set of the Nash equilibria
        self.Nash_equilibria = []
        
        if solverFlag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':
            
            optimExitflag = 'globallyOptimal'
    
            # Value of the objective function
            objValue = self.optModel.objective_rule()
    
            # Print the results on the screen 
            if self.stdout_msgs:
                print "\nsolver.status = ",optSoln.solver.termination_condition,"\n"
                print "objective value = ",objValue

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
                print "\nWARNING! No optimal solutions found (solution.solver.status = ",optSoln.Solution.status,", solver.status =",optSoln.solver.status,", solver.termination_condition = ",optSoln.solver.termination_condition,")\n"
    
            # Write the results into the output file
            if self.output_file != None: 
                pass    # *** To be completed ***
            else:
                pass
    
        # Time required to run 
        elapsed_run_pt = str(timedelta(seconds = time.clock() - start_run_pt))
        elapsed_run_wt = str(timedelta(seconds = time.time() - start_run_wt))
    
        if self.stdout_msgs:
           print 'NashEqFinder took (hh:mm:ss) (processing/wall) time: pyomo = {}/{}  ,  solver = {}/{}  ,  run = {}/{} for a game with {} cells in its payoff matrix\n'.format(elapsed_pyomo_pt,elapsed_pyomo_wt,elapsed_solver_pt,elapsed_solver_wt,elapsed_run_pt,elapsed_run_wt, len(self.game.payoff_matrix)) 
    
    def run(self):
        """
        Runs the Nash equilibrium finder
        """
        if self.NashEq_type.lower() == 'pure':
            self.findPure()
        else:
            pass 

        return [self.Nash_equilibria,self.exit_flag]

#--------- Sample implementation ------
if __name__ == "__main__":

    from game import *
    
    print "\n-- Prisoner's Dilemma ---"
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
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 
    
    
    print "\n-- Game of pure coordination ---"
    # Pure strategy Nash eq: (Left,Left) and (Right,Right)
    
    game_name = "Pure coordination"
    numberOfPlayers = 2
    players_names = ['row','column']
    
    players_strategies = {}
    players_strategies['row'] = ['Left','Right']
    players_strategies['column'] = ['Left','Right']
    
    payoff_matrix = {}
    payoff_matrix[(('row','Left'),('column','Left'))] = {'row':1,'column':1}
    payoff_matrix[(('row','Left'),('column','Right'))] = {'row':0,'column':0}
    payoff_matrix[(('row','Right'),('column','Left'))] = {'row':0,'column':0}
    payoff_matrix[(('row','Right'),('column','Right'))] = {'row':1,'column':1}
    
    # Define an instance of the game
    PC = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(PC, stdout_msgs = True)
    [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 
    
    print "\n-- Game of Battle of the sexes ---"
    # Pure strategy Nash eq: (B,B) and (F,F)
    
    game_name = "Battle of the sexes"
    numberOfPlayers = 2
    players_names = ['husband','wife']
    
    players_strategies = {}
    players_strategies['husband'] = ['B','F']
    players_strategies['wife'] = ['B','F']
    
    payoff_matrix = {}
    payoff_matrix[(('husband','B'),('wife','B'))] = {'husband':2,'wife':1}
    payoff_matrix[(('husband','B'),('wife','F'))] = {'husband':0,'wife':0}
    payoff_matrix[(('husband','F'),('wife','B'))] = {'husband':0,'wife':0}
    payoff_matrix[(('husband','F'),('wife','F'))] = {'husband':1,'wife':2}
    
    # Define an instance of the game
    BS = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(BS, stdout_msgs = True)
    [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 
    
    print "\n-- Game of Matching pennies ---"
    # Pure strategy Nash eq: None
    
    game_name = "Matching pennies"
    numberOfPlayers = 2
    players_names = ['row','column']
    
    players_strategies = {}
    players_strategies['row'] = ['Heads','Tails']
    players_strategies['column'] = ['Heads','Tails']
    
    payoff_matrix = {}
    payoff_matrix[(('row','Heads'),('column','Heads'))] = {'row':1,'column':-1}
    payoff_matrix[(('row','Heads'),('column','Tails'))] = {'row':-1,'column':1}
    payoff_matrix[(('row','Tails'),('column','Heads'))] = {'row':-1,'column':1}
    payoff_matrix[(('row','Tails'),('column','Tails'))] = {'row':1,'column':-1}
    
    # Define an instance of the game
    MP = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(MP, stdout_msgs = True)
    [Nash_equilibria,exit_flag]  = NashEqFinderInst.run()
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 
    
    print "\n-- Test game with two players and several strategies for each player ---"
    # This is a game with two players and multiple strategies
    # Pure strategy Nash eq: (c,y)
    game_name = "Hw1Prob4"
    numberOfPlayers = 2
    players_names = ['row','column']
    
    players_strategies = {}
    players_strategies['row'] = ['a','b','c','d']
    players_strategies['column'] = ['x','y','z']
    
    payoff_matrix = {}
    payoff_matrix[(('row','a'),('column','x'))] = {'row':1,'column':2}
    payoff_matrix[(('row','a'),('column','y'))] = {'row':2,'column':2}
    payoff_matrix[(('row','a'),('column','z'))] = {'row':5,'column':1}
    
    payoff_matrix[(('row','b'),('column','x'))] = {'row':4,'column':1}
    payoff_matrix[(('row','b'),('column','y'))] = {'row':3,'column':5}
    payoff_matrix[(('row','b'),('column','z'))] = {'row':3,'column':3}
    
    payoff_matrix[(('row','c'),('column','x'))] = {'row':5,'column':2}
    payoff_matrix[(('row','c'),('column','y'))] = {'row':4,'column':4}
    payoff_matrix[(('row','c'),('column','z'))] = {'row':7,'column':0}
    
    payoff_matrix[(('row','d'),('column','x'))] = {'row':2,'column':3}
    payoff_matrix[(('row','d'),('column','y'))] = {'row':0,'column':4}
    payoff_matrix[(('row','d'),('column','z'))] = {'row':3,'column':0}
    
    # Define an instance of the game
    Hw1Pb4 = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(Hw1Pb4, stdout_msgs = True)
    [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 
    
    
    print "\n-- Test game with three players ---"
    # This is a game with three players and two strategies
    # Pure strategy Nash eq: (c,y)
    
    game_name = "Pure coordination"
    numberOfPlayers = 3
    players_names = ['voter1','voter2','voter3']
    
    players_strategies = {}
    players_strategies['voter1'] = ['candidateA','candidateB']
    players_strategies['voter2'] = ['candidateA','candidateB']
    players_strategies['voter3'] = ['candidateA','candidateB']
    
    payoff_matrix = {}
    payoff_matrix[(('voter1','candidateA'),('voter2','candidateA'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    payoff_matrix[(('voter1','candidateA'),('voter2','candidateA'),('voter3','candidateB'))] = {'voter1':1,'voter2':0,'voter3':0}
    payoff_matrix[(('voter1','candidateA'),('voter2','candidateB'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    payoff_matrix[(('voter1','candidateA'),('voter2','candidateB'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    payoff_matrix[(('voter1','candidateB'),('voter2','candidateA'),('voter3','candidateA'))] = {'voter1':1,'voter2':0,'voter3':0}
    payoff_matrix[(('voter1','candidateB'),('voter2','candidateA'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    payoff_matrix[(('voter1','candidateB'),('voter2','candidateB'),('voter3','candidateA'))] = {'voter1':0,'voter2':1,'voter3':1}
    payoff_matrix[(('voter1','candidateB'),('voter2','candidateB'),('voter3','candidateB'))] = {'voter1':0,'voter2':1,'voter3':1}
    
    # Define an instance of the game
    Hw1Pb9 = game(game_name, players_names, players_strategies, payoff_matrix)
    
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(Hw1Pb9, stdout_msgs = True)
    [Nash_equilibria,exit_flag] = NashEqFinderInst.run()
    
    print 'exit_flag = ',exit_flag
    print 'Nash_equilibria = ',Nash_equilibria 

    
