# create a 100 strategies names S1, S2, S3, ..., S100
import random
import sys
sys.path.append('/Users/elieeshoa/Dropbox/Elie_Eshoa/Ali_codes/')
from pyomoSolverCreator import *
from pyomo.environ import *
from pyomo.opt import *
from game import *
from NashEqFinder import *


# strategies = ['S' + str(i) for i in range(1,101)]

# payoff_matrix = {}

# # pupulate the payoff matrix with random payoffs
# for i in range(1,101):
#     for j in range(1,101):
#         payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] = \
#             {'row': random.randint(-15,0), 'column': random.randint(-15,0)}


def generate_random_matrix(SIZE):
    payoff_matrix = {}  
    for i in range(1,SIZE+1):
        for j in range(1,SIZE+1):
            payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] = \
                {'row': random.randint(-15,0), 'column': random.randint(-15,0)}
    return payoff_matrix


    
# write the payoff matrix to a file
def write_matrix(payoff_matrix, SIZE):
    print("Writing matrix to file a matrix of size ", SIZE)
    with open(f"matrices/payoff_matrix_{SIZE}.txt", "w") as f:
        for i in range(1,SIZE+1):
            for j in range(1,SIZE+1):
                f.write(f"""{payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")]['row']} """)
                f.write(f"""{payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")]['column']} """)
            f.write("\n")

# read the payoff matrix from the file
def read_matrix(payoff_matrix, SIZE):
    print("Reading the payoff matrix from the file")
    with open(f"matrices/payoff_matrix_{SIZE}.txt", "r") as f:
        lines = f.readlines()
        # print(len(lines))
        # print(len(lines[0].split()))
        for i in range(SIZE):
            for j in range(SIZE):
                payoff_matrix[('row', f"S{i+1}"), ('column', f"S{j+1}")] = \
                    {'row': int(lines[i].split()[2*j]), 'column': int(lines[i].split()[2*j+1])}

def gen(sizes):
    # pupulate the payoff matrix with random payoffs upper triangle
    for SIZE in sizes:
        payoff_matrix = generate_random_matrix(SIZE)
        write_matrix(payoff_matrix, SIZE)







# create a payoff matrix with only one nash equilibrium
def generate_nash_equilibrium(SIZE):
    payoff_matrix = generate_random_matrix(SIZE)
    # set one random cell to (1, 1)
    nash_i = random.randint(1, SIZE)
    nash_j = random.randint(1, SIZE)
    # nash_i = random.randint(1,args.SIZE)
    # nash_j = random.randint(1,args.SIZE)
    payoff_matrix[('row', f"S{nash_i}"), ('column', f"S{nash_j}")] = \
        {'row': 1, 'column': 1}
    return payoff_matrix


# sizes = [5, 10, 20, 30, 50, 75, 100, 150, 250, 400, 700, 1000]
# gen(sizes)
# 315     99225
# 385     148225
# 445     198025
# 500     250000
# 550     302500
# 590     348100
# 630     396900
# 670     448900
# 705     497025
# 740     547600
# 775     600625
# 805     648025
# 835     697225
# 865     748225
# 895     801025
# 920     846400
# 950     902500
# 975     950625
# sizes = [315, 385, 445, 500, 550, 590, 630, 670, 705, 740, 775, 805, 835, 865, 895, 920, 950, 975]
# gen(sizes)
gen([5])








# SIZE = 4
# payoff_matrix = generate_random_matrix(SIZE)
# write_matrix(payoff_matrix, SIZE)
# def convert_to_payoffMatrix_key(i):
#     """
#     This function converts the elements of the set I in the pyomo model
#     (or elements of gameStatesForI) to the format of keys of the payoff matrix
#     of the game, i.e., ('p1','s1','p2','s2') is converted to (('p1','s1'),('p2','s2')) 
#     (see optModel.I for details)
#     """
#     gameState = []
#     done = 0
#     k1 = list(i)
#     while done == 0:
#         gameState.append(tuple(k1[0:2]))
#         # print('tuple(k1[0:2])', tuple(k1[0:2]))
#         del k1[0:2]
#         if len(k1) == 0:
#             done = 1
#     return tuple(gameState)



# payoff_matrix1 = {}
# # pick two random strategies
# i = random.sample(range(1,SIZE+1), 2)
# LB = 0
# UB = 20
# desired_equilibrium = (f"S{i[0]}", f"S{i[1]}")
# players_strategies = [f"S{i}" for i in range(1,SIZE+1)]
# # pupulate the payoff matrix with random payoffs
# for i in range(1,SIZE+1):
#     for j in range(1,SIZE+1):
#         payoff_matrix1[('row', f"S{i}"), ('column', f"S{j}")] = \
#             {'row': random.randint(-15,0), 'column': random.randint(-15,0)}
# players_names1 = ['row', 'column']
# #--- Create a pyomo model optModel ---
# optModel = ConcreteModel()
# #--- Define sets --- # Set of players
# optModel.P = Set(initialize = players_names1) 
# # optModel.I_row = Set(initialize = [tuple([k3 for k2 in k1 for k3 in k2] + ['row']) for k1 in payoff_matrix1.keys()])   
# # optModel.I_column = Set(initialize = [tuple([k3 for k2 in k1 for k3 in k2] + ['column']) for k1 in payoff_matrix1.keys()])
# # combine the sets I_row and I_column to get the set I
# optModel.I = Set(initialize = [tuple([k3 for k2 in k1 for k3 in k2] + ['row']) for k1 in payoff_matrix1.keys()] + [tuple([k3 for k2 in k1 for k3 in k2] + ['column']) for k1 in payoff_matrix1.keys()])
# #--- Define the variables --- 
# # optModel.y_row = Var(optModel.I_row, domain_type=IntegerSet, lb=LB, ub=UB)
# # optModel.y_column = Var(optModel.I_column, domain_type=IntegerSet, lb=LB, ub=UB)
# optModel.y = Var(optModel.I, domain=Integers, bounds=(LB, UB))
# #--- Define the objective function and constraints ----
# # Objective function
# optModel.objective_rule = Objective(rule = lambda optModel: 0, sense = maximize)
# # Constraint checking the best strategy of player p given the strategy of 
# # all other players 
# def NashCond_rule(optModel,*i):
#     # Convert the game state to the format of keys of the payoff matrix
#     l = list(i)
#     # l = ['row','s1','column','s2','row]
#     if (l[1], l[3]) == desired_equilibrium:
#         if l[-1] == 'row':
#             return optModel.y[i] >= 1
#             expr = [optModel.y[iterator] for iterator in [(l[0], row_strat, l[2], l[3], 'row') for row_strat in players_strategies] ]
#             return optModel.y[i] >= max(expr)
#             # constraints = [(optModel.y[i] >= optModel.y[iterator]) for iterator in [(l[0], row_strat, l[2], l[3], 'row') for row_strat in players_strategies]]
#             # return constraints
#         elif l[-1] == 'column':
#             return optModel.y[i] >= 1
#             return optModel.y[i] >= max([optModel.y[iterator] for iterator in [(l[0], l[1], l[2], column_strat, 'column') for column_strat in players_strategies]])
#             # constraints = [(optModel.y[i] >= optModel.y[iterator]) for iterator in [(l[0], l[1], l[2], column_strat, 'column') for column_strat in players_strategies]]
#             # return constraints

#     else:
#         if l[-1] == 'row':
#             return optModel.y[i] >= LB
#         elif l[-1] == 'column':
#             return optModel.y[i] >= LB
# optModel.NashCond = Constraint(optModel.I, rule=NashCond_rule)
# print('optModel.NashCond', optModel.NashCond)

# optimization_solver = 'gurobi'
# solverType = pyomoSolverCreator(optimization_solver)
# optModel.preprocess()

# #- Solve the optModel (tee=True shows the solver output) -
# try:
#     optSoln = solverType.solve(optModel,tee=False)
#     solverFlag = 'normal'
# except:
#     print ("WARNING! ", optimization_solver," failed. An alternative solver is tried")  

#     if optimization_solver.lower() == 'gurobi':
#         optimization_solver = 'cplex'
#     elif optimization_solver.lower() == 'cplex':
#         optimization_solver = 'gurobi'

#     # Try solving with the alternative solver
#     solverType = pyomoSolverCreator(optimization_solver)
#     try:
#         optSoln = solverType.solve(optModel,tee=False)
#         solverFlag = 'normal'
#     except:
#         solverFlag = 'solverError'
#         print('\nWARNING! The alternative solver failed. No solution was returned')








# import optlang
# import sympy
# model = optlang.Model(name=f'Original Model')
# variable_names = []

# # Add the variables
# for i in range(1,SIZE+1):
#     for j in range(1,SIZE+1):
#         var1 = optlang.Variable(f"a1_{i}_{j}", lb=0, ub=20, type='integer', problem=model)
#         var2 = optlang.Variable(f"a2_{i}_{j}", lb=0, ub=20, type='integer', problem=model) 
#         model.add(var1)
#         model.add(var2)
#         variable_names.append(f"a1_{i}_{j}")
#         variable_names.append(f"a2_{i}_{j}")

# # Add the objective function
# model.objective = optlang.Objective(expression=0, direction='min')


# # pick a random pair
# r1 = random.randint(1,SIZE)
# r2 = random.randint(1,SIZE)
# LB = 0



# for var in variable_names:
#     if var.endswith(f"_{r1}_{r2}"):
#         if var[:2] == 'a1':
#             # print("HIIII")
#             # constraint is a_ij >= max(a1_p,q) for all p 
#             # write a sympy expression for the constraint
#             # expression = model.variables[var] - sympy.Max(*[model.variables[f"{var[:2]}_{p}_{r2}"] for p in range(1,SIZE+1)])
#             # c = optlang.Constraint(expression, lb=0)
#             eps = 3
#             expressions = [(model.variables[var] - model.variables[f"{var[:2]}_{p}_{r2}"] - eps) for p in range(1,SIZE+1) if p != r1]
#             constraints = [optlang.Constraint(expression, lb=0) for expression in expressions]
#             for c in constraints:
#                 model.add(c)
#         elif var[:2] == 'a2':
#             expressions = [(model.variables[var] - model.variables[f"{var[:2]}_{r1}_{q}"] - eps) for q in range(1,SIZE+1) if q != r2]
#             constraints = [optlang.Constraint(expression, lb=0) for expression in expressions]
#             for c in constraints:
#                 model.add(c)
#     else:
#         # constraint is a_ij >= LB
#         c = optlang.Constraint(model.variables[var] - LB, lb=0)
#         model.add(c)


# # print model constraints
# for c in model.constraints:
#     print(c)


# model.optimize()
# print(model.status)
# payoff_matrix = {}
# for i in range(SIZE):
#     for j in range(SIZE):
#         # print the primal
#         print(f"a1_{i+1}_{j+1} = {model.variables[f'a1_{i+1}_{j+1}'].primal}")
#         payoff_matrix[('row', f"S{i+1}"), ('column', f"S{j+1}")] = \
#             {'row': model.variables[f"a1_{i+1}_{j+1}"].primal, 'column': model.variables[f"a2_{i+1}_{j+1}"].primal}

# print("r1, r2 = ", r1, r2)



# # payoff_matrix = generate_nash_equilibrium(SIZE)

# game_name = f"{SIZE} strategies"
# players_names = ['row','column']

# players_strategies = {}
# strategies = ['S' + str(i) for i in range(1,SIZE+1)]
# players_strategies['row'] = strategies
# players_strategies['column'] = strategies

# PD = game(game_name, players_names, players_strategies, payoff_matrix)

# # Define an instance of the NashEqFinder
# NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
# [Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
# print("Nash_equilibria optlangFindPure = ", Nash_equilibria)
# # order the equilibria
# Nash_equilibria.sort()