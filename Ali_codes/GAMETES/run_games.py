from NashEqFinder import *
from game import *
from generate_game import *
import datetime


for SIZE in [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]:
# for SIZE in [1000]:
    print(f"\n\n\n\n\n\n-------- SIZE: {SIZE} --------")
    print(f"Time before setting the game: {datetime.datetime.now()}")
    game_name = f"{SIZE} strategies"
    players_names = ['row','column']
    players_strategies = {}
    strategies = ['S' + str(i) for i in range(1,SIZE+1)]
    players_strategies['row'] = strategies
    players_strategies['column'] = strategies
    payoff_matrix = {}
    print(f"Time before reading the matrix: {datetime.datetime.now()}")
    read_matrix(payoff_matrix, SIZE)
    
    # Define an instance of the game
    print(f"Time before setting the game object: {datetime.datetime.now()}")
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    # Define an instance of the NashEqFinder
    print(f"Time before setting the first NashEqFinder: {datetime.datetime.now()}")
    NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    # print(f"Time before running the first optlangRun: {datetime.datetime.now()}")
    # [Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # print(f"Time before running the first show_matrix: {datetime.datetime.now()}")
    # show_matrix(game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Game called by optlangFindPure")
    # print("Nash_equilibria optlangFindPure = ", Nash_equilibria)

    # Pick a random pair of numbers
    i = random.randint(1, SIZE)
    j = random.randint(1, SIZE)

    print(f"Time before running newEquilibria: {datetime.datetime.now()}")
    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row',f'S{i}'), ('column',f'S{j}'))], strategies=strategies)




    