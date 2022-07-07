from NashEqFinder import *
from game import *
from generate_game import *


for SIZE in [5, 10, 20, 30, 50, 75, 100, 150, 250, 400, 700, 1000]:
    game_name = f"{SIZE} strategies"
    players_names = ['row','column']
    players_strategies = {}
    strategies = ['S' + str(i) for i in range(1,SIZE+1)]
    players_strategies['row'] = strategies
    players_strategies['column'] = strategies
    payoff_matrix = {}
    read_matrix(payoff_matrix, SIZE)
    
    # Define an instance of the game
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(PD, stdout_msgs = True)
    [Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    show_matrix(game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Game called by optlangFindPure")
    print("Nash_equilibria optlangFindPure = ", Nash_equilibria)

    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','S15'), ('column','S10'))], strategies=strategies)



    