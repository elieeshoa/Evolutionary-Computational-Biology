# os.chdir(dir_path)
from NashEqFinder_manuscript import *
from game_manuscript import *

# --------- Helper Functions ------
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

# --------- Sample implementation ------
def run_prisoners_dilemma():
    print("\n-- Prisoner's Dilemma ---")
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
    NashEqFinderInst = NashEqFinder(PD, stdout_msgs=False, stdout_timing=False)
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # (self, original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells)
    show_matrix_original(game_payoff_matrix, game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Payoff Matrix", [])
    # Nash_equilibria is a list of [[('row', 'C'), ('column', 'C')], [('row', 'D'), ('column', 'D')]
    # but we need it to be a list of [(('row', 'C'), ('column', 'C')), (('row', 'D'), ('column', 'D'))]
    original_nasheq_cells = [(cell[0], cell[1]) for cell in Nash_equilibria]
    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','C'), ('column','C'))],
                                   strategies=['C', 'D'],
                                   original_nasheq_cells=original_nasheq_cells)
    # find the new equilibria
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()


def run_snowdrift():
    print ("\n-- Snowdrift's Dilemma ---")
    # # Pure strategy Nash eq = [(C,C), (D,D)]
    game_name = "Snowdrift's Dilemma"
    numberOfPlayers = 2
    players_names = ['row','column']
    players_strategies = {}
    players_strategies['row'] = ['C','D']
    players_strategies['column'] = ['C','D']
    payoff_matrix = {}
    payoff_matrix[(('row','C'),('column','C'))] = {'row':3,'column':3}
    payoff_matrix[(('row','C'),('column','D'))] = {'row':1,'column':5}
    payoff_matrix[(('row','D'),('column','C'))] = {'row':5,'column':1}
    payoff_matrix[(('row','D'),('column','D'))] = {'row':0,'column':0}
    # Define an instance of the game
    SD = game(game_name, players_names, players_strategies, payoff_matrix)
    # Define an instance of the NashEqFinder
    NashEqFinderInst = NashEqFinder(SD, stdout_msgs=False, stdout_timing=False)
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # (self, original_payoff_matrix, payoff_matrix, nash_equilibria, strategies, method, changed_cells)
    show_matrix_original(game_payoff_matrix, game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Payoff Matrix", [])
    # Nash_equilibria is a list of [[('row', 'C'), ('column', 'C')], [('row', 'D'), ('column', 'D')]
    # but we need it to be a list of [(('row', 'C'), ('column', 'C')), (('row', 'D'), ('column', 'D'))]
    original_nasheq_cells = [(cell[0], cell[1]) for cell in Nash_equilibria]
    print("original_nasheq_cells:", original_nasheq_cells)
    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row','C'), ('column','C'))],
                                   strategies=['C', 'D'],
                                   original_nasheq_cells=original_nasheq_cells)
    # find the new equilibria
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()


def run_5_by_5():
    """ 5 by 5 """
    SIZE = 5
    game_name = f"{SIZE} strategies"
    players_names = ['row','column']
    players_strategies = {}
    strategies = ['S' + str(i) for i in range(1,SIZE+1)]
    players_strategies['row'] = strategies
    players_strategies['column'] = strategies
    payoff_matrix = {}
    read_matrix(payoff_matrix, SIZE)
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    NashEqFinderInst = NashEqFinder(PD, optimization_solver='gurobi', stdout_msgs=False, stdout_timing=False)
    [Nash_equilibria, exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    original_nasheq_cells = [(cell[0], cell[1]) for cell in Nash_equilibria]
    i = 5
    j = 2
    print("original_nasheq_cells:", original_nasheq_cells)
    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row',f'S{i}'), ('column',f'S{j}'))],
                                  strategies=strategies,
                                  original_nasheq_cells=original_nasheq_cells)
    
if __name__ == "__main__":
    run_5_by_5()
    # run_prisoners_dilemma()