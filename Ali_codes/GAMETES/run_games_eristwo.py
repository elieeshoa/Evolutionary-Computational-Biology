import argparse
from NashEqFinder import *
from game import *
import datetime


def generate_random_matrix(SIZE):
    payoff_matrix = {}
    for i in range(1,SIZE+1):
        for j in range(1,SIZE+1):
            payoff_matrix[('row', f"S{i}"), ('column', f"S{j}")] = \
                {'row': random.randint(-30,0), 'column': random.randint(-30,0)}
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
                
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, help="The size of the game.")
    args = parser.parse_args()

    print(f"\n\n\n\n\n\n-------- SIZE: {args.size} --------")
    print(f"Time before setting the game: {datetime.datetime.now()}")
    game_name = f"{args.size} strategies"
    players_names = ['row','column']
    players_strategies = {}
    strategies = ['S' + str(i) for i in range(1,args.size+1)]
    players_strategies['row'] = strategies
    players_strategies['column'] = strategies
    payoff_matrix = {}
    print(f"Time before writing the matrix: {datetime.datetime.now()}")
    payoff_matrix = generate_random_matrix(args.size)
    # set one random cell to (1, 1)
    nash_i = random.randint(1,args.size)
    nash_j = random.randint(1,args.size)
    payoff_matrix[('row', f"S{nash_i}"), ('column', f"S{nash_j}")] = \
        {'row': 1, 'column': 1}
    write_matrix(payoff_matrix, args.size)
    print(f"Time before reading the matrix: {datetime.datetime.now()}")
    read_matrix(payoff_matrix, args.size)

    # Define an instance of the game
    print(f"Time before setting the game object: {datetime.datetime.now()}")
    PD = game(game_name, players_names, players_strategies, payoff_matrix)
    # Define an instance of the NashEqFinder
    print(f"Time before setting the first NashEqFinder: {datetime.datetime.now()}")
    NashEqFinderInst = NashEqFinder(PD, optimization_solver='gurobi')
    # print(f"Time before running the first optlangRun: {datetime.datetime.now()}")
    #[Nash_equilibria,exit_flag, game_payoff_matrix] = NashEqFinderInst.optlangRun()
    # print(f"Time before running the first show_matrix: {datetime.datetime.now()}")
    # show_matrix(game_payoff_matrix, Nash_equilibria, players_strategies['row'], "Original Game called by optlangFindPure")
    # print("Nash_equilibria optlangFindPure = ", Nash_equilibria)

    # Pick a random pair of numbers
    i = random.randint(1, args.size)
    j = random.randint(1, args.size)


    original_nasheq_cells = [(('row', f'S{nash_i}'), ('column', f'S{nash_j}'))]

    print(f"Time before running newEquilibria: {datetime.datetime.now()}")
    NashEqFinderInst.newEquilibria(nasheq_cells=[(('row',f'S{i}'), ('column',f'S{j}'))], strategies=strategies, original_nasheq_cells=original_nasheq_cells)

if __name__ == "__main__":
    main()
            