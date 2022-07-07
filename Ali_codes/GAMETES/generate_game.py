# create a 100 strategies names S1, S2, S3, ..., S100
import random


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
    print("Writing matrix to file")
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
        print(len(lines))
        print(len(lines[0].split()))
        for i in range(SIZE):
            for j in range(SIZE):
                payoff_matrix[('row', f"S{i+1}"), ('column', f"S{j+1}")] = \
                    {'row': int(lines[i].split()[2*j]), 'column': int(lines[i].split()[2*j+1])}

def gen():
    # pupulate the payoff matrix with random payoffs upper triangle
    for SIZE in [5, 10, 20, 30, 50, 75, 100, 150, 250, 400, 700, 1000]:
        payoff_matrix = generate_random_matrix(SIZE)
        write_matrix(payoff_matrix, SIZE)