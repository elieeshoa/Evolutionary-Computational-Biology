from datetime import datetime as dt
import matplotlib.pyplot as plt
import numpy as np


string = \
"""size = 5, 2c iteration 0 with desired nach equilibria: [(('row', 'S2'), ('column', 'S1'))] || Time: 0:00:00.089607
size = 5, 2c iteration 1 with desired nach equilibria: [(('row', 'S2'), ('column', 'S1'))] || Time: 0:00:00.008921
size = 5, 2c iteration 2 with desired nach equilibria: [(('row', 'S2'), ('column', 'S1'))] || Time: 0:00:00.005190
size = 5, 2c iteration 3 with desired nach equilibria: [(('row', 'S2'), ('column', 'S1'))] || Time: 0:00:00.005113
size = 5, 2c iteration 4 with desired nach equilibria: [(('row', 'S2'), ('column', 'S1'))] || Time: 0:00:00.004982
size = 10, 2c iteration 0 with desired nach equilibria: [(('row', 'S5'), ('column', 'S5'))] || Time: 0:00:00.002025
size = 10, 2c iteration 1 with desired nach equilibria: [(('row', 'S5'), ('column', 'S5'))] || Time: 0:00:00.004903
size = 10, 2c iteration 2 with desired nach equilibria: [(('row', 'S5'), ('column', 'S5'))] || Time: 0:00:00.010714
size = 10, 2c iteration 3 with desired nach equilibria: [(('row', 'S5'), ('column', 'S5'))] || Time: 0:00:00.005935
size = 10, 2c iteration 4 with desired nach equilibria: [(('row', 'S5'), ('column', 'S5'))] || Time: 0:00:00.005876
size = 20, 2c iteration 0 with desired nach equilibria: [(('row', 'S2'), ('column', 'S13'))] || Time: 0:00:00.002486
size = 20, 2c iteration 1 with desired nach equilibria: [(('row', 'S2'), ('column', 'S13'))] || Time: 0:00:00.010239
size = 20, 2c iteration 2 with desired nach equilibria: [(('row', 'S2'), ('column', 'S13'))] || Time: 0:00:00.007542
size = 20, 2c iteration 3 with desired nach equilibria: [(('row', 'S2'), ('column', 'S13'))] || Time: 0:00:00.008141
size = 20, 2c iteration 4 with desired nach equilibria: [(('row', 'S2'), ('column', 'S13'))] || Time: 0:00:00.006896
size = 30, 2c iteration 0 with desired nach equilibria: [(('row', 'S26'), ('column', 'S18'))] || Time: 0:00:00.051712
size = 30, 2c iteration 1 with desired nach equilibria: [(('row', 'S26'), ('column', 'S18'))] || Time: 0:00:00.022506
size = 30, 2c iteration 2 with desired nach equilibria: [(('row', 'S26'), ('column', 'S18'))] || Time: 0:00:00.010083
size = 30, 2c iteration 3 with desired nach equilibria: [(('row', 'S26'), ('column', 'S18'))] || Time: 0:00:00.011068
size = 30, 2c iteration 4 with desired nach equilibria: [(('row', 'S26'), ('column', 'S18'))] || Time: 0:00:00.015029
size = 50, 2c iteration 0 with desired nach equilibria: [(('row', 'S27'), ('column', 'S7'))] || Time: 0:00:00.052816
size = 50, 2c iteration 1 with desired nach equilibria: [(('row', 'S27'), ('column', 'S7'))] || Time: 0:00:00.299637
size = 50, 2c iteration 2 with desired nach equilibria: [(('row', 'S27'), ('column', 'S7'))] || Time: 0:00:00.047638
size = 50, 2c iteration 3 with desired nach equilibria: [(('row', 'S27'), ('column', 'S7'))] || Time: 0:00:00.055767
size = 50, 2c iteration 4 with desired nach equilibria: [(('row', 'S27'), ('column', 'S7'))] || Time: 0:00:00.064685
size = 100, 2c iteration 0 with desired nach equilibria: [(('row', 'S33'), ('column', 'S88'))] || Time: 0:00:00.181271
size = 100, 2c iteration 1 with desired nach equilibria: [(('row', 'S33'), ('column', 'S88'))] || Time: 0:00:00.164494
size = 100, 2c iteration 2 with desired nach equilibria: [(('row', 'S33'), ('column', 'S88'))] || Time: 0:00:00.156647
size = 100, 2c iteration 3 with desired nach equilibria: [(('row', 'S33'), ('column', 'S88'))] || Time: 0:00:00.169929
size = 100, 2c iteration 4 with desired nach equilibria: [(('row', 'S33'), ('column', 'S88'))] || Time: 0:00:00.215340
size = 250, 2c iteration 0 with desired nach equilibria: [(('row', 'S77'), ('column', 'S149'))] || Time: 0:00:00.689703
size = 250, 2c iteration 1 with desired nach equilibria: [(('row', 'S77'), ('column', 'S149'))] || Time: 0:00:00.852490
size = 250, 2c iteration 2 with desired nach equilibria: [(('row', 'S77'), ('column', 'S149'))] || Time: 0:00:00.910095
size = 250, 2c iteration 3 with desired nach equilibria: [(('row', 'S77'), ('column', 'S149'))] || Time: 0:00:00.947000
size = 250, 2c iteration 4 with desired nach equilibria: [(('row', 'S77'), ('column', 'S149'))] || Time: 0:00:00.914404
size = 400, 2c iteration 0 with desired nach equilibria: [(('row', 'S123'), ('column', 'S332'))] || Time: 0:00:02.228152
size = 400, 2c iteration 1 with desired nach equilibria: [(('row', 'S123'), ('column', 'S332'))] || Time: 0:00:02.548280
size = 400, 2c iteration 2 with desired nach equilibria: [(('row', 'S123'), ('column', 'S332'))] || Time: 0:00:03.533587
size = 400, 2c iteration 3 with desired nach equilibria: [(('row', 'S123'), ('column', 'S332'))] || Time: 0:00:02.408298
size = 400, 2c iteration 4 with desired nach equilibria: [(('row', 'S123'), ('column', 'S332'))] || Time: 0:00:03.090571
size = 700, 2c iteration 0 with desired nach equilibria: [(('row', 'S622'), ('column', 'S495'))] || Time: 0:00:10.831137
size = 700, 2c iteration 1 with desired nach equilibria: [(('row', 'S622'), ('column', 'S495'))] || Time: 0:00:11.929054
size = 700, 2c iteration 2 with desired nach equilibria: [(('row', 'S622'), ('column', 'S495'))] || Time: 0:00:14.781851
size = 700, 2c iteration 3 with desired nach equilibria: [(('row', 'S622'), ('column', 'S495'))] || Time: 0:00:12.353810
size = 700, 2c iteration 4 with desired nach equilibria: [(('row', 'S622'), ('column', 'S495'))] || Time: 0:00:12.368400
size = 1000, 2c iteration 0 with desired nach equilibria: [(('row', 'S441'), ('column', 'S350'))] || Time: 0:00:27.779875
size = 1000, 2c iteration 1 with desired nach equilibria: [(('row', 'S441'), ('column', 'S350'))] || Time: 0:00:24.011225
size = 1000, 2c iteration 2 with desired nach equilibria: [(('row', 'S441'), ('column', 'S350'))] || Time: 0:00:22.359567
size = 1000, 2c iteration 3 with desired nach equilibria: [(('row', 'S441'), ('column', 'S350'))] || Time: 0:00:24.437762
size = 1000, 2c iteration 4 with desired nach equilibria: [(('row', 'S441'), ('column', 'S350'))] || Time: 0:00:36.910838"""


# For each size, we run 5 iterations of the algorithm, and we take the average time
# of the 5 iterations. We then plot the average time against the size of the game.
# We also plot the theoretical time complexity of the algorithm, which is O(n^2).
# We see that the algorithm is indeed O(n^2), as the plot is a straight line.

def get_averages():
    sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    averages = []
    for i in range(len(sizes)):
        # Get the average time for each size
        # Get the lines of the string 'string' defined above
        lines = string.split('\n')
        times = [(line.split('||')[1].split(': ')[1]) for line in lines[i*5:i*5+5]]
        # turn the string datetimes into seconds
        times = [dt.strptime(time, '%H:%M:%S.%f') for time in times]
        times = [time.second + time.microsecond/1000000 for time in times]
        print(times)
        average = sum(times)/len(times)
        averages.append(average)
    print(averages)
    return averages

def plot_averages():
    sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    averages = get_averages()
    
    # plt.plot(sizes, [size**2 for size in sizes])
    # fit a quadratic curve to the data
    fit = np.polyfit(sizes, averages, 2)
    fit_fn = np.poly1d(fit)
    # draw the fit curve fit_fn(sizes) for a range from 0 to 1000
    x = np.linspace(0, 1000, 1000)
    plt.plot(x, fit_fn(x), '--k', label='Quadratic fit', linewidth=1, color='dodgerblue', alpha=0.8)
    
    plt.plot(sizes, averages, label='Average Time', linewidth=1, color='red')
    # plot the data points
    plt.scatter(sizes, averages, color='red', s=15)



    plt.xlabel('Size of game ($n$)')
    plt.ylabel('Average CPU time (seconds)')
    plt.title('Average time to change payoffs to new Nash equilibrium')
    # add thin gridlines
    plt.grid(linewidth=0.2)

    # legends
    plt.legend(loc='upper left')


    plt.savefig('large_games.png', dpi=300)

    plt.show()

plot_averages()
