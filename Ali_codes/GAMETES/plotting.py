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

def get_averages(string):
    sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    # averages = []
    # for i in range(len(sizes)):
    #     # Get the average time for each size
    #     # Get the lines of the string
    #     lines = string.split('\n')
    #     times = [(line.split('||')[1].split(': ')[1]) for line in lines[i*5:i*5+5]]
    #     # turn the string datetimes into seconds
    #     times = [dt.strptime(time, '%H:%M:%S.%f') for time in times]
    #     times = [time.second + time.microsecond/1000000 for time in times]
    #     print(times)
    #     average = sum(times)/len(times)
    #     averages.append(average)
    # print(averages)
    # return averages

    lines = string_summer_2023.split('\n')

    # Create a dictionary of times for every size.
    times = {}
    
    for line in lines:
        # Get the size from the line.
        size = line.split('size = ')[1].split(',')[0]

        # Get the time from the line.
        time = dt.strptime(line.split('Time: ')[1], '%H:%M:%S.%f')

        # Format the time as seconds.
        time = time.second + time.microsecond/1000000 + time.minute*60

        print(f"for size {size}, time is {time}")

        # Add the time to the dictionary of times for the given size.
        if size in times:
            times[size].append(time)
        else:
            times[size] = [time]


    # Calculate the average time for every size.
    average_times = {}
    for size in times:
        average_times[size] = sum(times[size])/len(times[size])

    occurences = {size: len(times[size]) for size in times}
    return average_times, occurences


def get_all_times(string):
    # returns a dictionary with the size as the key and a list of times as the value
    lines = string.split('\n')
    times = {}
    for line in lines:
        # Get the size from the line.
        size = line.split('size = ')[1].split(',')[0]

        # Get the time from the line.
        time = dt.strptime(line.split('Time: ')[1], '%H:%M:%S.%f')

        # Format the time as seconds.
        time = time.second + time.microsecond/1000000 + time.minute*60

        # Add the time to the dictionary of times for the given size.
        if size in times:
            times[size].append(time)
        else:
            times[size] = [time]

    return times

def plot_averages(string, quadratic=True):
    # sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    # averages = get_averages(string)

    averages_dict, occurences_dict = get_averages(string)
    sizes = [int(size) for size in averages_dict.keys()]
    sizes.sort()
    averages = [averages_dict[str(size)] for size in sizes]

    for size, average in zip(sizes, averages):
        print(f"size {size} has square {size**2} and average {average} over {occurences_dict[str(size)]} occurences")




    
    # plt.plot(sizes, [size**2 for size in sizes])
    # fit a quadratic curve to the data
    if quadratic:
        fit = np.polyfit(sizes, averages, 2)
        fit_fn = np.poly1d(fit)
        # draw the fit curve fit_fn(sizes) for a range from 0 to 1000
        x = np.linspace(5, max(sizes), 1000)
        # plt.plot(x, fit_fn(x), '--k', label='Quadratic fit', linewidth=1, color='dodgerblue', alpha=0.8)
        plt.plot(x, fit_fn(x), '--k', label='Quadratic fit', linewidth=1.2, color='#00A6B7', alpha=0.8)
        plt.xticks(np.arange(0, max(sizes)+1, 100))
        # give them 45 degree angle
        plt.xticks(rotation=45)
        plt.yticks(np.arange(0, max(averages)+1, 5))

        # print r^2
        r2 = np.corrcoef(sizes, averages)[0, 1]**2
        print(f"R^2 = {r2}")
        
        # plt.plot(sizes, averages, label='Average Time', linewidth=1, color='red')
        # make the color FF6805
        # plt.plot(sizes, averages, label='Average Time', linewidth=1.5, color='#FF7F00')
        # label the points with how many times they occured

        # plot the data points
        # plt.scatter(sizes, averages, color='red', s=15)
        plt.scatter(sizes, averages, color='#FF7F00', s=25, zorder=3)
        plt.xlabel('Number of strategies ($n$)', fontsize=15)

        plt.ylabel('Average CPU time (seconds)', fontsize=15)
        plt.title('Average time to change payoffs to new Nash equilibrium')
        # make the y numbers integers and not 1.0e7
        plt.ticklabel_format(style='plain')
        # for size, average in zip(sizes, averages):
        #     plt.text(size, average, occurences_dict[str(size)], fontsize=8, color='black', alpha=0.8)
        # add thin gridlines
        plt.grid(linewidth=0.4)

        # legends
        plt.legend(loc='upper left')
        plt.gcf().set_size_inches(7/1.1, 6/1.1)
        plt.tight_layout()


        plt.savefig('large_games_quadratic.png', dpi=300)

        plt.show()

    else: # linear
        sizes = [size**2 for size in sizes]
        fit = np.polyfit(sizes, averages, 1)
        fit_fn = np.poly1d(fit)
        x = np.linspace(25, max(sizes), 1000)
        plt.xticks(np.arange(0, max(sizes)+1, 100000))
        # give them 45 degree angle
        plt.xticks(rotation=45)
        plt.yticks(np.arange(0, max(averages)+1, 5))

        # plt.plot(x, fit_fn(x), '--k', label='Linear fit', linewidth=1, color='dodgerblue', alpha=0.8)
        plt.plot(x, fit_fn(x), '--k', label='Linear fit', linewidth=1.2, color='#00A6B7', alpha=0.8)
        # plt.plot(sizes, averages, label='Average Time', linewidth=1, color='red')
        # make the color FF6805
        # plt.plot(sizes, averages, label='Average Time', linewidth=1.5, color='#FF7F00')
        # plt.scatter(sizes, averages, color='red', s=15)
        plt.scatter(sizes, averages, color='#FF7F00', s=25, zorder=3)
        plt.xlabel('Size of game ($n^2$)', fontsize=15)
        for size, average in zip(sizes, averages):
            plt.text(size, average, occurences_dict[str(int(np.sqrt(size)))], fontsize=10, color='black', alpha=0.8, zorder=4)
            
        # calculate R^2
        residuals = np.polyfit(sizes, averages, 1)
        r_squared = 1 - (np.var(residuals) / np.var(averages))
        r2 = np.corrcoef(sizes, averages)[0, 1]**2
        print(f"R^2 is {r_squared}")
        print(f"R^2 is {r2}")
        # Calculate residuals
        residuals = averages - fit_fn(sizes)

        # Calculate total sum of squares
        total_sum_squares = np.sum((averages - np.mean(averages))**2)

        # Calculate residual sum of squares
        residual_sum_squares = np.sum(residuals**2)

        # Calculate R^2
        r_squared = 1 - (residual_sum_squares / total_sum_squares)

        print("R^2:", r_squared)

       
        



        plt.ylabel('Average CPU time (seconds)', fontsize=15)
        plt.title('Average time to change payoffs to new Nash equilibrium')
        # make the y numbers integers and not 1.0e7
        plt.ticklabel_format(style='plain')
        # add thin gridlines
        plt.grid(linewidth=0.4)

        # legends
        plt.legend(loc='upper left')
        # set size of plot to be 7 by 6
        plt.gcf().set_size_inches(7/1.1, 6/1.1)
        plt.tight_layout()


        plt.savefig('large_games_linear.png', dpi=300)

        plt.show()

def plot_medians(string):
    # plot the values of the medians and with error bars
    # sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    d = get_all_times(string)
    sizes = [int(size) for size in d.keys()]
    sizes.sort()
    medians = np.array([np.median(d[str(size)]) for size in sizes])
    

    # calculate the standard deviation
    # stds = [np.std(d[str(size)]) for size in sizes]

    # calculate the 25th and 75th percentile
    q1s = np.array([np.percentile(d[str(size)], 25) for size in sizes])
    q3s = np.array([np.percentile(d[str(size)], 75) for size in sizes])

    for size, median, q1, q3 in zip(sizes, medians, q1s, q3s):
        print(f"size {size} has median {median} and q1 = {q1} and q3 = {q3}  over {len(d[str(size)])} occurences")

    sizes = [size**2 for size in sizes]
    # plot the fit 
    fit = np.polyfit(sizes, medians, 1)
    fit_fn = np.poly1d(fit)
    x = np.linspace(25, max(sizes), 1000)
    plt.xticks(np.arange(0, max(sizes)+1, 100000))
    # give them 45 degree angle
    plt.xticks(rotation=45)
    
    # draw the fit line
    plt.plot(x, fit_fn(x), '--k', label='Linear fit', linewidth=1.2, color='#00A6B7', alpha=0.8)
    # print r^2
    r2 = np.corrcoef(sizes, medians)[0, 1]**2
    print(f"R^2 = {r2}")





    # scatter the medians
    plt.scatter(sizes, medians, color='red', s=15)
    # plot the error bars
    # plt.errorbar(sizes, medians, yerr=stds, fmt='o', color='#FF7F00', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.errorbar(sizes, medians, yerr=[medians-q1s, q3s-medians], fmt='o', color='#FF7F00', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.xlabel('Size of game ($n^2$)', fontsize=15)
    plt.ylabel('CPU time (seconds)', fontsize=15)
    plt.title('Median time to change payoffs to new Nash equilibrium')
    # make the y numbers integers and not 1.0e7
    plt.ticklabel_format(style='plain')
    # add thin gridlines
    plt.grid(linewidth=0.4)

    # set size of plot to be 7 by 6
    plt.gcf().set_size_inches(7/1.1, 6/1.1)
    plt.tight_layout()

    plt.savefig('large_games_median.png', dpi=300)

    plt.show()




# load the file "log Summer 2023.txt" into a string
string_summer_2023 = open('log Summer 2023.txt', 'r').read()
# string_summer_2023 = open('log2.txt', 'r').read()
# string_summer_2023 = open('log Summer 2023 2.txt', 'r').read()
# string_summer_2023 = open('log_eristwo.txt', 'r').read()
# string_summer_2023 = open('log_original_timer.txt', 'r').read()
# string_summer_2023 = open('log_cplex.txt', 'r').read()
string_summer_2023 = open('log_eris/log_eristwo.txt', 'r').read()

# string_summer_2023_alt = open('log_eristwo.txt', 'r').read()
# # remove all lines with "Original" or "CPLEX" in them
# string_summer_2023_alt = '\n'.join([line for line in string_summer_2023_alt.split('\n') if 'Original' not in line and 'CPLEX' not in line])
# string_summer_2023_org = open('log_eris/log_eristwo.txt', 'r').read()
# # remove all lines with "CPLEX" in them
# string_summer_2023_org = '\n'.join([line for line in string_summer_2023_org.split('\n') if 'CPLEX' not in line])
# # combine the two strings
# string_summer_2023 = string_summer_2023_alt + "\n" + string_summer_2023_org
# # print(string_summer_2023)

# remove all lines with "Original" in them
# string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'Original' not in line])

# get all lines with "Original" in them
# string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'Original' in line])

# remove the lines with "CPLEX" in them
string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'CPLEX' not in line])

# remove lines without "RANDOMLY" in them
# string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'RANDOMLY' in line])

# remove the lines that don't have "CPLEX" in them
# string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'CPLEX' in line])

# write all lines with "Original" in them to a file called "log_original_timer.txt"
# with open('log_original_timer.txt', 'a') as f:
#     f.write('\n'.join([line for line in string_summer_2023.split('\n') if 'Original' in line]))

# plot_averages(string_summer_2023, quadratic=False)
# plot_averages(string_summer_2023, quadratic=False)
# plot_averages(string_summer_2023, quadratic=True)
plot_medians(string_summer_2023)


# 58627   ee869   RUN   normal     eris2n4     hn008       job_740    Jun 29 00:54
# 58632   ee869   PEND  normal     eris2n4                 job_895    Jun 29 00:54
# 58633   ee869   PEND  normal     eris2n4                 job_920    Jun 29 00:54
# 58634   ee869   PEND  normal     eris2n4                 job_950    Jun 29 00:54
# 58635   ee869   PEND  normal     eris2n4                 job_975    Jun 29 00:54
# 61369   ee869   PEND  normal     eris2n4                 job_20     Jun 29 12:22
# 61370   ee869   PEND  normal     eris2n4                 job_315    Jun 29 12:22
# 61371   ee869   PEND  normal     eris2n4                 job_385    Jun 29 12:22
# 61372   ee869   PEND  normal     eris2n4                 job_445    Jun 29 12:22
# 61373   ee869   PEND  normal     eris2n4                 job_500    Jun 29 12:22
# 61374   ee869   PEND  normal     eris2n4                 job_550    Jun 29 12:22
# 61375   ee869   PEND  normal     eris2n4                 job_590    Jun 29 12:22
# 61376   ee869   PEND  normal     eris2n4                 job_630    Jun 29 12:22
# 61377   ee869   PEND  normal     eris2n4                 job_670    Jun 29 12:22
# 61378   ee869   PEND  normal     eris2n4                 job_705    Jun 29 12:22
# 61379   ee869   PEND  normal     eris2n4                 job_740    Jun 29 12:22
# 61380   ee869   PEND  normal     eris2n4                 job_775    Jun 29 12:22
# 61381   ee869   PEND  normal     eris2n4                 job_805    Jun 29 12:22
# 61382   ee869   PEND  normal     eris2n4                 job_835    Jun 29 12:22
# 61383   ee869   PEND  normal     eris2n4                 job_865    Jun 29 12:22
# 61384   ee869   PEND  normal     eris2n4                 job_895    Jun 29 12:22
# 61385   ee869   PEND  normal     eris2n4                 job_920    Jun 29 12:22
# 61386   ee869   PEND  normal     eris2n4                 job_950    Jun 29 12:22
# 61387   ee869   PEND  normal     eris2n4                 job_975    Jun 29 12:22

# bkill 61369
# bkill 61370
# bkill 61371
# bkill 61372
# bkill 61373
# bkill 61374
# bkill 61375
# bkill 61376
# bkill 61377
# bkill 61378
# bkill 61379
# bkill 61380
# bkill 61381
# bkill 61382
# bkill 61383
# bkill 61384
# bkill 61385
# bkill 61386
# bkill 61387
