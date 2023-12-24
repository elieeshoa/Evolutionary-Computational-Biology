from datetime import datetime as dt
import matplotlib.pyplot as plt
import numpy as np



# For each size, we run 5 iterations of the algorithm, and we take the average time
# of the 5 iterations. We then plot the average time against the size of the game.
# We also plot the theoretical time complexity of the algorithm, which is O(n^2).
# We see that the algorithm is indeed O(n^2), as the plot is a straight line.

def get_averages(string):
    sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
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

    # Create a dictionary of times for every size.
    times = {}
    lines = string.split('\n')
    for line in lines:
        # Get the size from the line.
        size = line.split('size = ')[1].split(',')[0]

        # Get the time from the line.
        time = dt.strptime(line.split('Time: ')[1], '%H:%M:%S.%f')

        # Format the time as seconds.
        time = time.second + time.microsecond/1000000 + time.minute*60

        # print(f"for size {size}, time is {time}")

        # Add the time to the dictionary of times for the given size.
        if size in times:
            times[size].append(time)
        else:
            times[size] = [time]

    print(times)
    occurences = {size: len(times[size]) for size in times}
    

    # calculate the standard deviation
    # stds = [np.std(d[str(size)]) for size in sizes]

    # calculate the 25th and 75th percentile
    q1s = np.array([np.percentile(d[str(size)], 25) for size in sizes])
    q3s = np.array([np.percentile(d[str(size)], 75) for size in sizes])

    for size, median, q1, q3 in zip(sizes, medians, q1s, q3s):
        print(f"size {size} has median {median} and q1 = {q1} and q3 = {q3}  over {len(d[str(size)])} occurences")
        
    sizes = [int(size**2) for size in sizes]
    # plot the fit 
    fit = np.polyfit(sizes, medians, 1)
    fit_fn = np.poly1d(fit)
    x = np.linspace(25, max(sizes), 1000)
    plt.xticks(np.arange(0, max(sizes)+1, 100000))
    # give them 45 degree angle
    plt.xticks(rotation=45)
    
    # draw the fit line
    plt.plot(x, fit_fn(x), '--', label='Linear fit', linewidth=1.2, color='#00A6B7', alpha=1, zorder=2)
    # print r^2
    r2 = np.corrcoef(sizes, medians)[0, 1]**2
    print(f"R^2 = {r2}")

    # write r2 on the plot in the center
    plt.text(0.5, 0.85, f"$R^2 = {r2:.2f}$", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, fontsize=15)

    for size, med in zip(sizes, medians):
        plt.text(size, med, occurences[str(int(np.sqrt(size)))], fontsize=10, color='black', alpha=0.8, zorder=4)

    


    # scatter the medians
    plt.scatter(sizes, medians, color='red', s=15)
    # plot the error bars
    # plt.errorbar(sizes, medians, yerr=[medians-q1s, q3s-medians], fmt='o', color='#FF7F00', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.errorbar(sizes, medians, yerr=[medians-q1s, q3s-medians], fmt='o', color='r', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.xlabel('Size of the game ($n^2$)', fontsize=12)
    plt.ylabel('CPU runtime (seconds)', fontsize=12)
    plt.title('Median CPU Runtime over 10 Runs for Finding Alternative Solutions', fontsize=12, pad=10)
    # make the y numbers integers and not 1.0e7
    plt.ticklabel_format(style='plain')
    plt.ylim(bottom=-18)
    # add thin gridlines
    plt.grid(linewidth=0.4)

    # set size of plot to be 7 by 6
    plt.gcf().set_size_inches(7/1, 6/1)
    plt.tight_layout()
    plt.legend(loc='upper left', fontsize=12)

    plt.savefig('large_games_median.png', dpi=300)

    plt.show()




# load the file "log Summer 2023.txt" into a string
string_summer_2023 = open('log_eris/log_eristwo_november_2023.txt', 'r').read()

# remove the lines with "CPLEX" in them
string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'CPLEX' not in line])
# remove the lines with the word "mistake" in them  
string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'mistake' not in line])
# remove the empty lines
string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if line != ''])

string_summer_2023 = \
"""size = 775, date: 2023-11-21 00:07:13.136288, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:01:21.187628
size = 705, date: 2023-11-21 00:07:24.277801, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:01:08.561076
size = 740, date: 2023-11-21 00:07:59.834118, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:01:27.501620
size = 805, date: 2023-11-21 00:10:18.943677, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:01:36.778754
size = 835, date: 2023-11-21 00:11:05.775043, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:01:50.175007
size = 865, date: 2023-11-21 00:17:57.201387, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:01:55.291947
size = 895, date: 2023-11-21 00:18:34.293389, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:02:26.591452
size = 920, date: 2023-11-21 00:21:09.036552, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:02:41.974890
size = 950, date: 2023-11-21 00:25:53.534510, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:02:39.519503
size = 975, date: 2023-11-21 00:28:58.775918, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:03:11.112620
size = 1000, date: 2023-11-21 00:29:29.361624, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S612'), ('column', 'S226'))] || Time: 0:03:11.444593
size = 740, date: 2023-11-21 01:23:44.604753, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:00:54.517719
size = 705, date: 2023-11-21 01:24:25.370294, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:22.435107
size = 775, date: 2023-11-21 01:28:25.282209, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:40.022907
size = 805, date: 2023-11-21 01:37:41.040306, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:00:29.614124
size = 835, date: 2023-11-21 01:53:01.918296, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:54.318424
size = 865, date: 2023-11-21 01:59:26.801574, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:00:44.880167
size = 895, date: 2023-11-21 02:12:37.652715, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:00:46.648811
size = 705, date: 2023-11-21 02:14:45.576048, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:21.723061
size = 740, date: 2023-11-21 02:22:33.868433, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:01:00.624078
size = 920, date: 2023-11-21 02:23:22.466030, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:00:49.856399
size = 775, date: 2023-11-21 02:32:46.756732, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:37.422717
size = 950, date: 2023-11-21 02:34:32.850730, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:01:06.953861
size = 975, date: 2023-11-21 02:42:39.010178, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:02:01.500002
size = 805, date: 2023-11-21 02:45:52.060387, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:01:07.242313
size = 1000, date: 2023-11-21 02:55:06.760246, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S612'), ('column', 'S226'))] || Time: 0:01:36.546589
size = 705, date: 2023-11-21 03:14:01.046708, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:19.967409
size = 835, date: 2023-11-21 03:14:29.244446, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:40.014828
size = 865, date: 2023-11-21 03:17:38.246907, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:00:50.167225
size = 740, date: 2023-11-21 03:22:35.994228, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:00:32.406148
size = 775, date: 2023-11-21 03:36:35.094882, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:42.207329
size = 895, date: 2023-11-21 03:46:48.716957, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:00:52.366877
size = 805, date: 2023-11-21 04:00:10.614610, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:00:30.770553
size = 920, date: 2023-11-21 04:04:36.780679, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:01:53.562219
size = 705, date: 2023-11-21 04:07:23.108484, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:21.331959
size = 950, date: 2023-11-21 04:16:41.327635, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:00:56.397762
size = 740, date: 2023-11-21 04:25:01.215663, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:00:57.394584
size = 975, date: 2023-11-21 04:27:48.493550, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:01:13.938631
size = 835, date: 2023-11-21 04:27:53.555715, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:40.263209
size = 865, date: 2023-11-21 04:37:47.427966, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:00:40.932447
size = 775, date: 2023-11-21 04:41:41.167364, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:39.814275
size = 1000, date: 2023-11-21 04:42:38.932388, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S612'), ('column', 'S226'))] || Time: 0:01:11.174716
size = 705, date: 2023-11-21 04:58:42.174964, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:26.692385
size = 895, date: 2023-11-21 05:09:38.198202, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:00:53.611622
size = 805, date: 2023-11-21 05:12:25.226714, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:00:38.297397
size = 740, date: 2023-11-21 05:26:05.987463, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:00:27.743443
size = 920, date: 2023-11-21 05:38:27.044467, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:01:21.630702
size = 835, date: 2023-11-21 05:40:44.956687, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:51.349397
size = 775, date: 2023-11-21 05:48:40.023956, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:42.619258
size = 705, date: 2023-11-21 05:53:32.973721, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S383'), ('column', 'S327'))] || Time: 0:00:26.538133
size = 950, date: 2023-11-21 05:57:57.517502, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:00:49.094913
size = 865, date: 2023-11-21 06:00:47.241113, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:00:36.988614
size = 975, date: 2023-11-21 06:18:19.217432, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:01:51.603950
size = 805, date: 2023-11-21 06:20:32.109487, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:00:40.835082
size = 740, date: 2023-11-21 06:24:14.602642, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S244'), ('column', 'S249'))] || Time: 0:00:32.566461
size = 1000, date: 2023-11-21 06:36:35.291301, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S612'), ('column', 'S226'))] || Time: 0:01:08.361059
size = 895, date: 2023-11-21 06:43:35.361416, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:00:56.596568
size = 775, date: 2023-11-21 06:49:47.698592, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S341'), ('column', 'S222'))] || Time: 0:00:43.436877
size = 835, date: 2023-11-21 07:03:40.119551, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:42.966533
size = 865, date: 2023-11-21 07:20:50.533242, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:01:02.062527
size = 920, date: 2023-11-21 07:22:37.369669, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:01:16.597332
size = 805, date: 2023-11-21 07:30:30.377932, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S521'), ('column', 'S194'))] || Time: 0:00:46.244112
size = 950, date: 2023-11-21 07:43:41.011209, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:00:49.910528
size = 975, date: 2023-11-21 08:09:59.624919, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:01:00.221096
size = 895, date: 2023-11-21 08:18:26.857021, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:00:35.866531
size = 835, date: 2023-11-21 08:25:54.754230, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S771'), ('column', 'S724'))] || Time: 0:00:37.703223
size = 1000, date: 2023-11-21 08:32:40.440067, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S612'), ('column', 'S226'))] || Time: 0:00:18.690110
size = 865, date: 2023-11-21 08:42:35.115246, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S341'), ('column', 'S782'))] || Time: 0:00:45.058543
size = 920, date: 2023-11-21 09:00:58.795123, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:01:39.059494
size = 950, date: 2023-11-21 09:26:23.166980, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:00:52.715593
size = 895, date: 2023-11-21 09:51:44.524010, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S820'), ('column', 'S682'))] || Time: 0:01:09.132952
size = 975, date: 2023-11-21 09:59:18.262786, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:01:18.064225
size = 920, date: 2023-11-21 10:34:46.956498, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S395'), ('column', 'S316'))] || Time: 0:01:09.699644
size = 950, date: 2023-11-21 11:19:09.667286, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S587'), ('column', 'S645'))] || Time: 0:02:05.983756
size = 630, date: 2023-11-21 11:45:50.236135, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:35.808531
size = 670, date: 2023-11-21 11:47:23.354237, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:41.914625
size = 975, date: 2023-11-21 11:50:33.389166, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S60'), ('column', 'S67'))] || Time: 0:01:08.594534
size = 630, date: 2023-11-21 12:28:08.689562, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:12.984660
size = 670, date: 2023-11-21 12:39:36.157577, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:17.849373
size = 590, date: 2023-11-21 12:41:39.379023, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:29.135543
size = 630, date: 2023-11-21 13:02:00.660589, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:20.662306
size = 670, date: 2023-11-21 13:18:50.249334, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:16.370873
size = 500, date: 2023-11-21 13:19:40.837520, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:16.132622
size = 590, date: 2023-11-21 13:20:16.638592, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:12.488336
size = 550, date: 2023-11-21 13:21:34.629183, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:23.134206
size = 630, date: 2023-11-21 13:35:52.699689, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:17.890369
size = 500, date: 2023-11-21 13:46:28.616920, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:06.554154
size = 590, date: 2023-11-21 13:50:06.656274, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:12.343849
size = 550, date: 2023-11-21 13:53:58.766410, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:09.822371
size = 670, date: 2023-11-21 13:58:57.141260, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:14.851036
size = 500, date: 2023-11-21 14:07:45.311100, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:07.012779
size = 630, date: 2023-11-21 14:09:40.257097, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:18.165974
size = 590, date: 2023-11-21 14:19:00.687278, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:12.073055
size = 550, date: 2023-11-21 14:20:26.128156, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:11.062114
size = 500, date: 2023-11-21 14:29:14.902244, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:06.475512
size = 670, date: 2023-11-21 14:37:39.234520, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:13.987160
size = 630, date: 2023-11-21 14:43:34.536858, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:18.962212
size = 550, date: 2023-11-21 14:46:20.393877, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:10.015081
size = 590, date: 2023-11-21 14:49:04.739923, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:22.399107
size = 500, date: 2023-11-21 14:50:41.113496, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:05.907074
size = 500, date: 2023-11-21 15:12:05.462040, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:06.459396
size = 550, date: 2023-11-21 15:12:16.252465, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:10.857030
size = 630, date: 2023-11-21 15:17:40.781933, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S214'), ('column', 'S370'))] || Time: 0:00:18.763257
size = 670, date: 2023-11-21 15:18:01.273532, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:15.743821
size = 590, date: 2023-11-21 15:19:16.718750, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:11.389720
size = 500, date: 2023-11-21 15:33:26.451550, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S452'), ('column', 'S273'))] || Time: 0:00:06.224044
size = 550, date: 2023-11-21 15:39:02.411383, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:10.045326
size = 590, date: 2023-11-21 15:50:51.151902, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S576'), ('column', 'S564'))] || Time: 0:00:12.052044
size = 670, date: 2023-11-21 15:58:59.302797, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S129'), ('column', 'S60'))] || Time: 0:00:15.191864
size = 550, date: 2023-11-21 16:05:48.248910, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S15'), ('column', 'S454'))] || Time: 0:00:11.550698
size = 20, date: 2023-11-21 16:35:50.445632, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.268616
size = 20, date: 2023-11-21 16:36:02.660440, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.189101
size = 20, date: 2023-11-21 16:36:12.728586, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.118702
size = 20, date: 2023-11-21 16:36:22.192023, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.075698
size = 20, date: 2023-11-21 16:36:32.245060, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.074648
size = 20, date: 2023-11-21 16:36:41.856296, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.074980
size = 20, date: 2023-11-21 16:36:51.980892, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S13'), ('column', 'S5'))] || Time: 0:00:00.086571
size = 445, date: 2023-11-21 16:59:11.816066, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:11.253368
size = 385, date: 2023-11-21 17:00:00.503261, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:07.544171
size = 315, date: 2023-11-21 17:01:11.295819, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:04.771158
size = 315, date: 2023-11-21 17:12:59.917210, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:05.584122
size = 385, date: 2023-11-21 17:15:48.419328, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:03.496705
size = 445, date: 2023-11-21 17:19:55.671856, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:06.201066
size = 250, date: 2023-11-21 17:20:22.825866, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:02.655238
size = 315, date: 2023-11-21 17:22:14.283613, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:02.682017
size = 250, date: 2023-11-21 17:26:52.096365, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:01.668490
size = 385, date: 2023-11-21 17:27:51.679587, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:03.414698
size = 315, date: 2023-11-21 17:31:34.978822, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:02.320369
size = 250, date: 2023-11-21 17:32:09.160297, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:02.813439
size = 445, date: 2023-11-21 17:36:49.536375, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:05.959988
size = 250, date: 2023-11-21 17:37:29.792824, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:01.735602
size = 385, date: 2023-11-21 17:40:14.841497, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:03.085341
size = 315, date: 2023-11-21 17:40:46.167315, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:05.736214
size = 250, date: 2023-11-21 17:42:52.040178, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:00.736305
size = 315, date: 2023-11-21 17:49:50.981792, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:05.846864
size = 385, date: 2023-11-21 17:52:40.371101, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:03.677880
size = 445, date: 2023-11-21 17:53:25.846149, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:05.457671
size = 315, date: 2023-11-21 17:59:10.695591, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S91'), ('column', 'S129'))] || Time: 0:00:07.327931
size = 385, date: 2023-11-21 18:04:46.095826, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:03.000533
size = 445, date: 2023-11-21 18:10:02.500407, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:05.393151
size = 385, date: 2023-11-21 18:17:10.529347, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S19'), ('column', 'S320'))] || Time: 0:00:02.612280
size = 250, date: 2023-11-21 18:20:54.280602, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:01.865402
size = 445, date: 2023-11-21 18:26:06.555073, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S274'), ('column', 'S13'))] || Time: 0:00:05.509850
size = 250, date: 2023-11-21 18:26:10.533607, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S26'), ('column', 'S165'))] || Time: 0:00:00.705443""" \
+ \
"""
size = 920, date: 2023-11-22 22:29:01.283597, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S376'), ('column', 'S215'))] || Time: 0:00:57.618471
size = 975, date: 2023-11-22 22:53:38.333955, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:00:54.021789
size = 920, date: 2023-11-22 23:51:25.202257, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S376'), ('column', 'S215'))] || Time: 0:01:04.332915
size = 975, date: 2023-11-23 00:29:39.893335, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:01:29.250583
size = 920, date: 2023-11-23 01:18:50.863347, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S376'), ('column', 'S215'))] || Time: 0:00:57.101857
size = 975, date: 2023-11-23 02:12:42.589422, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:01:22.510179
size = 920, date: 2023-11-23 02:54:30.694050, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S376'), ('column', 'S215'))] || Time: 0:00:17.162087
size = 975, date: 2023-11-23 03:57:56.486052, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:01:21.225831
size = 975, date: 2023-11-23 05:27:47.720900, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:01:19.528786
size = 975, date: 2023-11-23 06:57:02.908337, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S714'), ('column', 'S335'))] || Time: 0:01:20.684903""" \
+ \
"""
size = 740, date: 2023-11-24 01:28:19.618082, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:29.517366
size = 740, date: 2023-11-24 02:23:07.408577, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:28.999845
size = 740, date: 2023-11-24 03:18:56.138035, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:29.345782
size = 740, date: 2023-11-24 04:16:21.301497, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:29.598085
size = 740, date: 2023-11-24 05:14:50.114677, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:25.477152
size = 740, date: 2023-11-24 06:10:53.366062, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S416'), ('column', 'S138'))] || Time: 0:00:27.773408""" 


# remove the lines with "Original" in them
# string_summer_2023 = '\n'.join([line for line in string_summer_2023.split('\n') if 'Original' not in line])


# plot_averages(string_summer_2023, quadratic=False)
# plot_averages(string_summer_2023, quadratic=False)
# plot_averages(string_summer_2023, quadratic=True)
# without original solution, and upto date: 2023-11-21 17:22:14.283613, R^2 was 0.91
plot_medians(string_summer_2023)