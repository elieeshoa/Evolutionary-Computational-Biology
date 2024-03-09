from datetime import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint



# For each size, we run 5 iterations of the algorithm, and we take the average time
# of the 5 iterations. We then plot the average time against the size of the game.
# We also plot the theoretical time complexity of the algorithm, which is O(n^2).
# We see that the algorithm is indeed O(n^2), as the plot is a straight line.

def get_averages(string):
    sizes = [5, 10, 20, 30, 50, 100, 250, 400, 700, 1000]
    lines = string.split('\n')

    # Create a dictionary of times for every size.
    times = {}
    
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

    # remove previous graph
    plt.clf()

    d = get_all_times(string)
    sizes = [int(size) for size in d.keys()]
    print('sizes are', sizes)
    sizes.sort()
    averages = [np.mean(d[str(size)]) for size in sizes]

    averages_dict, occurences_dict = get_averages(string)
    

    # for size, average in zip(sizes, averages):
    #     print(f"size {size} has square {size**2} and average {average} over {occurences_dict[str(size)]} occurences")
    
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
        plt.xlabel('Size of the game ($n^2$)', fontsize=12)

        plt.ylabel('CPU runtime (seconds)', fontsize=12)
        plt.title('Average CPU Runtime for Finding Unique Solutions')
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
        # plt.ylim(bottom=-18)
       

        # plt.plot(x, fit_fn(x), '--k', label='Linear fit', linewidth=1, color='dodgerblue', alpha=0.8)
        plt.plot(x, fit_fn(x), '--k', label='Linear fit', linewidth=1.2, color='#00A6B7', alpha=0.8)
        # plt.plot(sizes, averages, label='Average Time', linewidth=1, color='red')
        # make the color FF6805
        # plt.plot(sizes, averages, label='Average Time', linewidth=1.5, color='#FF7F00')
        # plt.scatter(sizes, averages, color='red', s=15)
        # plot standard deviation as error bars
        print('sizes are', sizes)
        
        # plt.title('Average CPU Runtime for Finding Unique Solutions')

        # for size, average in zip(sizes, averages):
        #     plt.text(size, average, occurences_dict[(str(int(np.sqrt(size))))], fontsize=10, color='black', alpha=0.8, zorder=4)

        
            
        # calculate R^2
        r2 = np.corrcoef(sizes, averages)[0, 1]**2
        # write r2 on the plot in the center
        plt.text(0.5, 0.85, f"$R^2 = {r2:.2f}$", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, fontsize=15)
        print("R^2:", r2)

        # scatter the averages
        plt.scatter(sizes, averages, color='#FF7F00', s=25, zorder=3)

        plt.errorbar(sizes, averages, yerr=[np.std(d[str(int(np.sqrt(size)))]) for size in sizes], fmt='o', color='r', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
        plt.xlabel('Size of the game ($n^2$)', fontsize=18)
        plt.ylabel('CPU runtime (seconds)', fontsize=18)

        # make the y numbers integers and not 1.0e7
        plt.ticklabel_format(style='plain')
        plt.ylim(bottom=-18)
        # add thin gridlines
        plt.grid(linewidth=0.4)

        # set size of plot to be 7 by 6
        plt.gcf().set_size_inches(7/1, 6/1)
        plt.tight_layout()
        

        plt.savefig('large_games_linear.png', dpi=300)

        

def plot_medians(string, display_text=False):
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

    pprint(times)
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

    if display_text:
        for size, med in zip(sizes, medians):
            plt.text(size, med, occurences[str(int(np.sqrt(size)))], fontsize=10, color='black', alpha=0.8, zorder=4)
            plt.text(size, med-10, int(np.sqrt(size)), fontsize=10, color='black', alpha=0.8, zorder=4, rotation=45, horizontalalignment='center', verticalalignment='center')
    


    # scatter the medians
    plt.scatter(sizes, medians, color='red', s=15)
    # plot the error bars
    # plt.errorbar(sizes, medians, yerr=[medians-q1s, q3s-medians], fmt='o', color='#FF7F00', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.errorbar(sizes, medians, yerr=[medians-q1s, q3s-medians], fmt='o', color='r', alpha=0.8, capsize=3, capthick=1, elinewidth=1, zorder=3)
    plt.xlabel('Size of the game ($n^2$)', fontsize=18)
    plt.ylabel('CPU runtime (seconds)', fontsize=18)
    # plt.title('Median CPU Runtime for Finding Unique Solutions', fontsize=12, pad=10)
    # make the y numbers integers and not 1.0e7
    plt.ticklabel_format(style='plain')
    plt.ylim(bottom=-18)
    # add thin gridlines
    plt.grid(linewidth=0.4)

    # set size of plot to be 7 by 6
    plt.gcf().set_size_inches(7/1, 6/1)
    plt.tight_layout()
    # plt.legend(loc='upper left', fontsize=12)

    plt.savefig('large_games_median.png', dpi=300)

    # plt.show()
    # plt.close()




string_january_2024 = \
"""size = 315, date: 2024-01-08 17:32:35.425971, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:06.706547
size = 385, date: 2024-01-08 17:34:41.159164, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:09.707680
size = 20, date: 2024-01-08 17:35:58.943128, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.060522
size = 20, date: 2024-01-08 17:36:01.375171, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.028552
size = 20, date: 2024-01-08 17:36:03.309950, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.024578
size = 20, date: 2024-01-08 17:36:05.298478, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.060327
size = 20, date: 2024-01-08 17:36:07.276189, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.026634
size = 20, date: 2024-01-08 17:36:09.221535, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.026036
size = 20, date: 2024-01-08 17:36:11.208361, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.025396
size = 20, date: 2024-01-08 17:36:13.087477, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.027711
size = 20, date: 2024-01-08 17:36:15.060302, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.023159
size = 20, date: 2024-01-08 17:36:17.022758, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.026518
size = 20, date: 2024-01-08 17:36:19.018239, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.023984
size = 20, date: 2024-01-08 17:36:20.969785, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S14'), ('column', 'S10'))] || Time: 0:00:00.026651
size = 445, date: 2024-01-08 17:37:14.702683, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:15.603963
size = 500, date: 2024-01-08 17:40:05.221058, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:22.448880
size = 315, date: 2024-01-08 17:43:07.933075, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:02.932770
size = 385, date: 2024-01-08 17:49:02.285635, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.042440
size = 315, date: 2024-01-08 17:51:12.043666, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:02.909902
size = 550, date: 2024-01-08 17:52:42.289523, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S265'), ('column', 'S123'))] || Time: 0:00:33.653375
size = 445, date: 2024-01-08 17:58:12.902565, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:14.459509
size = 315, date: 2024-01-08 17:59:38.259758, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.246974
size = 385, date: 2024-01-08 18:00:48.973099, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:09.992490
size = 500, date: 2024-01-08 18:04:57.605788, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:12.320781
size = 315, date: 2024-01-08 18:07:50.004568, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.832607
size = 385, date: 2024-01-08 18:12:35.354437, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:04.713171
size = 445, date: 2024-01-08 18:14:45.815745, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:06.898345
size = 315, date: 2024-01-08 18:16:19.485626, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:04.978142
size = 550, date: 2024-01-08 18:23:34.851078, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S265'), ('column', 'S123'))] || Time: 0:00:14.083012
size = 385, date: 2024-01-08 18:24:39.085828, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:08.378625
size = 315, date: 2024-01-08 18:24:48.256132, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.768201
size = 500, date: 2024-01-08 18:25:04.665405, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:12.462405
size = 445, date: 2024-01-08 18:31:06.078659, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:08.628038
size = 315, date: 2024-01-08 18:33:17.413955, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.578038
size = 385, date: 2024-01-08 18:36:51.337988, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.524979
size = 315, date: 2024-01-08 18:41:38.347281, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.385838
size = 500, date: 2024-01-08 18:44:52.413599, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.345870
size = 550, date: 2024-01-08 18:47:57.323541, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S265'), ('column', 'S123'))] || Time: 0:00:14.471733
size = 445, date: 2024-01-08 18:48:06.178953, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:09.212628
size = 385, date: 2024-01-08 18:48:57.872782, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.551635
size = 315, date: 2024-01-08 18:50:15.033019, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.447088
size = 315, date: 2024-01-08 18:58:40.165457, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.307659
size = 385, date: 2024-01-08 19:01:02.831586, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.433269
size = 445, date: 2024-01-08 19:05:00.879795, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:08.968120
size = 500, date: 2024-01-08 19:06:00.765754, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.639333
size = 315, date: 2024-01-08 19:07:01.765444, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S213'), ('column', 'S54'))] || Time: 0:00:03.314418
size = 385, date: 2024-01-08 19:13:16.888455, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:13.569503
size = 445, date: 2024-01-08 19:21:28.148599, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:07.984414
size = 385, date: 2024-01-08 19:25:46.215453, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.534671
size = 500, date: 2024-01-08 19:25:57.780376, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.653362
size = 445, date: 2024-01-08 19:36:59.708187, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:08.458189
size = 385, date: 2024-01-08 19:37:23.105250, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:05.899590
size = 500, date: 2024-01-08 19:46:20.145902, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.589333
size = 385, date: 2024-01-08 19:49:09.228452, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S278'), ('column', 'S308'))] || Time: 0:00:06.127056
size = 445, date: 2024-01-08 19:53:17.896527, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:08.981697
size = 500, date: 2024-01-08 20:05:53.641669, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.776652
size = 550, date: 2024-01-08 20:08:18.367234, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:29.380926
size = 445, date: 2024-01-08 20:08:56.501419, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:08.446014
size = 590, date: 2024-01-08 20:09:54.658171, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:36.478635
size = 630, date: 2024-01-08 20:12:08.993551, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:48.103014
size = 445, date: 2024-01-08 20:26:20.563543, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:09.305761
size = 500, date: 2024-01-08 20:26:45.205626, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:12.508281
size = 550, date: 2024-01-08 20:40:27.058180, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:16.290482
size = 445, date: 2024-01-08 20:43:28.825885, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S100'), ('column', 'S430'))] || Time: 0:00:09.128182
size = 590, date: 2024-01-08 20:45:13.802603, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:15.873635
size = 500, date: 2024-01-08 20:47:19.875035, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.411091
size = 630, date: 2024-01-08 20:52:50.922001, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:19.013558
size = 550, date: 2024-01-08 21:05:19.960785, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:15.528126
size = 500, date: 2024-01-08 21:07:31.510067, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:12.176220
size = 590, date: 2024-01-08 21:13:12.341818, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:15.286172
size = 630, date: 2024-01-08 21:23:45.326082, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:18.347655
size = 500, date: 2024-01-08 21:26:12.754034, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S249'), ('column', 'S395'))] || Time: 0:00:11.523627
size = 550, date: 2024-01-08 21:29:17.159297, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:14.763383
size = 250, date: 2024-01-08 21:31:44.039778, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:03.371686
size = 250, date: 2024-01-08 21:37:47.650461, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:02.011051
size = 590, date: 2024-01-08 21:39:32.966336, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:15.177081
size = 250, date: 2024-01-08 21:42:39.210941, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.730249
size = 250, date: 2024-01-08 21:47:40.290638, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.487794
size = 670, date: 2024-01-08 21:52:07.349530, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S574'), ('column', 'S114'))] || Time: 0:00:58.827056
size = 250, date: 2024-01-08 21:52:43.969646, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.597629
size = 550, date: 2024-01-08 21:54:25.179520, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:14.610844
size = 630, date: 2024-01-08 21:55:09.816755, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:16.307407
size = 250, date: 2024-01-08 21:57:47.122644, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.494297
size = 250, date: 2024-01-08 22:02:45.043684, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.551046
size = 590, date: 2024-01-08 22:07:43.127810, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:19.734510
size = 250, date: 2024-01-08 22:07:47.663929, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.426561
size = 250, date: 2024-01-08 22:12:50.135664, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:01.717876
size = 250, date: 2024-01-08 22:17:39.963146, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:02.006230
size = 550, date: 2024-01-08 22:19:30.961299, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:14.889261
size = 250, date: 2024-01-08 22:22:41.563620, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:02.027338
size = 630, date: 2024-01-08 22:26:28.125750, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:12.260842
size = 250, date: 2024-01-08 22:27:42.372582, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S94'), ('column', 'S100'))] || Time: 0:00:02.452737
size = 590, date: 2024-01-08 22:35:26.801625, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:19.459644
size = 670, date: 2024-01-08 22:39:02.116121, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S574'), ('column', 'S114'))] || Time: 0:00:27.419109
size = 550, date: 2024-01-08 22:43:46.089761, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:15.163104
size = 630, date: 2024-01-08 22:57:34.465451, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S432'), ('column', 'S331'))] || Time: 0:00:21.901735
size = 590, date: 2024-01-08 23:04:06.331249, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:18.539993
size = 705, date: 2024-01-08 23:06:43.311214, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:01:41.897332
size = 550, date: 2024-01-08 23:08:26.604235, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:15.024376
size = 590, date: 2024-01-08 23:31:27.198680, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:18.927397
size = 550, date: 2024-01-08 23:33:01.131571, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:14.551615
size = 630, date: 2024-01-08 23:50:55.285671, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:01:17.123281
size = 670, date: 2024-01-08 23:51:48.912328, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:01:20.404055
size = 550, date: 2024-01-08 23:57:31.069239, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:15.499621
size = 590, date: 2024-01-08 23:57:46.237929, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:18.981028
size = 550, date: 2024-01-09 00:21:03.698105, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:11.525186
size = 705, date: 2024-01-09 00:22:07.299690, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:41.586054
size = 590, date: 2024-01-09 00:23:54.014502, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:17.712515
size = 630, date: 2024-01-09 00:43:54.817084, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:28.051844
size = 550, date: 2024-01-09 00:45:36.655966, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S59'), ('column', 'S447'))] || Time: 0:00:14.846629
size = 590, date: 2024-01-09 00:51:32.571418, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:17.512144
size = 670, date: 2024-01-09 01:01:04.233577, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:38.206861
size = 590, date: 2024-01-09 01:17:49.933847, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S581'), ('column', 'S7'))] || Time: 0:00:17.659243
size = 705, date: 2024-01-09 01:31:48.906622, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:44.473255
size = 630, date: 2024-01-09 01:32:45.924597, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:26.315685
size = 740, date: 2024-01-09 01:46:36.199430, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:02:11.388739
size = 670, date: 2024-01-09 02:04:13.478189, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:42.268539
size = 630, date: 2024-01-09 02:18:05.115646, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:20.991004
size = 705, date: 2024-01-09 02:42:41.320625, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:31.805895
size = 1000, date: 2024-01-09 02:46:14.062248, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:04:29.193188
size = 670, date: 2024-01-09 03:03:45.479707, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:38.760375
size = 630, date: 2024-01-09 03:04:42.375577, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:26.706504
size = 740, date: 2024-01-09 03:14:10.247625, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:31.454534
size = 630, date: 2024-01-09 03:46:12.086278, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:17.382706
size = 670, date: 2024-01-09 04:00:53.311365, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:39.232598
size = 705, date: 2024-01-09 04:01:38.251172, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:36.926525
size = 740, date: 2024-01-09 04:25:48.721515, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:01:09.372173
size = 630, date: 2024-01-09 04:29:35.062774, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:21.963590
size = 670, date: 2024-01-09 04:58:18.479614, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:30.668782
size = 705, date: 2024-01-09 05:02:44.498343, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:35.121113
size = 630, date: 2024-01-09 05:14:05.099488, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:25.010922
size = 1000, date: 2024-01-09 05:21:13.878494, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:02:28.956935
size = 740, date: 2024-01-09 05:34:36.458271, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:01:10.952830
size = 670, date: 2024-01-09 05:57:00.642927, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:37.698436
size = 630, date: 2024-01-09 05:58:28.017490, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:19.040882
size = 705, date: 2024-01-09 06:06:17.126113, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:34.299418
size = 630, date: 2024-01-09 06:33:37.001916, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:14.459739
size = 740, date: 2024-01-09 06:39:39.857468, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:45.560789
size = 670, date: 2024-01-09 06:47:12.352429, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:29.917138
size = 705, date: 2024-01-09 07:05:44.547975, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:27.462255
size = 630, date: 2024-01-09 07:06:08.032791, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:14.825862
size = 1000, date: 2024-01-09 07:19:52.428563, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:02:31.977807
size = 670, date: 2024-01-09 07:32:10.176619, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:30.977862
size = 740, date: 2024-01-09 07:36:13.712397, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:42.624193
size = 630, date: 2024-01-09 07:39:28.957304, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S152'), ('column', 'S598'))] || Time: 0:00:13.901954
size = 705, date: 2024-01-09 07:55:40.881787, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:26.225707
size = 670, date: 2024-01-09 08:17:08.369905, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:25.640478
size = 740, date: 2024-01-09 08:29:27.603720, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:52.862194
size = 705, date: 2024-01-09 08:43:46.718615, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:27.661468
size = 670, date: 2024-01-09 09:01:07.961359, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:29.019561
size = 1000, date: 2024-01-09 09:05:11.251963, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:03:01.772232
size = 740, date: 2024-01-09 09:26:12.173348, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:44.312189
size = 705, date: 2024-01-09 09:34:06.520102, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:25.939037
size = 670, date: 2024-01-09 09:49:18.222416, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S180'), ('column', 'S559'))] || Time: 0:00:29.437707
size = 740, date: 2024-01-09 10:21:44.582323, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:00:40.153329
size = 705, date: 2024-01-09 10:24:04.160952, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S400'), ('column', 'S21'))] || Time: 0:00:24.897065
size = 1000, date: 2024-01-09 10:57:01.164206, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:04:14.310231
size = 740, date: 2024-01-09 11:35:04.626187, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:01:14.391450
size = 740, date: 2024-01-09 13:00:53.282315, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:01:08.302911
size = 1000, date: 2024-01-09 13:47:19.031629, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:02:42.733594
size = 740, date: 2024-01-09 14:32:29.402815, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S562'), ('column', 'S63'))] || Time: 0:03:25.872533
size = 975, date: 2024-01-09 14:52:08.887054, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:03:23.531864
size = 950, date: 2024-01-09 14:54:07.033658, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:03:05.207473
size = 920, date: 2024-01-09 15:21:50.915259, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:04:55.529414
size = 895, date: 2024-01-09 16:08:34.661445, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S289'), ('column', 'S847'))] || Time: 0:04:02.432458
size = 1000, date: 2024-01-09 16:22:12.806196, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:03:51.269050
size = 950, date: 2024-01-09 16:31:21.513665, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:16.635525
size = 975, date: 2024-01-09 16:36:34.876600, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:00:59.720394
size = 950, date: 2024-01-09 17:46:33.860476, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:33.730560
size = 920, date: 2024-01-09 17:54:31.144070, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:01:04.526773
size = 975, date: 2024-01-09 17:59:54.694256, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:00:43.046337
size = 895, date: 2024-01-09 18:24:09.403065, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S289'), ('column', 'S847'))] || Time: 0:01:43.693224
size = 1000, date: 2024-01-09 18:47:47.392953, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:03:28.991797
size = 950, date: 2024-01-09 19:04:33.353961, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:06.093426
size = 975, date: 2024-01-09 19:22:27.385200, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:31.506884
size = 920, date: 2024-01-09 19:58:20.558691, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:31.731912
size = 895, date: 2024-01-09 20:20:24.489166, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S289'), ('column', 'S847'))] || Time: 0:02:15.117804
size = 950, date: 2024-01-09 20:23:10.367942, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:30.568441
size = 975, date: 2024-01-09 20:47:54.744814, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:32.041330
size = 1000, date: 2024-01-09 21:22:43.082681, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S163'), ('column', 'S928'))] || Time: 0:03:29.774575
size = 950, date: 2024-01-09 21:38:34.895553, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:04.631413
size = 920, date: 2024-01-09 22:01:14.343551, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:27.286472
size = 975, date: 2024-01-09 22:13:19.701082, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:27.989449
size = 895, date: 2024-01-09 22:21:11.872010, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S289'), ('column', 'S847'))] || Time: 0:01:52.210682
size = 950, date: 2024-01-09 22:56:58.784450, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:27.094520
size = 975, date: 2024-01-09 23:38:23.250490, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:28.514240
size = 920, date: 2024-01-10 00:06:34.703655, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:31.204213
size = 950, date: 2024-01-10 00:13:02.970741, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:31.727380
size = 895, date: 2024-01-10 00:13:44.014166, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S289'), ('column', 'S847'))] || Time: 0:00:16.470964
size = 975, date: 2024-01-10 01:02:49.610710, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:00:58.782053
size = 950, date: 2024-01-10 01:28:57.346273, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:17.685677
size = 920, date: 2024-01-10 02:12:10.652125, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:01:29.080630
size = 975, date: 2024-01-10 02:27:43.281068, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:26.752814
size = 950, date: 2024-01-10 02:45:12.419942, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:09.885887
size = 775, date: 2024-01-10 03:11:25.094431, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:02:51.760617
size = 975, date: 2024-01-10 03:53:02.430478, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:27.150285
size = 950, date: 2024-01-10 04:04:06.099064, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:26.369452
size = 920, date: 2024-01-10 04:15:55.357412, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:13.759648
size = 775, date: 2024-01-10 05:15:41.337348, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:17.615488
size = 975, date: 2024-01-10 05:18:09.528419, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:25.716011
size = 950, date: 2024-01-10 05:20:16.346001, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S805'), ('column', 'S783'))] || Time: 0:01:02.453292
size = 920, date: 2024-01-10 06:11:01.118616, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:01:54.488535
size = 975, date: 2024-01-10 06:43:23.995194, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S126'), ('column', 'S171'))] || Time: 0:01:26.652972
size = 775, date: 2024-01-10 07:01:31.648316, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:17.829196
size = 805, date: 2024-01-10 07:21:16.355801, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:02:49.525137
size = 920, date: 2024-01-10 08:05:41.131141, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:31.246305
size = 775, date: 2024-01-10 08:50:56.047287, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:27.017735
size = 805, date: 2024-01-10 09:10:57.700686, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:13.820825
size = 920, date: 2024-01-10 10:09:39.170992, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:29.499635
size = 775, date: 2024-01-10 10:38:21.074756, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:17.018243
size = 805, date: 2024-01-10 10:45:31.644997, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:35.830973
size = 805, date: 2024-01-10 12:03:38.470713, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:12.118538
size = 920, date: 2024-01-10 12:05:12.402732, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:02:32.054896
size = 775, date: 2024-01-10 12:26:33.547978, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:44.492001
size = 805, date: 2024-01-10 13:14:33.594307, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:00:51.528099
size = 920, date: 2024-01-10 13:35:22.693176, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S218'), ('column', 'S558'))] || Time: 0:01:09.351035
size = 835, date: 2024-01-10 13:39:35.643313, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S488'), ('column', 'S810'))] || Time: 0:02:54.862286
size = 775, date: 2024-01-10 14:10:45.020202, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:57.221713
size = 805, date: 2024-01-10 14:22:54.676266, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:08.845961
size = 865, date: 2024-01-10 15:01:40.404779, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:02:15.689492
size = 835, date: 2024-01-10 15:07:07.491540, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S488'), ('column', 'S810'))] || Time: 0:00:28.531580
size = 805, date: 2024-01-10 15:32:33.808877, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:06.082906
size = 775, date: 2024-01-10 15:59:33.880869, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:58.110847
size = 835, date: 2024-01-10 16:16:49.669083, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S488'), ('column', 'S810'))] || Time: 0:01:06.962881
size = 865, date: 2024-01-10 16:22:33.778661, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:00:53.818852
size = 805, date: 2024-01-10 16:40:26.041016, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:01.608070
size = 865, date: 2024-01-10 17:27:08.468775, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:09.411668
size = 835, date: 2024-01-10 17:44:00.543699, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S488'), ('column', 'S810'))] || Time: 0:01:00.598864
size = 775, date: 2024-01-10 17:48:09.582966, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:02:16.914589
size = 805, date: 2024-01-10 18:02:13.439307, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:00:58.718450
size = 865, date: 2024-01-10 18:30:21.591202, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:05.435923
size = 835, date: 2024-01-10 18:57:51.151060, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S488'), ('column', 'S810'))] || Time: 0:00:18.906754
size = 805, date: 2024-01-10 19:12:00.538312, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:20.558287
size = 865, date: 2024-01-10 19:33:57.378699, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:16.200672
size = 775, date: 2024-01-10 19:55:36.283885, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:49.689829
size = 805, date: 2024-01-10 20:24:17.020307, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:03.128887
size = 865, date: 2024-01-10 20:37:24.339572, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:05.808356
size = 865, date: 2024-01-10 21:42:40.343135, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:14.805699
size = 805, date: 2024-01-10 21:47:01.505261, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S452'), ('column', 'S708'))] || Time: 0:01:14.921923
size = 775, date: 2024-01-10 21:55:34.954524, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:02:01.520044
size = 865, date: 2024-01-10 22:49:51.249158, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:11.252177
size = 775, date: 2024-01-10 23:27:50.449589, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S54'), ('column', 'S540'))] || Time: 0:01:34.811397
size = 865, date: 2024-01-10 23:59:18.405452, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:06.826108
size = 865, date: 2024-01-11 01:05:24.385285, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:09.502685
size = 865, date: 2024-01-11 02:10:16.444179, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:01:08.613016
size = 865, date: 2024-01-11 03:12:25.325438, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S547'), ('column', 'S791'))] || Time: 0:00:41.862925

# on Jan 13 starting 2:50pm. This is doing 5 iterations at a time
size = 920, date: 2024-01-13 16:50:35.829957, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:02:46.025240
size = 950, date: 2024-01-13 16:59:35.700094, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:04:15.810210
size = 1000, date: 2024-01-13 16:59:37.778701, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:04:14.453289
size = 895, date: 2024-01-13 17:14:05.001386, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:03:43.734650
size = 975, date: 2024-01-13 17:28:14.533076, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:04:19.255488
size = 920, date: 2024-01-13 18:31:14.169415, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:00:41.026778
size = 1000, date: 2024-01-13 18:56:03.536590, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:01:12.544320
size = 950, date: 2024-01-13 18:57:06.732924, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:01:19.633157
size = 895, date: 2024-01-13 19:46:30.731573, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:01:32.063871
size = 920, date: 2024-01-13 19:51:15.611797, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:00:44.581476
size = 1000, date: 2024-01-13 20:29:27.827516, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:01:16.521217
size = 975, date: 2024-01-13 20:31:51.627095, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:02:13.579116
size = 950, date: 2024-01-13 20:34:34.705984, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:01:33.482032
size = 920, date: 2024-01-13 21:15:18.407602, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:00:57.336121
size = 895, date: 2024-01-13 21:27:36.620948, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:01:12.479798
size = 1000, date: 2024-01-13 22:01:13.393427, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:01:13.195515
size = 975, date: 2024-01-13 22:05:45.400017, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:00:35.113307
size = 950, date: 2024-01-13 22:07:26.712479, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:01:27.690617
size = 895, date: 2024-01-13 22:38:02.960686, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:01:14.710081
size = 920, date: 2024-01-13 22:38:49.394545, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:00:43.806246
size = 975, date: 2024-01-13 23:31:25.820261, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:01:22.867332
size = 1000, date: 2024-01-13 23:38:48.641635, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:01:13.940139
size = 950, date: 2024-01-13 23:43:25.873658, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:01:34.283593
size = 895, date: 2024-01-13 23:47:43.523610, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:01:08.076717
size = 920, date: 2024-01-14 00:01:31.652642, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S541'), ('column', 'S816'))] || Time: 0:01:09.630661
size = 895, date: 2024-01-14 00:56:02.954694, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S545'), ('column', 'S619'))] || Time: 0:01:22.451521
size = 975, date: 2024-01-14 00:56:51.296003, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:01:39.437151
size = 1000, date: 2024-01-14 01:15:44.701009, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S16'), ('column', 'S477'))] || Time: 0:01:15.432131
size = 950, date: 2024-01-14 01:21:29.949570, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S441'), ('column', 'S262'))] || Time: 0:01:03.000638
size = 975, date: 2024-01-14 02:22:27.267061, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S574'), ('column', 'S546'))] || Time: 0:01:30.318088
size = 895, date: 2024-01-14 02:40:38.987843, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:02:21.828617
size = 920, date: 2024-01-14 02:47:14.339349, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:02:44.130557
size = 1000, date: 2024-01-14 02:56:20.043751, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:03:41.819627
size = 950, date: 2024-01-14 03:00:37.276256, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:04:06.892220
size = 975, date: 2024-01-14 03:26:34.227529, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:03:18.749595
size = 895, date: 2024-01-14 04:05:30.419619, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:01:02.176233
size = 920, date: 2024-01-14 04:18:41.874018, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:00:48.189057
size = 1000, date: 2024-01-14 04:44:05.357711, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:01:25.169701
size = 950, date: 2024-01-14 05:03:13.682431, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:01:18.773692
size = 975, date: 2024-01-14 05:07:46.947711, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:00:54.297657
size = 895, date: 2024-01-14 05:13:04.994778, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:01:16.252475
size = 920, date: 2024-01-14 05:31:29.863458, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:01:03.027654
size = 1000, date: 2024-01-14 06:09:58.921009, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:01:30.260276
size = 895, date: 2024-01-14 06:21:16.807791, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:01:13.464079
size = 975, date: 2024-01-14 06:28:11.440292, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:01:05.031811
size = 950, date: 2024-01-14 06:42:12.921742, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:01:58.218358
size = 920, date: 2024-01-14 06:43:52.426000, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:01:00.545721
size = 895, date: 2024-01-14 07:30:13.972324, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:01:10.471744
size = 1000, date: 2024-01-14 07:34:17.565981, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:01:29.260886
size = 975, date: 2024-01-14 07:49:43.886896, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:01:09.108623
size = 920, date: 2024-01-14 07:54:01.426698, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:01:04.539781
size = 950, date: 2024-01-14 08:22:31.550938, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:01:56.773619
size = 895, date: 2024-01-14 08:36:21.412766, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S716'), ('column', 'S551'))] || Time: 0:01:00.263519
size = 1000, date: 2024-01-14 09:00:47.271518, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:00:55.518124
size = 920, date: 2024-01-14 09:06:27.104173, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S464'), ('column', 'S449'))] || Time: 0:00:53.783915
size = 975, date: 2024-01-14 09:08:32.238517, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:01:21.209728
size = 950, date: 2024-01-14 09:59:28.953957, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:01:48.733285
size = 1000, date: 2024-01-14 10:24:22.297808, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S737'), ('column', 'S477'))] || Time: 0:01:18.638808
size = 975, date: 2024-01-14 10:27:42.873172, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S838'), ('column', 'S33'))] || Time: 0:01:16.425239
size = 950, date: 2024-01-14 11:37:58.385323, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S246'), ('column', 'S277'))] || Time: 0:01:35.945733
size = 775, date: 2024-01-14 12:17:02.604001, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:01:29.483144
size = 835, date: 2024-01-14 12:22:49.194168, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:01:52.626242
size = 865, date: 2024-01-14 12:26:50.255211, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:02:35.842510
size = 805, date: 2024-01-14 12:27:33.374107, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:02:42.417056
size = 1000, date: 2024-01-14 12:38:09.154233, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:03:31.506146
size = 775, date: 2024-01-14 13:26:07.080575, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:00:48.033980
size = 835, date: 2024-01-14 13:41:58.476145, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:01:08.171558
size = 805, date: 2024-01-14 13:57:52.709398, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:01:12.398265
size = 865, date: 2024-01-14 13:58:13.478898, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:01:00.546452
size = 775, date: 2024-01-14 14:28:28.797286, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:00:55.330860
size = 1000, date: 2024-01-14 14:35:05.569426, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:44.778124
size = 835, date: 2024-01-14 14:49:41.721079, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:01:08.530778
size = 805, date: 2024-01-14 15:11:46.185154, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:01:20.013484
size = 865, date: 2024-01-14 15:16:59.251675, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:00:59.457231
size = 775, date: 2024-01-14 15:32:19.547878, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:00:49.345367
size = 835, date: 2024-01-14 15:55:16.600292, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:00:30.297877
size = 1000, date: 2024-01-14 16:06:30.015218, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:49.732119
size = 805, date: 2024-01-14 16:27:07.778075, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:01:23.707101
size = 775, date: 2024-01-14 16:34:06.284042, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:01:00.079943
size = 865, date: 2024-01-14 16:35:35.582150, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:01:17.047320
size = 835, date: 2024-01-14 17:04:12.695333, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:01:10.230385
size = 775, date: 2024-01-14 17:36:35.784675, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S404'), ('column', 'S451'))] || Time: 0:00:49.500722
size = 1000, date: 2024-01-14 17:42:06.938663, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:37.266955
size = 805, date: 2024-01-14 17:44:20.120455, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:00:53.417779
size = 865, date: 2024-01-14 17:51:30.597980, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:01:02.884268
size = 835, date: 2024-01-14 18:12:24.413940, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S289'), ('column', 'S633'))] || Time: 0:01:12.160096
size = 775, date: 2024-01-14 18:56:14.240928, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:02:00.826675
size = 805, date: 2024-01-14 19:00:15.373736, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S101'), ('column', 'S48'))] || Time: 0:01:19.326408
size = 865, date: 2024-01-14 19:10:09.616222, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S410'), ('column', 'S10'))] || Time: 0:00:53.586682
size = 835, date: 2024-01-14 19:10:47.535571, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:02:48.290926
size = 1000, date: 2024-01-14 19:20:20.618240, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:56.757039
size = 775, date: 2024-01-14 20:12:38.388019, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:00:39.712680
size = 835, date: 2024-01-14 20:52:49.982116, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:01:15.985347
size = 1000, date: 2024-01-14 20:54:47.655751, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:37.924833
size = 775, date: 2024-01-14 21:17:37.499682, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:00:44.285930
size = 835, date: 2024-01-14 22:20:34.581444, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:01:53.575367
size = 775, date: 2024-01-14 22:26:50.057697, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:00:42.643420
size = 1000, date: 2024-01-14 22:28:13.959941, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:53.896504
# R^2 = 0.92 up till here
size = 775, date: 2024-01-14 23:31:47.027236, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:00:52.047324
size = 835, date: 2024-01-14 23:46:29.383014, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:00:33.884876
size = 1000, date: 2024-01-15 00:06:26.264879, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:37.809203
size = 775, date: 2024-01-15 00:34:51.986228, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S461'), ('column', 'S467'))] || Time: 0:00:37.354775
size = 835, date: 2024-01-15 01:08:01.008547, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:01:32.265043
size = 1000, date: 2024-01-15 01:40:06.450713, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:42.254093
size = 835, date: 2024-01-15 02:28:53.200456, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S686'), ('column', 'S220'))] || Time: 0:01:22.418443
size = 1000, date: 2024-01-15 03:22:25.138602, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:44.969537
size = 865, date: 2024-01-15 04:18:15.584045, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:03:12.483836
size = 1000, date: 2024-01-15 05:01:39.326546, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S337'), ('column', 'S300'))] || Time: 0:00:43.957785
size = 865, date: 2024-01-15 06:15:05.571110, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:16.151133
size = 865, date: 2024-01-15 07:48:29.400832, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:32.186143
size = 865, date: 2024-01-15 09:22:39.512190, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:31.080727
size = 805, date: 2024-01-15 10:31:43.876554, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:02:34.600832
size = 865, date: 2024-01-15 10:53:10.739727, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:44.747149
size = 775, date: 2024-01-15 11:03:19.697293, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:02:17.505599
size = 805, date: 2024-01-15 12:03:35.934299, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:09.420677
size = 1000, date: 2024-01-15 12:08:05.161177, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:04:25.887365
size = 865, date: 2024-01-15 12:21:20.770623, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:21.460274
size = 775, date: 2024-01-15 12:23:10.156278, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:35.789954
size = 805, date: 2024-01-15 13:19:18.941762, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:26.698222
size = 775, date: 2024-01-15 13:29:00.485173, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:37.088177
size = 865, date: 2024-01-15 13:49:21.914327, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:32.870980
size = 1000, date: 2024-01-15 14:29:51.368895, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:01:20.905477
size = 805, date: 2024-01-15 14:33:14.402569, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:09.148608
size = 775, date: 2024-01-15 14:34:24.028823, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:50.404823
size = 865, date: 2024-01-15 15:18:47.944859, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:14.470753
size = 775, date: 2024-01-15 15:38:51.862098, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:01:09.378441
size = 805, date: 2024-01-15 15:47:31.141425, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:05.187847
size = 740, date: 2024-01-15 16:01:43.760429, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:01:27.019369
size = 1000, date: 2024-01-15 16:22:58.526867, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:01:35.179328
size = 775, date: 2024-01-15 16:46:35.984528, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:36.771058
size = 865, date: 2024-01-15 16:47:15.098380, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S276'), ('column', 'S515'))] || Time: 0:01:27.256686
size = 805, date: 2024-01-15 17:03:32.256142, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:14.131054
size = 740, date: 2024-01-15 17:08:20.779559, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:00:52.484632
size = 775, date: 2024-01-15 17:51:55.174292, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:45.753898
size = 740, date: 2024-01-15 18:01:39.278987, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:00:48.420606
size = 1000, date: 2024-01-15 18:19:32.237026, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:01:41.870212
size = 805, date: 2024-01-15 18:20:45.594172, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:22.556681
size = 740, date: 2024-01-15 18:55:36.179874, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:00:45.927914
size = 775, date: 2024-01-15 18:56:58.319753, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:33.999697
size = 805, date: 2024-01-15 19:40:20.930609, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:23.807409
size = 740, date: 2024-01-15 19:49:11.588638, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:00:57.053584
# R^2 = 0.9498 up till here
size = 775, date: 2024-01-15 20:09:59.670518, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S53'), ('column', 'S288'))] || Time: 0:00:38.408857
size = 1000, date: 2024-01-15 20:17:52.247963, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:02:54.278921
size = 740, date: 2024-01-15 20:43:54.751518, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:01:00.243496
size = 740, date: 2024-01-15 20:52:51.917356, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:01:57.300572
size = 805, date: 2024-01-15 21:01:22.263080, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S707'), ('column', 'S240'))] || Time: 0:01:18.420649
size = 740, date: 2024-01-15 21:38:05.441738, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:01:11.309986
size = 705, date: 2024-01-15 22:18:49.450091, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:01:17.775478
size = 740, date: 2024-01-15 22:18:53.756752, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:00:27.810156
size = 805, date: 2024-01-15 22:25:45.512291, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:02:00.446818
size = 1000, date: 2024-01-15 22:27:08.306430, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:02:10.457428
size = 740, date: 2024-01-15 22:32:31.820936, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:01:22.523408
# R^2 = 0.95119 up till here
size = 740, date: 2024-01-15 22:35:42.092951, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:01:53.738777
size = 740, date: 2024-01-15 23:13:16.068288, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:00:45.500975
size = 705, date: 2024-01-15 23:22:25.820137, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:37.805178
size = 740, date: 2024-01-15 23:24:10.544267, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S192'), ('column', 'S404'))] || Time: 0:01:00.425391
size = 805, date: 2024-01-15 23:55:24.916338, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:41.638998
size = 740, date: 2024-01-16 00:00:11.393956, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:00:33.569345
size = 740, date: 2024-01-16 00:17:16.705724, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:01:05.650086
size = 705, date: 2024-01-16 00:19:40.594739, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:18.031637
size = 1000, date: 2024-01-16 00:23:50.009799, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:02:38.550422
size = 740, date: 2024-01-16 01:07:02.826239, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:01:00.803533
size = 805, date: 2024-01-16 01:11:33.765309, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:55.468635
size = 705, date: 2024-01-16 01:14:52.692807, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:34.617151
size = 740, date: 2024-01-16 01:21:08.452227, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:01:13.163704
size = 705, date: 2024-01-16 02:10:01.853512, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:28.619129
size = 740, date: 2024-01-16 02:14:14.490111, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:00:29.076848
size = 740, date: 2024-01-16 02:23:46.951351, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:00:52.769580
size = 1000, date: 2024-01-16 02:24:49.671936, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:02:27.036659
size = 805, date: 2024-01-16 02:26:48.036058, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:01:06.697608
# R^2 = 0.9528 up till here
size = 705, date: 2024-01-16 03:05:11.021674, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:20.416212
size = 740, date: 2024-01-16 03:21:02.449581, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:00:59.173157
size = 740, date: 2024-01-16 03:26:31.490684, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:00:42.925877
size = 805, date: 2024-01-16 03:41:45.526877, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:46.444251
size = 705, date: 2024-01-16 04:00:33.324693, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:19.885558
size = 1000, date: 2024-01-16 04:25:01.369502, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S286'), ('column', 'S296'))] || Time: 0:02:25.161144
size = 740, date: 2024-01-16 04:30:07.140313, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:01:48.934401
size = 740, date: 2024-01-16 04:32:46.713210, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:01:14.052371
size = 805, date: 2024-01-16 04:58:33.950178, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:38.748360
size = 705, date: 2024-01-16 04:59:28.914574, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:24.999826
size = 740, date: 2024-01-16 05:36:22.304008, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S392'), ('column', 'S288'))] || Time: 0:01:36.885103
size = 740, date: 2024-01-16 05:48:28.534092, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:00:56.044151
size = 705, date: 2024-01-16 05:56:38.923238, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S498'), ('column', 'S35'))] || Time: 0:00:31.368234
size = 805, date: 2024-01-16 06:18:17.544031, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:01:23.297902
size = 740, date: 2024-01-16 07:04:07.662695, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:00:55.217787
size = 805, date: 2024-01-16 07:36:11.675700, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:33.543598
size = 740, date: 2024-01-16 08:21:49.409560, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S365'), ('column', 'S662'))] || Time: 0:01:02.990044
size = 805, date: 2024-01-16 08:56:16.429757, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S84'), ('column', 'S222'))] || Time: 0:00:36.074502
size = 670, date: 2024-01-16 13:13:14.683012, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:01:41.232827
size = 740, date: 2024-01-16 13:21:31.328081, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:02:04.352122
size = 630, date: 2024-01-16 13:28:13.514104, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:01:36.679938
size = 805, date: 2024-01-16 13:28:45.616247, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:02:37.210338
size = 705, date: 2024-01-16 13:41:01.400623, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:02:26.802495
size = 670, date: 2024-01-16 14:18:46.745725, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:27.045002
size = 740, date: 2024-01-16 14:39:50.558296, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:34.681392
size = 630, date: 2024-01-16 14:49:36.050703, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:29.306077
size = 805, date: 2024-01-16 15:01:31.054007, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:00:53.581838
size = 670, date: 2024-01-16 15:09:40.161022, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:37.273440
size = 705, date: 2024-01-16 15:33:44.908968, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:00:29.929346
size = 740, date: 2024-01-16 15:41:12.678182, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:42.751402
size = 630, date: 2024-01-16 15:57:55.637035, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:31.384860
size = 670, date: 2024-01-16 16:00:43.704125, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:36.315469
size = 805, date: 2024-01-16 16:18:28.689857, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:01:08.673636
size = 740, date: 2024-01-16 16:49:27.776005, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:40.662264
# R^2 = 0.95579 up till here
size = 670, date: 2024-01-16 17:02:00.282923, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:37.827149
size = 630, date: 2024-01-16 17:04:26.627671, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:33.027081
size = 705, date: 2024-01-16 17:04:40.534952, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:00:25.413830
size = 805, date: 2024-01-16 17:43:54.849193, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:01:00.934286
size = 670, date: 2024-01-16 18:01:30.778863, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:33.404745
size = 740, date: 2024-01-16 18:01:47.956384, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:39.513428
size = 630, date: 2024-01-16 18:06:47.077154, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:34.005208
size = 705, date: 2024-01-16 18:25:36.876867, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:07.986573
# R^2 = 0.958089 up till here
size = 670, date: 2024-01-16 18:50:34.349338, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:44.351546
size = 805, date: 2024-01-16 18:53:44.325061, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:00:48.269959
size = 740, date: 2024-01-16 19:01:04.784437, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:31.910681
size = 630, date: 2024-01-16 19:15:01.864602, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:42.519987
size = 670, date: 2024-01-16 19:39:28.343468, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:43.382519
# R^2 = 0.96138 up till here
size = 705, date: 2024-01-16 19:48:09.857911, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:06.640850
size = 740, date: 2024-01-16 20:00:30.146771, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:37.734127
size = 805, date: 2024-01-16 20:04:45.597255, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:00:58.204532
size = 630, date: 2024-01-16 20:20:00.491526, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:38.335755
size = 670, date: 2024-01-16 20:28:26.313967, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S80'), ('column', 'S260'))] || Time: 0:00:42.625703
size = 740, date: 2024-01-16 21:06:57.046800, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:01:15.649167
size = 705, date: 2024-01-16 21:12:08.677764, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:14.053341
size = 630, date: 2024-01-16 21:24:31.126329, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:46.539465
size = 805, date: 2024-01-16 21:33:16.078935, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:01:49.443976
# R^2 = 0.96326 up till here
size = 590, date: 2024-01-16 22:05:08.166903, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:45.191503
size = 740, date: 2024-01-16 22:32:57.261231, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S160'), ('column', 'S508'))] || Time: 0:00:37.492914
size = 630, date: 2024-01-16 22:37:26.389973, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S202'), ('column', 'S339'))] || Time: 0:00:48.564049
size = 705, date: 2024-01-16 22:39:11.475158, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:26.954667
size = 590, date: 2024-01-16 22:43:28.907806, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:12.292499
size = 805, date: 2024-01-16 23:03:11.043453, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:01:40.559216
size = 590, date: 2024-01-16 23:22:54.925315, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:21.488972
# R^2 = 0.96426 up till here
size = 590, date: 2024-01-17 00:03:40.548254, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:17.048064
size = 705, date: 2024-01-17 00:06:07.740107, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:25.065344
size = 500, date: 2024-01-17 00:06:38.275897, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:31.271863
size = 805, date: 2024-01-17 00:16:46.288790, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S557'), ('column', 'S338'))] || Time: 0:01:24.763690
size = 445, date: 2024-01-17 00:21:41.946566, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:24.212343
size = 550, date: 2024-01-17 00:22:22.037373, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:43.752314
size = 590, date: 2024-01-17 00:43:51.293178, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:15.012908
size = 500, date: 2024-01-17 00:45:06.896973, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:10.173632
size = 445, date: 2024-01-17 00:50:38.987300, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:08.066434
size = 500, date: 2024-01-17 01:06:28.152207, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:09.611104
size = 445, date: 2024-01-17 01:08:10.038781, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:10.018802
size = 550, date: 2024-01-17 01:12:42.901183, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:35.128752
size = 590, date: 2024-01-17 01:17:19.364279, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:25.640374
size = 445, date: 2024-01-17 01:32:35.594767, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:11.252849
size = 500, date: 2024-01-17 01:36:17.659666, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:11.738129
size = 705, date: 2024-01-17 01:40:15.417989, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S664'), ('column', 'S529'))] || Time: 0:01:19.619619
size = 445, date: 2024-01-17 01:52:29.020465, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:08.428453
size = 590, date: 2024-01-17 01:53:57.341486, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:18.457441
size = 500, date: 2024-01-17 01:57:56.677928, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:08.883524
size = 550, date: 2024-01-17 02:05:32.178503, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:39.097697
size = 445, date: 2024-01-17 02:10:54.551629, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:09.340339
size = 500, date: 2024-01-17 02:21:27.640632, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:11.200575
size = 590, date: 2024-01-17 02:26:58.673248, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:22.802642
size = 445, date: 2024-01-17 02:32:02.074632, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:11.436166
size = 550, date: 2024-01-17 02:38:34.877896, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:18.217620
size = 500, date: 2024-01-17 02:46:47.164388, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:08.527615
size = 445, date: 2024-01-17 02:51:06.630610, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:08.413146
size = 590, date: 2024-01-17 02:59:28.946066, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S495'), ('column', 'S174'))] || Time: 0:00:19.683462
size = 500, date: 2024-01-17 03:08:12.489664, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:08.187292
size = 445, date: 2024-01-17 03:08:35.074837, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S338'), ('column', 'S335'))] || Time: 0:00:08.182074
size = 550, date: 2024-01-17 03:15:49.352160, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:37.677187
size = 500, date: 2024-01-17 03:30:27.746137, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S236'), ('column', 'S265'))] || Time: 0:00:12.112546
size = 550, date: 2024-01-17 04:12:00.908754, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:29.250125
size = 550, date: 2024-01-17 05:06:56.943108, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:31.580918
size = 550, date: 2024-01-17 05:42:03.160762, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:16.783561
size = 550, date: 2024-01-17 06:09:50.097218, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S162'), ('column', 'S241'))] || Time: 0:00:16.021064
# R^2 = 0.96459 up till here
size = 445, date: 2024-01-17 12:23:24.335463, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:17.341007
size = 550, date: 2024-01-17 12:28:25.545386, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:30.800783
size = 500, date: 2024-01-17 12:28:28.749300, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:24.947467
size = 590, date: 2024-01-17 12:30:40.356298, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:40.228692
size = 445, date: 2024-01-17 12:45:11.358852, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:05.189267
size = 500, date: 2024-01-17 12:57:22.051478, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:11.408032
size = 550, date: 2024-01-17 13:00:54.551925, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:13.921643
size = 445, date: 2024-01-17 13:01:38.816842, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:05.971833
size = 590, date: 2024-01-17 13:08:23.033370, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:16.647444
size = 445, date: 2024-01-17 13:18:45.077158, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:06.340460
size = 500, date: 2024-01-17 13:21:59.183902, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:16.173506
size = 550, date: 2024-01-17 13:26:25.895729, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:13.594930
size = 445, date: 2024-01-17 13:35:33.509092, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:05.854367
size = 590, date: 2024-01-17 13:37:21.588934, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:15.470711
size = 500, date: 2024-01-17 13:44:19.618802, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:11.634830
size = 550, date: 2024-01-17 13:52:20.508613, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:13.748687
size = 445, date: 2024-01-17 13:52:42.112531, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:06.234206
size = 590, date: 2024-01-17 14:07:40.712030, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:16.338387
size = 500, date: 2024-01-17 14:08:56.337631, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:19.747620
size = 445, date: 2024-01-17 14:09:44.922670, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:07.163878
size = 705, date: 2024-01-17 14:12:16.477458, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:01:21.075959
size = 550, date: 2024-01-17 14:19:29.592793, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:15.537419
size = 445, date: 2024-01-17 14:28:02.527726, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:05.863853
size = 500, date: 2024-01-17 14:31:31.911502, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:11.215296
size = 590, date: 2024-01-17 14:39:41.171713, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:15.235826
size = 445, date: 2024-01-17 14:44:09.402463, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S102'), ('column', 'S131'))] || Time: 0:00:05.377791
size = 550, date: 2024-01-17 14:44:59.084413, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:14.702003
size = 500, date: 2024-01-17 14:56:18.298426, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:09.840100
size = 590, date: 2024-01-17 15:08:46.372505, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:15.069421
size = 550, date: 2024-01-17 15:10:42.878968, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:14.182057
size = 705, date: 2024-01-17 15:11:32.077491, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:28.636197
size = 500, date: 2024-01-17 15:18:46.740594, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:10.104785
size = 550, date: 2024-01-17 15:36:15.638782, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:13.866278
size = 590, date: 2024-01-17 15:36:19.467657, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:17.706004
size = 500, date: 2024-01-17 15:43:15.408207, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S192'), ('column', 'S93'))] || Time: 0:00:12.746794
size = 705, date: 2024-01-17 15:58:30.515318, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:16.935491
size = 550, date: 2024-01-17 16:01:50.143086, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S150'), ('column', 'S295'))] || Time: 0:00:13.250565
size = 590, date: 2024-01-17 16:04:36.658812, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:00:18.716363
size = 705, date: 2024-01-17 16:46:54.776818, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:27.119095
size = 20, date: 2024-01-17 17:05:43.802298, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.386741
size = 20, date: 2024-01-17 17:06:01.202139, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.196912
size = 20, date: 2024-01-17 17:06:16.219749, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.131120
size = 20, date: 2024-01-17 17:06:32.633165, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.101849
size = 20, date: 2024-01-17 17:06:46.874774, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.103037
size = 20, date: 2024-01-17 17:07:00.717205, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.085312
size = 20, date: 2024-01-17 17:07:15.133767, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.124134
size = 20, date: 2024-01-17 17:07:29.975328, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.119238
size = 20, date: 2024-01-17 17:07:45.147591, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S13'), ('column', 'S6'))] || Time: 0:00:00.143476
size = 590, date: 2024-01-17 17:07:47.521772, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S575'), ('column', 'S385'))] || Time: 0:02:13.649622
size = 250, date: 2024-01-17 17:29:10.156683, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:46.638467
size = 705, date: 2024-01-17 17:33:25.401165, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:28.056312
size = 315, date: 2024-01-17 17:56:00.081168, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:01:00.593262
size = 385, date: 2024-01-17 18:16:57.742646, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:01:18.240566
size = 705, date: 2024-01-17 18:20:53.559592, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:16.381395
size = 250, date: 2024-01-17 18:36:33.672697, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:23.636155
size = 705, date: 2024-01-17 19:09:05.968527, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:34.610789
size = 250, date: 2024-01-17 19:30:55.023830, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:11.389286
size = 315, date: 2024-01-17 19:40:17.288635, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:34.108412
size = 705, date: 2024-01-17 20:00:57.080026, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:36.021363
size = 250, date: 2024-01-17 20:31:10.338348, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:12.131736
size = 385, date: 2024-01-17 20:54:09.403090, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:47.137183
size = 705, date: 2024-01-17 20:59:29.445194, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S252'), ('column', 'S329'))] || Time: 0:00:32.537930
size = 315, date: 2024-01-17 21:04:19.454733, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:26.601154
size = 250, date: 2024-01-17 21:13:04.816870, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:08.854785
size = 250, date: 2024-01-17 21:37:23.587825, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:08.762389
size = 315, date: 2024-01-17 21:43:59.059040, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:11.134169
size = 250, date: 2024-01-17 21:53:20.516625, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:06.698671
size = 385, date: 2024-01-17 21:53:39.495850, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:17.532634
size = 315, date: 2024-01-17 22:07:54.375831, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:10.855349
size = 250, date: 2024-01-17 22:08:01.052636, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:06.536523
size = 250, date: 2024-01-17 22:22:48.439195, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S96'), ('column', 'S142'))] || Time: 0:00:07.120637
size = 385, date: 2024-01-17 22:29:06.450264, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:15.969346
size = 315, date: 2024-01-17 22:31:23.199974, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:09.194015
size = 315, date: 2024-01-17 22:53:15.726700, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:11.090688
size = 385, date: 2024-01-17 23:03:10.448961, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:15.953154
size = 315, date: 2024-01-17 23:15:50.504368, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:10.916497
size = 385, date: 2024-01-17 23:36:08.179963, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:13.994047
size = 315, date: 2024-01-17 23:38:11.652902, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S215'), ('column', 'S39'))] || Time: 0:00:10.602560
size = 385, date: 2024-01-18 00:08:34.894520, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:14.289565
size = 385, date: 2024-01-18 00:39:51.362493, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:13.714797
size = 385, date: 2024-01-18 01:11:45.918407, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S121'), ('column', 'S19'))] || Time: 0:00:12.647496
size = 705, date: 2024-01-18 03:09:33.279253, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:01:15.704877
size = 835, date: 2024-01-18 03:37:59.723523, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:03:02.873987
size = 705, date: 2024-01-18 04:07:29.900435, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:39.044046
size = 705, date: 2024-01-18 04:52:55.787055, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:28.633031
size = 835, date: 2024-01-18 05:16:39.879521, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:05.416715
size = 705, date: 2024-01-18 05:36:36.681815, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:36.243618
size = 705, date: 2024-01-18 06:20:29.147539, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:35.283827
size = 835, date: 2024-01-18 06:31:43.128502, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:12.284479
size = 705, date: 2024-01-18 07:05:42.365762, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:28.427851
size = 705, date: 2024-01-18 07:49:20.714257, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:34.724114
size = 835, date: 2024-01-18 07:49:42.272129, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:07.863070
size = 705, date: 2024-01-18 08:32:32.966465, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:30.022477
size = 835, date: 2024-01-18 09:07:37.293338, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:07.863340
size = 705, date: 2024-01-18 09:17:43.127407, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S250'), ('column', 'S109'))] || Time: 0:00:34.812747
size = 835, date: 2024-01-18 10:25:41.145159, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:17.175151
size = 835, date: 2024-01-18 11:41:12.801822, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:21.256637
size = 835, date: 2024-01-18 12:56:37.239099, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:16.885753
size = 835, date: 2024-01-18 14:12:20.381018, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S354'), ('column', 'S213'))] || Time: 0:01:22.692231
"""

# after january 10, 2024 meeting
string_january_2024_ws_equal_1 = """
size = 1000, date: 2024-01-11 03:17:16.580444, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:05:23.648106
size = 975, date: 2024-01-11 03:27:43.153072, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:04:43.314642
size = 950, date: 2024-01-11 03:52:41.346127, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:04:44.711591
size = 920, date: 2024-01-11 04:02:37.476196, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:02:40.624743
size = 1000, date: 2024-01-11 05:35:09.570824, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:01:52.454912
size = 920, date: 2024-01-11 05:36:31.882341, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:00:55.724640
size = 975, date: 2024-01-11 05:40:14.233378, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:01:55.709836
size = 950, date: 2024-01-11 05:59:32.853625, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:01:27.381948
size = 920, date: 2024-01-11 06:50:38.612988, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:02.709156
size = 1000, date: 2024-01-11 07:25:48.253843, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:19.490106
size = 975, date: 2024-01-11 07:27:45.398085, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:01:50.666809
size = 950, date: 2024-01-11 07:36:57.713201, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:01:47.705630
size = 920, date: 2024-01-11 08:05:43.513338, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:09.922733
size = 975, date: 2024-01-11 09:15:24.412445, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:51.019034
size = 950, date: 2024-01-11 09:15:45.396741, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:01:59.373325
size = 920, date: 2024-01-11 09:23:32.411405, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:18.792308
size = 1000, date: 2024-01-11 09:24:03.973781, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:27.333804
size = 920, date: 2024-01-11 10:46:02.996288, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:19.942417
size = 950, date: 2024-01-11 10:55:38.780312, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:11.573664
size = 975, date: 2024-01-11 10:59:29.269693, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:01:12.727676
size = 1000, date: 2024-01-11 11:22:21.828305, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:03:38.670827
size = 920, date: 2024-01-11 12:07:36.982832, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:13.964590
size = 950, date: 2024-01-11 12:33:48.303024, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:01:38.846374
size = 975, date: 2024-01-11 12:41:39.794018, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:50.789102
size = 1000, date: 2024-01-11 13:18:24.728203, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:11.075606
size = 920, date: 2024-01-11 13:25:30.110225, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:17.449267
size = 950, date: 2024-01-11 14:15:09.763675, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:18.250788
size = 975, date: 2024-01-11 14:23:21.305926, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:53.799585
size = 920, date: 2024-01-11 14:44:44.552457, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:14.029542
size = 1000, date: 2024-01-11 15:13:29.978758, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:03:11.752764
size = 950, date: 2024-01-11 15:57:55.571717, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:09.661615
size = 920, date: 2024-01-11 16:06:54.278708, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:20.842082
size = 975, date: 2024-01-11 16:10:15.060340, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:01:57.084659
size = 895, date: 2024-01-11 16:19:01.212067, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:02:33.861742
size = 1000, date: 2024-01-11 17:03:56.254599, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:42.665936
size = 920, date: 2024-01-11 17:22:21.321143, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:08.090315
size = 950, date: 2024-01-11 17:40:25.973650, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:12.707688
size = 895, date: 2024-01-11 17:47:19.387967, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:00:29.851500
size = 975, date: 2024-01-11 17:51:57.816994, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:46.765612
size = 920, date: 2024-01-11 18:35:54.035866, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S876'), ('column', 'S716'))] || Time: 0:01:07.054219
size = 895, date: 2024-01-11 18:59:49.634041, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:19.510324
size = 1000, date: 2024-01-11 18:59:49.812891, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:03:14.679056
size = 950, date: 2024-01-11 19:22:17.094701, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:33.604107
size = 975, date: 2024-01-11 19:37:47.683061, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:52.141980
size = 895, date: 2024-01-11 20:17:37.751807, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:30.328465
size = 1000, date: 2024-01-11 20:58:22.067625, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:07:35.182795
size = 950, date: 2024-01-11 21:17:14.991769, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:12.838859
size = 895, date: 2024-01-11 21:33:32.639758, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:00:39.317492
size = 975, date: 2024-01-11 21:34:01.426797, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:49.548641
size = 895, date: 2024-01-11 22:44:19.317358, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:00:28.701510
size = 1000, date: 2024-01-11 22:54:39.242910, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:38.868691
size = 950, date: 2024-01-11 22:57:30.910456, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S18'), ('column', 'S122'))] || Time: 0:02:12.342116
size = 865, date: 2024-01-11 23:06:38.926263, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:03:12.505853
size = 975, date: 2024-01-11 23:17:18.966869, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S370'), ('column', 'S218'))] || Time: 0:00:49.924112
size = 895, date: 2024-01-12 00:00:10.734794, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:47.987209
size = 805, date: 2024-01-12 00:17:21.089837, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:02:56.015450
size = 835, date: 2024-01-12 00:19:28.927772, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:03:05.828073
size = 865, date: 2024-01-12 00:52:55.721669, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:03.912021
size = 1000, date: 2024-01-12 00:54:46.189418, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S642'), ('column', 'S529'))] || Time: 0:02:48.920383
size = 895, date: 2024-01-12 01:12:11.811062, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:17.109924
size = 805, date: 2024-01-12 01:49:55.213568, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:59.941412
size = 835, date: 2024-01-12 01:58:29.425526, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:03.940285
size = 865, date: 2024-01-12 02:16:27.658233, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:13.903997
size = 895, date: 2024-01-12 02:26:39.250599, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:15.371431
size = 775, date: 2024-01-12 02:35:39.358742, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:01:44.476059
size = 805, date: 2024-01-12 03:01:51.952336, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:01:02.475006
size = 835, date: 2024-01-12 03:19:02.807545, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:21.561147
size = 895, date: 2024-01-12 03:38:38.515029, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:19.345892
size = 775, date: 2024-01-12 03:41:35.053880, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:37.091615
size = 865, date: 2024-01-12 03:42:04.875487, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:18.081655
size = 805, date: 2024-01-12 04:15:30.951912, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:58.545510
size = 775, date: 2024-01-12 04:33:42.243029, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:45.294072
size = 835, date: 2024-01-12 04:37:46.407439, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:07.778143
size = 895, date: 2024-01-12 04:52:17.521906, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:17.383650
size = 865, date: 2024-01-12 05:04:10.788336, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:15.650692
size = 775, date: 2024-01-12 05:23:06.137709, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:41.220426
size = 805, date: 2024-01-12 05:26:37.449442, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:58.445079
size = 835, date: 2024-01-12 05:52:41.932610, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:09.852565
size = 895, date: 2024-01-12 06:02:43.465429, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S424'), ('column', 'S22'))] || Time: 0:01:02.485517
size = 775, date: 2024-01-12 06:14:29.288371, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:39.955261
size = 865, date: 2024-01-12 06:21:19.416494, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:01.778085
size = 805, date: 2024-01-12 06:32:01.637283, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:55.038193
size = 835, date: 2024-01-12 07:05:11.256671, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:04.463220
size = 775, date: 2024-01-12 07:06:31.366567, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:29.049570
size = 805, date: 2024-01-12 07:39:06.408418, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:01:05.003764
size = 865, date: 2024-01-12 07:39:49.700562, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:29.073790
size = 775, date: 2024-01-12 07:58:01.962500, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:39.112635
size = 835, date: 2024-01-12 08:17:25.184927, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:04.307772
size = 805, date: 2024-01-12 08:43:54.333483, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:01:08.000259
size = 775, date: 2024-01-12 08:49:33.331278, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:37.420275
size = 865, date: 2024-01-12 08:58:22.261043, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:27.463963
size = 835, date: 2024-01-12 09:27:52.869134, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:01.313148
size = 775, date: 2024-01-12 09:41:22.790915, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:56.791741
size = 805, date: 2024-01-12 09:50:34.913123, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:46.543705
size = 865, date: 2024-01-12 10:16:16.293570, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:27.788617
size = 775, date: 2024-01-12 10:31:42.809392, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:45.979343
size = 835, date: 2024-01-12 10:40:41.218017, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:05.286798
size = 805, date: 2024-01-12 10:59:57.904576, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:01:11.203720
size = 775, date: 2024-01-12 11:27:01.324142, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:37.880534
size = 865, date: 2024-01-12 11:55:15.472965, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:35.921804
size = 835, date: 2024-01-12 12:16:53.640205, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:08.065051
size = 775, date: 2024-01-12 12:18:50.372951, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S388'), ('column', 'S440'))] || Time: 0:00:33.169682
size = 805, date: 2024-01-12 12:26:43.799668, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:54.009739
size = 865, date: 2024-01-12 13:19:44.992928, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:36.049233
size = 835, date: 2024-01-12 13:33:04.755098, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:03.776201
size = 805, date: 2024-01-12 13:39:03.118162, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S363'), ('column', 'S32'))] || Time: 0:00:56.251772
size = 865, date: 2024-01-12 14:41:42.420080, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S263'), ('column', 'S852'))] || Time: 0:01:21.275256
size = 835, date: 2024-01-12 14:46:10.853401, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S650'), ('column', 'S151'))] || Time: 0:01:00.695706
size = 590, date: 2024-01-12 16:09:22.973105, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:36.747830
size = 590, date: 2024-01-12 16:44:15.277417, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:11.362817
size = 590, date: 2024-01-12 17:12:00.583928, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:11.240877
size = 590, date: 2024-01-12 17:39:56.059313, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:10.892242
size = 590, date: 2024-01-12 18:12:20.582289, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:13.606747
size = 590, date: 2024-01-12 18:44:57.535544, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:16.954185
size = 670, date: 2024-01-12 18:57:48.222435, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:01:06.065946
size = 740, date: 2024-01-12 19:04:14.158968, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:02:27.236164
size = 590, date: 2024-01-12 19:16:55.811920, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:13.668484
size = 630, date: 2024-01-12 19:31:13.082277, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:01:50.758250
size = 705, date: 2024-01-12 19:39:59.633765, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:01:48.914164
size = 590, date: 2024-01-12 19:48:57.012221, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:13.676338
size = 670, date: 2024-01-12 19:52:31.626884, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:21.212252
size = 590, date: 2024-01-12 20:18:13.266841, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:10.837042
size = 740, date: 2024-01-12 20:22:00.573346, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:00:39.097378
size = 630, date: 2024-01-12 20:29:54.388526, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:19.903750
size = 670, date: 2024-01-12 20:33:03.826442, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:26.046281
size = 705, date: 2024-01-12 20:44:56.052424, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:25.318838
size = 590, date: 2024-01-12 20:45:06.848026, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:11.075693
size = 630, date: 2024-01-12 21:08:50.422312, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:25.305821
size = 590, date: 2024-01-12 21:12:12.256011, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:09.537872
size = 670, date: 2024-01-12 21:17:53.184755, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:20.471357
size = 740, date: 2024-01-12 21:19:10.174613, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:00:57.077934
size = 705, date: 2024-01-12 21:32:20.446503, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:33.774811
size = 590, date: 2024-01-12 21:39:50.988771, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S330'), ('column', 'S366'))] || Time: 0:00:10.532304
size = 630, date: 2024-01-12 21:48:20.964613, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:23.241687
size = 670, date: 2024-01-12 22:01:52.820972, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:34.805962
size = 740, date: 2024-01-12 22:17:34.424871, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:00:58.008223
size = 705, date: 2024-01-12 22:19:02.982312, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:32.711010
size = 630, date: 2024-01-12 22:26:59.877945, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:23.547287
size = 550, date: 2024-01-12 22:41:30.636477, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:30.681754
size = 670, date: 2024-01-12 22:43:03.980018, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:29.494668
size = 630, date: 2024-01-12 23:06:28.898400, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:21.958189
size = 705, date: 2024-01-12 23:07:51.535146, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:31.252906
size = 740, date: 2024-01-12 23:14:05.019717, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:00:58.905468
size = 550, date: 2024-01-12 23:15:40.715415, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:11.861722
size = 670, date: 2024-01-12 23:23:51.484873, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:30.680071
size = 550, date: 2024-01-12 23:43:46.162964, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:16.293517
size = 630, date: 2024-01-12 23:45:42.056852, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:20.924882
size = 705, date: 2024-01-12 23:55:29.663245, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:31.550666
size = 670, date: 2024-01-13 00:02:44.208757, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:28.147573
size = 550, date: 2024-01-13 00:11:56.555702, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:14.609894
size = 740, date: 2024-01-13 00:12:53.377614, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S195'), ('column', 'S236'))] || Time: 0:00:45.342666
size = 630, date: 2024-01-13 00:24:28.762477, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:20.241889
size = 550, date: 2024-01-13 00:38:26.221156, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:14.291949
size = 705, date: 2024-01-13 00:41:02.516381, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:31.161942
size = 670, date: 2024-01-13 00:41:31.232192, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:34.679405
size = 630, date: 2024-01-13 01:02:37.546170, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:19.663878
size = 550, date: 2024-01-13 01:06:10.240485, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:12.970364
size = 670, date: 2024-01-13 01:20:24.190284, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:34.198855
size = 705, date: 2024-01-13 01:25:57.462101, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:30.670611
size = 550, date: 2024-01-13 01:31:34.884344, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:14.011506
size = 630, date: 2024-01-13 01:40:15.567512, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:19.122184
size = 550, date: 2024-01-13 01:57:13.208259, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:13.386505
size = 670, date: 2024-01-13 01:58:02.274374, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:32.815520
size = 705, date: 2024-01-13 02:11:16.767194, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:28.803158
size = 630, date: 2024-01-13 02:18:08.868108, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:19.044612
size = 550, date: 2024-01-13 02:24:58.412275, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:13.494045
size = 670, date: 2024-01-13 02:38:00.577252, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S302'), ('column', 'S474'))] || Time: 0:00:39.831423
size = 550, date: 2024-01-13 02:51:44.785556, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:12.849198
size = 630, date: 2024-01-13 02:54:41.353112, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S536'), ('column', 'S474'))] || Time: 0:00:19.042056
size = 705, date: 2024-01-13 02:56:14.751292, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:28.945473
size = 550, date: 2024-01-13 03:19:04.131484, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:14.186876
size = 705, date: 2024-01-13 03:42:36.602739, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:29.057060
size = 550, date: 2024-01-13 03:45:21.746220, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S117'), ('column', 'S439'))] || Time: 0:00:13.007850
size = 385, date: 2024-01-13 03:47:25.671116, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:09.834388
size = 445, date: 2024-01-13 03:50:00.793156, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:16.502870
size = 500, date: 2024-01-13 03:51:19.960403, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:20.815223
size = 315, date: 2024-01-13 03:56:23.519263, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:05.983036
size = 385, date: 2024-01-13 04:01:48.379672, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:05.189800
size = 315, date: 2024-01-13 04:05:52.926072, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.443599
size = 445, date: 2024-01-13 04:11:29.295331, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.455443
size = 385, date: 2024-01-13 04:13:34.470273, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:05.585057
size = 315, date: 2024-01-13 04:13:45.472904, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.567257
size = 500, date: 2024-01-13 04:15:58.820481, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:10.186175
size = 315, date: 2024-01-13 04:21:44.348502, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.228828
size = 385, date: 2024-01-13 04:25:06.367127, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.620736
size = 705, date: 2024-01-13 04:28:51.479344, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S46'), ('column', 'S538'))] || Time: 0:00:30.714192
size = 315, date: 2024-01-13 04:29:09.927825, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.268089
size = 445, date: 2024-01-13 04:29:15.208530, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:09.194639
size = 500, date: 2024-01-13 04:35:38.501465, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:11.635150
size = 385, date: 2024-01-13 04:36:19.120983, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.910109
size = 315, date: 2024-01-13 04:36:54.013467, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.647472
size = 315, date: 2024-01-13 04:44:47.996819, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.234954
size = 445, date: 2024-01-13 04:46:36.324721, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.221844
size = 385, date: 2024-01-13 04:47:24.504009, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.559604
size = 315, date: 2024-01-13 04:52:24.110161, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.478112
size = 500, date: 2024-01-13 04:55:28.263981, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:10.306595
size = 385, date: 2024-01-13 04:58:47.352016, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.534568
size = 315, date: 2024-01-13 05:00:15.335383, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.031739
size = 445, date: 2024-01-13 05:03:24.798218, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.700113
size = 315, date: 2024-01-13 05:08:08.391079, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.051837
size = 385, date: 2024-01-13 05:10:12.423919, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.420108
size = 500, date: 2024-01-13 05:15:15.290238, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:09.810449
size = 315, date: 2024-01-13 05:16:00.410865, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:02.959489
size = 445, date: 2024-01-13 05:21:15.901918, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.346985
size = 385, date: 2024-01-13 05:21:24.764150, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.423164
size = 315, date: 2024-01-13 05:23:59.952948, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S133'), ('column', 'S137'))] || Time: 0:00:03.002084
size = 385, date: 2024-01-13 05:32:30.252347, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.153274
size = 500, date: 2024-01-13 05:34:12.008453, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:09.479652
size = 445, date: 2024-01-13 05:38:11.131317, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:10.588200
size = 385, date: 2024-01-13 05:43:57.685006, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:04.192253
size = 500, date: 2024-01-13 05:53:37.810363, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:09.743278
size = 385, date: 2024-01-13 05:54:59.597556, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S67'), ('column', 'S7'))] || Time: 0:00:05.779353
size = 445, date: 2024-01-13 05:55:33.588470, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:10.586860
size = 445, date: 2024-01-13 06:12:44.724849, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.251387
size = 500, date: 2024-01-13 06:13:12.757749, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:12.601689
size = 445, date: 2024-01-13 06:30:07.625721, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.528175
size = 500, date: 2024-01-13 06:32:13.096923, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:12.454536
size = 445, date: 2024-01-13 06:47:43.259901, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.051396
size = 500, date: 2024-01-13 06:51:04.700831, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:11.900562
size = 445, date: 2024-01-13 07:05:19.779078, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S239'), ('column', 'S304'))] || Time: 0:00:08.243440
size = 500, date: 2024-01-13 07:10:33.907246, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:11.874638
size = 500, date: 2024-01-13 07:29:54.849657, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S266'), ('column', 'S39'))] || Time: 0:00:13.336473
size = 20, date: 2024-01-13 15:48:48.389301, Original Solution on RANDOMLY GENERATED MATRIX with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.073559
size = 20, date: 2024-01-13 15:48:50.682009, RANDOMLY GENERATED MATRIX iteration 0 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.028058
size = 20, date: 2024-01-13 15:48:52.584437, RANDOMLY GENERATED MATRIX iteration 1 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.022924
size = 20, date: 2024-01-13 15:48:54.442869, RANDOMLY GENERATED MATRIX iteration 2 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.011979
size = 20, date: 2024-01-13 15:48:56.309127, RANDOMLY GENERATED MATRIX iteration 3 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.023782
size = 20, date: 2024-01-13 15:48:58.136595, RANDOMLY GENERATED MATRIX iteration 4 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.023567
size = 20, date: 2024-01-13 15:49:00.016078, RANDOMLY GENERATED MATRIX iteration 5 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.025432
size = 20, date: 2024-01-13 15:49:01.871118, RANDOMLY GENERATED MATRIX iteration 6 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.024718
size = 20, date: 2024-01-13 15:49:03.712684, RANDOMLY GENERATED MATRIX iteration 7 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.022255
size = 20, date: 2024-01-13 15:49:05.511765, RANDOMLY GENERATED MATRIX iteration 8 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.026364
size = 20, date: 2024-01-13 15:49:07.326116, RANDOMLY GENERATED MATRIX iteration 9 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.022590
size = 20, date: 2024-01-13 15:49:09.178532, RANDOMLY GENERATED MATRIX iteration 10 with desired nash equilibria: [(('row', 'S5'), ('column', 'S11'))] || Time: 0:00:00.023283
"""

# string = string_january_2024_ws_equal_1
string = string_january_2024 

# remove the lines with "Original" in them
# string = '\n'.join([line for line in string.split('\n') if 'Original' not in line])

# remove empty lines
string = '\n'.join([line for line in string.split('\n') if line != ''])
# remove lines starting with #
string = '\n'.join([line for line in string.split('\n') if not line.startswith('#')])

# set font to calibri
# plt.rcParams["font.family"] = "Calibri"

# plot_averages(string, quadratic=False)
# plot_averages(string, quadratic=False)
# plot_averages(string, quadratic=True)
# without original solution, and upto date: 2023-11-21 17:22:14.283613, R^2 was 0.91
# plt.clf()
# plot_medians(string, display_text=True)
plot_medians(string)
# plot_averages(string, quadratic=False)

# # for every size, order all of the lines by date and print the final string
# final_string = ''
# for size in [20, 250, 315, 385, 445, 500, 550, 590, 630, 670, 705, 740, 775, 805, 835, 865, 895, 920, 950, 975, 1000]:
#     lines = [line for line in string.split('\n') if 'size = {}'.format(size) in line]
#     lines = sorted(lines, key=lambda x: x.split(', date: ')[1])
#     final_string += '\n'.join(lines) + '\n\n'
# print(final_string)


