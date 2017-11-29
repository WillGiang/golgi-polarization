def orientation_percent(angles, want_ratio=False, dtheta=60):
    """Returns either a ratio or percentage
    
    angles: [array] of floats or ints
    want_ratio: [bool] Default False
    dtheta : [float] or [int]
    """
    is_oriented = np.abs(angles) <= dtheta
    n_oriented_angles = np.sum(is_oriented)
    
    ratio_oriented = n_oriented_angles/len(angles)
    
    if want_ratio:
        return ratio_oriented
    else:
        return ratio_oriented*100

def sample_proportion_error(p, n):
    """Returns the error for a sample proportion
    
    p: [float] 
        Proportion
    n: [int]
        Number of datapoints
    """
    return np.sqrt(p*(1-p) / n)

## Hypothesis test: difference between two proportions
## http://stattrek.com/hypothesis-test/difference-in-proportions.aspx

def pooled_sample_proportion(p1, p2, n1, n2):
    """Returns the pooled sample proportion between two proportions
    
    p1: [float] or [int]
        First proportion
    p2: [float] or [int]
        Second proportion
    n1: [int]
        Sample size associated with p1
    n2: [int]
        Sample size associated with p2
    """
    return (p1*n1 + p2*n2) / (n1 + n2)

def standard_error_between_two_proportions(p, n1, n2):
    """Returns the standard error between two proportions
    
    p: [float] or [int]
        Pooled sample proportion
    n1: [int]
        Sample size in first experiment
    n2: [int]
        Sample size in second experiment
    """
    return np.sqrt(p * (1 - p) * (1./n1 + 1./n2) )

def z_score_between_two_proportions(p1, p2, SE):
    """Returns the z-score between two proportions
    
    p1: [float] or [int]
        First proportion
    p2: [float] or [int]
        Second proportion
    SE: [float] or [int]
        Standard error between two proportions
    """
    return (p1 - p2) / SE



def plot_barchart(wild_orientation_percent, mutant_orientation_percent,
                            wild_errors, mutant_errors, wild_n, mutant_n,
                            colors, experiment, orientation_angle, plots_dir,
                            times, ax):
    
    n_groups = len(times) # 3
    index = np.arange(n_groups)
    bar_width = 0.4

    kwargs = {'alpha': 0.75,
              'capsize': 5,
              'edgecolor': 'k',
              'linewidth': 1
             }
    
    def autolabel(rects):
        """ Attach a text label on each bar displaying its height
    
         Shamelessly adapted from
         https://matplotlib.org/examples/api/barchart_demo.html
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.7*height,
                    '{}'.format(int(height)),
                    ha = 'center', va = 'bottom',
                    size = 15)
    
    rects_wild = ax.bar(index, wild_orientation_percent,
                        bar_width,
                        color = colors['Wild Type'],
                        label = 'Wild Type',
                        yerr = wild_errors,
                        **kwargs
                       )

    rects_mutant = ax.bar(index + bar_width, mutant_orientation_percent,
                          bar_width,
                          color = colors['Mutant DEE/DEE'],
                          label ='Mutant DEE/DEE',
                          yerr = mutant_errors,
                          **kwargs
                         )

    autolabel(rects_wild)
    autolabel(rects_mutant)
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Orientation %")
    ax.set_title("Golgi Polarization {} with {} deg".format(
        experiment, orientation_angle))
    ax.set_xticks(index + bar_width/2)
    ax.set_xticklabels(tuple(times))
    ax.set_ylim([0, 100])
    ax.legend(fontsize = 15)
    
    return None

def plot_table(wild_orientation_percent, wild_errors,
               mutant_orientation_percent, mutant_errors,
               wild_n, mutant_n, ax):
    
    table_columns = ('2 hrs', '4 hrs', '6 hrs', '2 hrs', '4 hrs', '6 hrs')
    rows = ['N', r'Oriented [%]', r'std err [%]']
    
    percentages = wild_orientation_percent + mutant_orientation_percent
    errors = wild_errors + mutant_errors
    table_data = [wild_n + mutant_n,
                  [int(x) for x in percentages],
                  [int(x) for x in errors]]
    ## We can cast the percentages above to ints since the errors are around 3% minimum
    
    table = ax.table(cellText = table_data,
                        colWidths = [0.1]*6,
                        colLabels = table_columns,
                        rowLabels = rows,
                        loc='center'
                       )
    
    table.set_fontsize(17)
    table.scale(1.5, 1.5)
    ax.axis('off')
    
    return None

def make_and_save_barchart_with_table(wild_orientation_percent, mutant_orientation_percent,
                            wild_errors, mutant_errors, wild_n, mutant_n, colors,
                            experiment, orientation_angle, plots_dir, times):
    plt.close('all')
    fig, axes = plt.subplots(2,1, figsize=(12,12))
    
    plot_table(wild_orientation_percent, wild_errors, mutant_orientation_percent, mutant_errors,
              wild_n, mutant_n, axes[0])
    
    
    plot_barchart(wild_orientation_percent, mutant_orientation_percent,
                            wild_errors, mutant_errors, wild_n, mutant_n, colors,
                            experiment, orientation_angle, plots_dir, times, axes[1])
    
    plt.savefig('{}{}_{}deg.png'.format(plots_dir, experiment, orientation_angle))
    plt.close(fig)
    return None


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.style.use('ggplot')
plt.rc('axes', titlesize = 22)
plt.rc('axes', labelsize = 22)
plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)

from glob import glob
import os

def main():
    ## Initialize experiment-related knowledge
    exp_names = ['exp 1', 'exp 3']

    unblinding_key = {'exp 1': {'Wild Type':'Animal 1',
                                'Mutant DEE/DEE':'Animal 2'}
                      ,
                      'exp 3': {'Wild Type':'Animal 2',
                                'Mutant DEE/DEE':'Animal 1'}
                     }

    csv_dir = r'./csvfiles/'
    plots_dir = r'./plots/'

    times = ['2hrs', '4hrs', '6hrs']

    colors = {'Wild Type': 'purple', 
              'Mutant DEE/DEE': 'blue'}

    ## Choice of angle here
    orientation_angle_list = np.arange(35, 65, 5)
    ##

    ## Initialize lists that will become dataframes later
    all_p_value_list = []
    all_times_list = []
    all_orientation_list = []
    all_exp_list = []

    ## Loop over multiple experiments
    for experiment in exp_names:
        ## Combine csv files into one pandas DataFrame
        csv_files = glob(os.path.join(csv_dir, experiment, "data*.csv"))
        df_in_each_file = (pd.read_csv(f) for f in csv_files)
        df = pd.concat(df_in_each_file, ignore_index=True)

        ## Add new columns to be used later
        df['animal_id'] = df.fileID.str.split('_', expand=True).get(0)
        df['time'] = df.fileID.str.split('_', expand=True).get(1).str.lower()

        new_df = df.filter(['golgiAreaInSquareMicrons',
                       'anglesInDeg',
                       'nucleiAreaInSquareMicrons',
                       'time',
                       'animal_id'],
                       axis=1
                          )

        ## area cut - keep data with golgi area between (10, 100) microns^2
        area_above_10_microns_squared = new_df['golgiAreaInSquareMicrons'] > 10
        area_below_100_microns_squared = new_df['golgiAreaInSquareMicrons'] < 100

        golgi_area_cut = area_above_10_microns_squared & area_below_100_microns_squared

        ## Split into wild type and mutant dataframes
        is_wild_type = new_df.animal_id == unblinding_key[experiment]['Wild Type']
        is_mutant = new_df.animal_id == unblinding_key[experiment]['Mutant DEE/DEE']

        wild_type = new_df[golgi_area_cut & is_wild_type]
        mutant = new_df[golgi_area_cut & is_mutant]

        for orientation_angle in orientation_angle_list:
            # Initialize list containers
            wild_n_list = []
            mutant_n_list = []
            wild_orientation_percent_list = []
            mutant_orientation_percent_list = []
            wild_error_list = []
            mutant_error_list = []

            for time in times:
                # Split by time
                wild_at_this_time = wild_type[wild_type.time == time]
                mutant_at_this_time = mutant[mutant.time == time]

                n_wild = len(wild_at_this_time)
                n_mutant = len(mutant_at_this_time)

                # Calculate the percent oriented for wild/mutant
                orientation_percent_wild = orientation_percent(wild_at_this_time['anglesInDeg'],
                                                              dtheta=orientation_angle)
                orientation_percent_mutant = orientation_percent(mutant_at_this_time['anglesInDeg'],
                                                                dtheta=orientation_angle)

                ratio_wild = orientation_percent_wild/100.
                ratio_mutant = orientation_percent_mutant/100.

                ## Rest of loop section is appending to lists
                wild_n_list.append(n_wild)
                mutant_n_list.append(n_mutant)

                wild_orientation_percent_list.append(orientation_percent_wild)
                mutant_orientation_percent_list.append(orientation_percent_mutant)

                wild_error_list.append(sample_proportion_error(ratio_wild, n_wild)*100)
                mutant_error_list.append(sample_proportion_error(ratio_mutant, n_mutant)*100)

                all_times_list.append(time)

            ## Generate plots for different angle orientations
            make_and_save_barchart_with_table(wild_orientation_percent_list,
                                              mutant_orientation_percent_list,
                                              wild_error_list, mutant_error_list,
                                              wild_n_list, mutant_n_list, colors,
                                              experiment, orientation_angle,
                                              plots_dir, times
                                             )

            ## Get numpy arrays of ratios and sample sizes
            p_wild = np.array(wild_orientation_percent_list)/100.
            p_mutant = np.array(mutant_orientation_percent_list)/100.

            n_wild = np.array(wild_n_list)
            n_mutant = np.array(mutant_n_list)

            ## We satisfy the criteria for Z tests
            p = pooled_sample_proportion(p_wild, p_mutant, n_wild, n_mutant)
            SE = standard_error_between_two_proportions(p, n_wild, n_mutant)
            z_scores = z_score_between_two_proportions(p_wild, p_mutant, SE)

            from scipy import special
            p_values = 1 - special.ndtr(z_scores)

            ## Rest is appending
            all_p_value_list.extend(p_values)
            all_orientation_list.extend([orientation_angle]*3)
            all_exp_list.extend([experiment]*3)


    ## Create pandas DataFrame
    d = {"p_value":all_p_value_list,
        "Orientation Angle":all_orientation_list,
        "Times":all_times_list,
        "Experiment":all_exp_list
        }

    df = pd.DataFrame(d)

    ## Split by experiment
    exp1 = df[df['Experiment'] == 'exp 1']
    exp3 = df[df['Experiment'] == 'exp 3']

    ## Change index to Orientation Angle
    ## Orientation Angle will be our 'x-axis' in our p-value plot
    exp1.set_index('Orientation Angle', inplace=True)
    exp3.set_index('Orientation Angle', inplace=True)
    
    ## Generate p-value plots
    plt.close('all')
    fig, axes = plt.subplots(2,1,figsize=(12,12), sharex=True)
    ax = axes.ravel()

    exp1.groupby('Times')['p_value'].plot(legend=True, ax=ax[0])

    ax[0].set_yscale('log')
    ax[0].set_title('exp 1')
    ax[0].set_xticks(np.arange(35, 65, 5))
    ax[0].set_ylabel('p-value')
    ax[0].axhline(0.05, ls="--", label='0.05')
    ax[0].legend()

    exp3.groupby('Times')['p_value'].plot(legend=True, ax=ax[1])

    ax[1].set_yscale('linear')
    ax[1].set_title('exp 3')
    ax[1].set_xlabel('Orientation Angle [deg]')
    ax[1].set_ylabel('p-value')
    ax[1].axhline(0.05, ls="--", label='0.05')
    ax[1].legend()

    plt.savefig('{}p_values_vs_angles.png'.format(plots_dir))