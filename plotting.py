"""
Part of contactanalysis scripts: advanced feature analysis for plotting

@author: Markus Koerbel
@email: mk988@cam.ac.uk

Contact zones analysis of cells interacting with a SLB. Based on 3-colour TIR fluoresence microscopy tiff stacks. Plotting part, which derives further measures from 
the analysed images and outputs a csv which can be directly plotted with saeborn. 

"""

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import signal

def CA2D_plotting(cell_results, cell_summary, contact_results, QC, frame_smooth, pixel_size):
    """
    Takes analysed segmented data from contactanalysis and extracts parameters for plotting and comparison between conditions.
    Input:  cell_results: DataFrame for one condition

    """
    frame_smooth = np.int(frame_smooth)

    #assign unique cell numbering
    cell_IDs = np.arange(len(cell_summary))+1
    cell_summary['cell_ID'] = cell_IDs
    condlist = [(cell_results['cell'] == cell_summary['cell'][i]) & (cell_results['file'] == cell_summary['file'][i]) for i in range(len(cell_summary))]
    choicelist = [cell_summary['cell_ID'][i] for i in range(len(cell_summary))]
    cell_results['cell_ID'] = np.select(condlist, choicelist)
    condlist = [(contact_results['cell'] == cell_summary['cell'][i]) & (contact_results['file'] == cell_summary['file'][i]) for i in range(len(cell_summary))]
    contact_results['cell_ID'] = np.select(condlist, choicelist)

    cell_summary_good = cell_summary[(cell_summary['QC'] == QC)].copy()
    cell_results_good = cell_results[(cell_results['QC'] == QC)].copy()
    contact_results_good = contact_results[(contact_results['QC'] == QC)].copy() 
    cell_results_good['time_to_CZ_max [s]'] = -1
    cell_results_good['CZ_area_change [um2s-1]'] = -1
    cell_results_good['CCZ_area_change [um2s-1]'] = -1
    cell_results_good['centripedality'] = 0
    cell_summary_good['time_length_total [s]'] = -1

    for _, df in cell_summary_good.iterrows(): 
        time_interval = cell_results_good['time [s]'].iloc[-1] - cell_results_good['time [s]'].iloc[-2]
        selector = (cell_results_good['cell_ID'] == df.cell_ID)
        # smooth timetraces
        cell_results_good.loc[selector, 'CZ_area_total_smooth [um2]']   = cell_results_good.loc[selector, 'CZ_area_total [um2]'].rolling(frame_smooth, center=True).mean()
        cell_results_good.loc[selector, 'CCZ_area_total_smooth [um2]']  = cell_results_good.loc[selector, 'CCZ_area_total [um2]'].rolling(frame_smooth, center=True).mean()
        cell_results_good.loc[selector, 'CCZ_area_cnt_smooth']          = cell_results_good.loc[selector, 'CCZ_area_cnt'].rolling(frame_smooth, center=True).mean()
        cell_results_good.loc[selector, 'CZ_area_smooth_change [um2s-1]']  = cell_results_good.loc[selector, 'CZ_area_total_smooth [um2]'].diff() / time_interval
        cell_results_good.loc[selector, 'CCZ_area_smooth_change [um2s-1]'] = cell_results_good.loc[selector, 'CCZ_area_total_smooth [um2]'].diff() / time_interval

        i_res = cell_results_good[cell_results_good['cell_ID'] == df.cell_ID]
        # create DataFrame for begin and max timepoints
        CZ_smooth_max = i_res[i_res['CZ_area_total_smooth [um2]'] == np.amax(i_res['CZ_area_total_smooth [um2]'])].iloc[0]
        CCZ_smooth_max = i_res[i_res['CCZ_area_total_smooth [um2]'] == np.amax(i_res['CCZ_area_total_smooth [um2]'])].iloc[0]

        if CCZ_smooth_max['CCZ_area_cnt'] > 1:
            CCZ_smooth_begin = i_res[i_res['CCZ_area_cnt_smooth'] >= 1].iloc[0]

            # calculate CCZ growth rate
            CCZ_duration = CCZ_smooth_max['time [s]'] - CCZ_smooth_begin['time [s]']
            CCZ_growth = CCZ_smooth_max['CCZ_area_total_smooth [um2]'] - CCZ_smooth_begin['CCZ_area_total_smooth [um2]']
            CCZ_rate = CCZ_growth/CCZ_duration
            CCZ_max_rate = np.amax(i_res['CCZ_area_total_smooth [um2]'].diff()) / time_interval
            subwindow = i_res[(i_res['time [s]'] >= CCZ_smooth_begin['time [s]']) & (i_res['time [s]'] <= CCZ_smooth_max['time [s]'])]
            subwindow.loc[:,'lin_fit'] = subwindow['CCZ_area_total_smooth [um2]'].iloc[0] + np.arange(len(subwindow)) * CCZ_rate
            linearity = ((subwindow['lin_fit'] - subwindow['CCZ_area_total_smooth [um2]']) ** 2).mean() ** 0.5
        else:
            CCZ_smooth_begin = i_res.iloc[-1]
            CCZ_duration = -1
            CCZ_growth = -1
            CCZ_rate = -1
            CCZ_max_rate = -1
            linearity = -1

        # calculate timings for interaction stages
        cell_results_good.loc[selector, 'time_to_CZ_max [s]']   = cell_results_good.loc[selector, 'time_to_Ca [s]'] - CZ_smooth_max['time_to_Ca [s]']

        # add cell results to summary
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_stage_1 [s]'] = cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'],'time_CCZ_first [s]'] - cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'],'time_CZ_first [s]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_stage_2 [s]'] = cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'],'time_Ca [s]'] - cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'],'time_CCZ_first [s]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_stage_3 [s]'] = CZ_smooth_max['time_to_Ca [s]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_stage_4 [s]'] = cell_results_good.loc[cell_results_good['cell_ID']==df['cell_ID'], 'time_to_CZ_max [s]'].iloc[-1] #np.max(cell_results_good.loc[cell_results_good['cell_ID']==df['cell_ID'], 'time_to_CZ [s]']) - CZ_smooth_max['time_to_Ca [s]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_total [s]']   = cell_results_good.loc[selector, 'time_to_CZ [s]'].iloc[-1] #cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], ['time_length_stage_1 [s]', 'time_length_stage_2 [s]', 'time_length_stage_3 [s]', 'time_length_stage_4 [s]']].sum(axis=1)
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CZ_area_total_max [um2]'] = CZ_smooth_max['CZ_area_total [um2]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CCZ_area_total_max [um2]'] = CCZ_smooth_max['CCZ_area_total [um2]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CCZ_growth [um2]'] = CCZ_growth
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CCZ_growth_duration [s]'] = CCZ_duration
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_CCZ_growth_start [s]'] = CCZ_smooth_begin['time [s]']
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CCZ_rate [um2s-1]'] = CCZ_rate
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'CCZ_max_rate [um2s-1]'] = CCZ_max_rate
        cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'growth_linearity_rmsd'] = linearity
        
        # add time references to all_contacts
        contact_results_good.loc[contact_results_good['cell_ID'] == df.cell_ID, 'time_to_CCZ_first [s]'] = contact_results_good.loc[contact_results_good['cell_ID'] == df.cell_ID, 'time_to_Ca [s]'] + cell_summary_good.loc[cell_summary_good['cell_ID']==df['cell_ID'], 'time_length_stage_2 [s]'].iloc[0]
        contact_results_good.loc[contact_results_good['cell_ID'] == df.cell_ID, 'speed [ums-1]'] = np.sqrt(contact_results_good.loc[contact_results_good['cell_ID'] == df.cell_ID, 'x-displ [um]'] ** 2 + contact_results_good.loc[contact_results_good['cell_ID'] == df.cell_ID, 'y-displ [um]'] ** 2)/time_interval

        contact_results_good['radius [um]'] = np.sqrt(contact_results_good['area [um2]']/np.pi)


        # calculate centripedality
        def cosT(v1, v2):
            return np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))
        for i_frame in cell_results_good.loc[selector, 'frame'].values:
            i_contacts = contact_results_good.loc[(contact_results_good['cell_ID'] == df.cell_ID) & (contact_results_good['frame'] == i_frame) & (contact_results_good['contact'] == 'CCZ') & (contact_results_good['contact_time [s]']  > time_interval)]
            score = []
            center = cell_results_good.loc[(cell_results_good['cell_ID'] == df.cell_ID) & (cell_results_good['frame'] == i_frame), ['cell_center_x [px]', 'cell_center_y [px]']].values[0]*pixel_size
            for _, contact in i_contacts.iterrows():
                displ = np.array([contact['x-displ [um]'], contact['y-displ [um]']])
                c2c = np.array([contact['x-coord [um]'], contact['y-coord [um]']]) - displ - center #contact to center vector
                if not ((np.linalg.norm(displ) == 0) or (np.linalg.norm(c2c) == 0)):
                    score.append(cosT(displ, c2c))
                #else:
                #    print(i_contacts.cell_ID)
            cell_results_good.loc[(cell_results_good['cell_ID'] == df.cell_ID) & (cell_results_good['frame'] == i_frame), 'centripedality'] = np.median(score)

    condlist = [
        (cell_results_good['time_to_CCZ_first [s]'] < 0),
        ((cell_results_good['QC'] >= 'good') & (cell_results_good['time_to_CCZ_first [s]'] >= 0) & (cell_results_good['time_to_Ca [s]'] < 0)) | ((cell_results_good['QC'] >= 'good_noCa') & (cell_results_good['time_to_CCZ_first [s]'] >= 0)),
        (cell_results_good['QC'] >= 'good') & (cell_results_good['time_to_Ca [s]'] >= 0) & (cell_results_good['time_to_CZ_max [s]'] < 0), 
        (cell_results_good['QC'] >= 'good') & (cell_results_good['time_to_CZ_max [s]'] >= 0)]
    choicelist = ['scanning', 'searching', 'spreading', 'synapsing']
    cell_results_good['stage'] = np.select(condlist, choicelist)
    return cell_summary_good, cell_results_good, contact_results_good


def CA2D_plot_Ca_traces(summary, results, s_filter):
    """
    Take triggering and non-triggering cells and make Ca-trace and total CZ area plot
    Input: cell_summary, cell_results, and s_filter: number of timepoints to use for smoothing. 
    """

    n_trig = len(summary[summary.QC == 'good'])
    n_notrig = len(summary[summary.QC =='good_noCa'])
    n = np.max([n_trig, n_notrig])

    fig, ax = plt.subplots(n, 2, figsize=(20, 3*n), sharex=True)
    fig.suptitle('Condition {}'.format(summary.condition.iloc[0]), fontsize = 24)

    ax[0,0].set_title('Triggering cells')
    cell_good = summary.loc[summary.QC == 'good', 'cell_ID']
    for idx, i in enumerate(cell_good): #summary[summary.QC == 'good'].iterrows():
        cada = results.loc[results.cell_ID == i,:]
        sig_trace = results.loc[results.cell_ID == i, 'Ca_signal [AU]']
        smooth_trace = signal.wiener(sig_trace, s_filter)
        time_ax = results.loc[results.cell_ID == i, 'time_to_CCZ_first [s]']
        # plot signal trace
        sns.scatterplot(data = cada.loc[cada['time_to_Ca [s]'] == 0], x = 'time_to_CCZ_first [s]', y = 'Ca_signal [AU]', color='black', ax=ax[idx,0], s=75)
        sns.lineplot(data = cada, x = 'time_to_CCZ_first [s]', y = 'Ca_signal [AU]', legend = '', color = 'blue', ax = ax[idx,0], alpha = 0.7, lw=1)
        sns.lineplot(x = time_ax, y = smooth_trace, legend = '', color = 'blue', ax = ax[idx,0], lw=2)

        # plot CZ trace
        ax2 = ax[idx,0].twinx()
        sns.lineplot(data = cada, x = 'time_to_CCZ_first [s]', y = 'CZ_area_total [um2]', legend = '', color = 'green', ax = ax2, lw=2)
        ax[idx,0].legend(labels=['cell_ID {}'.format(i),])
    
    ax[0,1].set_title('Non triggering cells')
    cell_good_noCa = summary.loc[summary.QC == 'good_noCa', 'cell_ID']
    for idx, i in enumerate(cell_good_noCa): #summary[summary.QC == 'good_noCa'].iterrows():
        cada = results.loc[results.cell_ID == i,:]
        sig_trace = results.loc[results.cell_ID == i, 'Ca_signal [AU]']
        smooth_trace = signal.wiener(sig_trace, s_filter)
        time_ax = results.loc[results.cell_ID == i, 'time_to_CCZ_first [s]']
        # plot signal trace
        sns.lineplot(data = cada, x = 'time_to_CCZ_first [s]', y = 'Ca_signal [AU]', legend = '', color = 'blue', ax = ax[idx,1], alpha=0.7, lw=1)
        sns.lineplot(x = time_ax, y = smooth_trace, legend = '', color = 'blue', ax = ax[idx,1], lw=2)

        # plot CZ trace
        ax2 = ax[idx,1].twinx()
        sns.lineplot(data = cada, x = 'time_to_CCZ_first [s]', y = 'CZ_area_total [um2]', legend = '', color = 'green', ax = ax2, lw=2)
        ax[idx,1].legend(labels=['cell_ID {}'.format(i),])

    return fig