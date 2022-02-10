"""
Created on Jul 18 2021

@author: Markus Koerbel
@email: mk988@cam.ac.uk

Uses output of contactanalysis_2D to estimate the probability of triggering as a function of a feature (contact radius or duration) to compare with modelling.
    Two probability functions will be calculated, based on the assumption which contact is responsible for triggering: 
        The two simplified assumptions are that first, the biggest contact is responsible for triggering, or 
        second, the contact with the longest duration is responsible. 

    If two cells, one triggering and one non triggering, are considered, the probability is calculated as follows:
    Bins of 'feature' are created. At the timepopint of triggering the biggest(longest lasting) contact is considered, and entries added to all bins greater than the given feature size.
    If the cell does not trigger within the observation period, the last frame is considered and only counts towards the total number of cells for bins smaller than the feature size. 
    For each bin the probability of triggering is then the (sum of entries) / (number of entries). That way cells are only considered within their observation period.  

    bin     |    1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   |
    Cell 1       0       0       0       1       1       1       1       1      
    Cell 2       0       0       0       0       -       -       -       -

    Probability  0       0       0      0.5      1       1       1       1  

INPUT 
    - Folder with the same folder structure as contactanalysis_2D
        . Includes output of contactanalysis_2D 'Plotting'

OUTPUT
    - csv file with triggering probability

"""

# USER SETTINGS

# Input Settings
feature = 'radius [um]' # choose between 'radius [um]' and 'contact_time [s]' as the x-axis
bins = [0,2000, 50] # either contact radius in nm OR contact time in s (depends on chosen feature), [start, end, bin_width]



###########################################
# IMPORT MODULES
# general modules
import numpy as np 
import pandas as pd

import tkinter as tk
from tkinter import Tk, Frame, Scrollbar, Listbox, Button, filedialog, Checkbutton, Label
import os

# Custom functions
from gui import get_filenames

def contactanalysisGui():
    gui_list = []
    plots = {}
    window = Tk()
    window.title('Select input folders')

    def addFolder():
        dirname = filedialog.askdirectory()
        listbox.insert(tk.END, dirname)

    def deleteFolder():
        listbox.delete(listbox.curselection())

    def Start():
        gui_list.append(listbox.get(0, tk.END))
        plots['triggering_probability'] = b_trigprob.get()
        plots['analysis'] = b_analysis_cell.get()
        plots['plotting'] = b_analysis_contact.get()
        print(gui_list)
        window.destroy()

    frame0 = Frame(window)
    frame0.pack(padx = 5, pady = 5, fill = tk.X)
    capt = Label(frame0, text = 'Contactanalysis 2D - Modelling Plots')
    capt.config(font =('Calibri', 14))
    capt.pack(side = tk.LEFT, padx = 10, pady = 5)

    frame1 = Frame(window, relief = tk.RAISED, borderwidth = 1)
    frame1.pack(padx = 10, pady = 10, fill = tk.BOTH)

    btn_del = Button(frame1, text = 'DEL', command=deleteFolder, font =('Calibri', 10))
    btn_del.pack(side = tk.RIGHT, padx=5)
    btn_add = Button(frame1, text = 'ADD', command=addFolder, font =('Calibri', 10))
    btn_add.pack(side = tk.RIGHT, padx=5)

    scrollbar = Scrollbar(frame1, orient = tk.HORIZONTAL)
    listbox = Listbox(frame1, xscrollcommand = scrollbar.set)
    scrollbar.config(command=listbox.xview)
    scrollbar.pack(side = tk.BOTTOM, fill=tk.X)
    listbox.pack(side = tk.LEFT, fill = tk.BOTH, expand = tk.TRUE)
    
    b_trigprob = tk.BooleanVar()
    b_analysis_cell = tk.BooleanVar()
    b_analysis_contact = tk.BooleanVar()

    frame2 = Frame(window, relief = tk.RAISED, borderwidth = 1)
    frame2.pack(padx = 10, pady = 10, fill = tk.X)
    btn_start = Button(frame2, text = 'START', command=Start, font =('Calibri', 10))
    btn_start.pack(side = tk.RIGHT, pady = 5, padx = 5)
    check_segment = Checkbutton(frame2, text = 'Triggering Probability', variable = b_trigprob, font =('Calibri', 10))
    check_segment.pack(side = tk.LEFT, padx=5)
    check_cell = Checkbutton(frame2, text = 'TBD', variable = b_analysis_cell, font =('Calibri', 10))
    check_cell.pack(side = tk.LEFT, padx=5)
    check_contact = Checkbutton(frame2, text = 'TBD', variable = b_analysis_contact, font =('Calibri', 10))
    check_contact.pack(side = tk.LEFT, padx=5)
    window.mainloop()

    return gui_list[0], plots

def contact_probability(summary, results, contacts, bins, feature):
    """
    Calculate triggering probability for feature (contact radius or duration) to compare with modelling.
    Two probability functions will be calculated, based on the assumption which contact is responsible for triggering: 
        The two simplified assumptions are that first, the biggest contact is responsible for triggering, or 
        second, the contact with the longest duration is responsible. 

    If two cells, one triggering and one non triggering, are considered, the probability is calculated as follows:
    Bins of 'feature' are created. At the timepopint of triggering the biggest(longest lasting) contact are considered, and entries added to all bins after that.
    If the cell does not trigger within the observation period, the last frame is considered and only counts towards the total number of cells for that bin. 
    For each bin the probability of triggering is then the (sum of entries) / (number of entries). That way cells are only considered within their observation period.  

    bin     |    1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   |
    Cell 1       0       0       0       1       1       1       1       1      
    Cell 2       0       0       0       0       -       -       -       -

    Probability  0       0       0      0.5      1       1       1       1  

    """

    # prepare sampling interval size
    bin_width = bins[2]
    n_bins = np.int(np.ceil(bins[1]/bin_width))
    # bin bounds
    end_bins = np.arange(bin_width,bin_width*(n_bins+1), bin_width)
    # prepare arrays
    # two counts will be made, for each assumption. r (contact radius), t (contact time)
    activated_counts_r = np.zeros(n_bins, dtype=np.int)
    activated_counts_t = np.zeros(n_bins, dtype=np.int)
    total_counts_r = np.zeros(n_bins, dtype=np.int)
    total_counts_t = np.zeros(n_bins, dtype=np.int)
    n_cells = len(summary.cell_ID)

    for i, cell in summary.iterrows():
        if cell.QC == 'good':
            trig_contacts = contacts.loc[(contacts.cell_ID == cell.cell_ID) & (contacts['time_to_Ca [s]'] == 0) & (contacts['contact'] == 'CCZ'), :]
            if len(trig_contacts) > 0:
                max_r_ = np.max(trig_contacts['radius [um]'].values)
                max_r = np.max(trig_contacts.loc[trig_contacts['radius [um]'] == max_r_, feature])
                max_t_ = np.max(trig_contacts['contact_time [s]'].values)
                max_t = np.max(trig_contacts.loc[trig_contacts['contact_time [s]'] == max_t_, feature])
            else:
                max_r = 0
                max_t = 0
            if feature == 'radius [um]':
                max_r = max_r*1000
                max_t = max_t*1000
            activated_counts_r += end_bins >= max_r
            activated_counts_t += end_bins >= max_t
            total_counts_r += 1
            total_counts_t += 1

        elif cell.QC == 'good_noCa':
            last_frame = np.amax(results.loc[results.cell_ID == cell.cell_ID, 'frame'])
            last_contacts = contacts.loc[(contacts.cell_ID == cell.cell_ID) & (contacts['frame'] == last_frame) & (contacts['contact'] == 'CCZ'), :]
            max_r_ = np.max(last_contacts['radius [um]'].values)
            max_r = np.max(last_contacts.loc[last_contacts['radius [um]'] == max_r_, feature])
            max_t_ = np.max(last_contacts['contact_time [s]'].values)
            max_t = np.max(last_contacts.loc[last_contacts['contact_time [s]'] == max_t_, feature])
            if feature == 'radius [um]':
                max_r = max_r*1000
                max_t = max_t*1000
            total_counts_r += end_bins-bin_width < max_r
            total_counts_t += end_bins-bin_width < max_t 

    return end_bins, activated_counts_r, activated_counts_t, total_counts_r, total_counts_t

#####################################
# MAIN

# have user input list of folder and tick which part of the program to run
[folders, plots] = contactanalysisGui()
print(str(len(folders)) + ' condition(s) selected. Lets plots this!')

for dirname in folders:        # dirname of inmput folder. One folder is one condition
    condition = os.path.basename(dirname)
    print(' Plotting condition ' + str(condition))

    if plots['triggering_probability']:
        cell_summary_plot = pd.read_csv(dirname + '/cell_summary_plot.csv')
        cell_results_plot = pd.read_csv(dirname + '/cell_results_plot.csv')
        contact_results_plot = pd.read_csv(dirname + '/contact_results_plot.csv')

        end_bins, activated_counts_r, activated_counts_t, total_counts_r, total_counts_t = contact_probability(cell_summary_plot, cell_results_plot, contact_results_plot, bins, feature)

        if feature == 'radius [um]':
            output = pd.DataFrame({
                'radius_bin_end [nm]': end_bins,
                'triggering_max_r': activated_counts_r,
                'triggering_max_t': activated_counts_t,
                'ncells_max_r': total_counts_r,
                'ncells_max_t': total_counts_t,
                'probability_max_r': activated_counts_r/total_counts_r,
                'probability_max_t': activated_counts_t/total_counts_t})
            output.to_csv(dirname + '\contact_probability_max_r.csv')
        elif feature == 'contact_time [s]':
            output = pd.DataFrame({
                'time_bin_end [s]': end_bins,
                'triggering_max_r': activated_counts_r,
                'triggering_max_t': activated_counts_t,
                'ncells_max_r': total_counts_r,
                'ncells_max_t': total_counts_t,
                'probability_max_r': activated_counts_r/total_counts_r,
                'probability_max_t': activated_counts_t/total_counts_t})
            output.to_csv(dirname + '\contact_probability_max_t.csv')
        else:
            print('No valid \'feature\' input!')