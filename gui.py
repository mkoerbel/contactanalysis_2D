import tkinter as tk
from tkinter import Tk, Frame, Scrollbar, Listbox, Button, filedialog, Checkbutton, Label
import os 
import pandas as pd

def get_filenames(dirname, folder_ids):
    sel_files = pd.DataFrame()
    with os.scandir(dirname) as entries:
        for dir in entries:
            image_path = 0
            # look for directories that have 'timelapse' in name; they should contain the tif file
            if folder_ids[0] in dir.name and dir.is_dir():
                #print(dir.name) #dir.path for full path
                # look for files in subfolder that match criteria (if)
                with os.scandir(dir) as files:
                    #print(len(list(files)))
                    for file in files:
                        if file.is_file() and file.name[-4:] == '.tif':
                            #print(file.name)
                            if image_path == 0:
                                image_path = file.path
                        if file.is_file() and file.name[-4:] == '.txt' and 'metadata' in file.name:
                            #print('Corresponding Metadata')
                            #print(file.name)
                            image_metadata = file.path
                        else:
                            image_metadata = []
                # get corresponding illumination file, to correct flat field
                illumination_dir = dir.path.replace(folder_ids[0], folder_ids[1])
                illumination_path = ''
                if os.path.isdir(illumination_dir):
                    #print(illumination_dir)
                    with os.scandir(illumination_dir) as files:
                        for file in files:
                            #print(file.name)
                            if file.is_file() and file.name[-4:] == '.tif':
                                #  print(file.name)
                                illumination_path = file.path
                #store filepaths in dataframe to close io before analysis
                #print(image_path, image_metadata, illumination_path)
                sel_files = sel_files.append({'timelapse': image_path, 
                                    'metadata': image_metadata, 
                                    'illumination': illumination_path, 
                                    'folder': dirname + '/' + dir.name}, ignore_index=True)
    return sel_files

def CA2D_gui(folder_ids):
    gui_list = []
    folder_list = []
    programs = {}
    window = Tk()
    window.title('Select input folders')

    def addFolder():
        dirname = filedialog.askdirectory()
        listbox.insert(tk.END, dirname)

    def deleteFolder():
        listbox.delete(listbox.curselection())

    def Start():
        gui_list.append(listbox.get(0, tk.END))
        programs['segmentation'] = b_segmentation.get()
        programs['analysis'] = b_analysis_cell.get()
        programs['plotting'] = b_analysis_contact.get()
        print('\nFolders selected:')
        for i in range(len(gui_list[0])):
            print('    {}. {}'.format(i+1, gui_list[0][i]))
        print('\n')
        window.destroy()

    frame0 = Frame(window)
    frame0.pack(padx = 5, pady = 5, fill = tk.X)
    capt = Label(frame0, text = 'Contactanalysis 2D')
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
    
    b_segmentation = tk.BooleanVar()
    b_analysis_cell = tk.BooleanVar()
    b_analysis_contact = tk.BooleanVar()

    frame2 = Frame(window, relief = tk.RAISED, borderwidth = 1)
    frame2.pack(padx = 10, pady = 10, fill = tk.X)
    btn_start = Button(frame2, text = 'START', command=Start, font =('Calibri', 10))
    btn_start.pack(side = tk.RIGHT, pady = 5, padx = 5)
    check_segment = Checkbutton(frame2, text = 'Contact Segmentation', variable = b_segmentation, font =('Calibri', 10))
    check_segment.pack(side = tk.LEFT, padx=5)
    check_cell = Checkbutton(frame2, text = 'Contact Analysis', variable = b_analysis_cell, font =('Calibri', 10))
    check_cell.pack(side = tk.LEFT, padx=5)
    check_contact = Checkbutton(frame2, text = 'Contact Plots', variable = b_analysis_contact, font =('Calibri', 10))
    check_contact.pack(side = tk.LEFT, padx=5)
    window.mainloop()

    for i_folder in gui_list[0]:
        folder_list.append([i_folder, get_filenames(i_folder, folder_ids)])

    return folder_list, programs