# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:20:58 2019

@author: ATPs
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

#filename = r"C:\Users\ATPs\OneDrive\Lab\Rutgers\other\20191025GaryClass\F20 01-447 Instructor Obligations - 102419.xlsx"

def plotAll(filename):
    '''
    filename is a excel file in the same format as the example we get
    '''

    # read excel file
    df = pd.read_excel(filename)
    # get list of Netid
    netids = list(df['Netid'])
    # generate empty sheet to store the time for each person
    
    weekdays = 'Sunday Monday Tuesday Wednesday Thursday Friday Saturday'.split()
    time_breaks = []
    for i in range(24):
        hour = '00' + str(i)
        hour = hour[-2:]
        time_breaks.append(hour+':00-'+hour+':29')
        time_breaks.append(hour+':30-'+hour+':59')
    time_breaks = time_breaks[15:-5]# not include all times in a day
    
    tdf = pd.DataFrame(np.zeros([len(time_breaks),len(weekdays)], dtype=int))
    tdf.index = time_breaks
    tdf.columns = weekdays
    
    # for each netid, get a dataframe to store the preference data.
    def getTimeForEach(preference):
        '''
        given a preference of a person, return a dataframe
        '''
        df_each = tdf.copy()
        times_avail = preference.strip(',').split(',')
        for t in times_avail:
            day, time_break = t.split('-', 1)
            day = weekdays[int(day) - 1]
            time_from, time_to = time_break.split('-')
            if len(time_from.split(':')[0]) != 2:
                time_from = '0' + time_from
            if len(time_to.split(':')[0]) != 2:
                time_to = '0' + time_to
            time_break = time_from + '-' + time_to
            df_each.loc[time_break, day] = 1
        return df_each
    
    
    dc_each = {}
    for row, r in df.iterrows():
        netid = r['Netid']
        preference = r['Preference data']
        dc_each[netid] = getTimeForEach(preference)
    
    # print a figure for each person
    def plotForEach(netid, outfolder=None):
        '''
        plot for each netid
        save pdf file in outfolder. if outfolder is None, it will be the folder of input excel file + 'pdf'
        '''
        if outfolder is None:
            outfolder = os.path.dirname(filename)
        if outfolder == '':
            outfolder = os.path.join('.', 'pdf')
        else:
            outfolder = os.path.join(outfolder, 'pdf')
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        outfile = os.path.join(outfolder,netid+'.pdf')
    #    netid = 'heiman'
        df_each = dc_each[netid]
        fig, ax = plt.subplots()
        im = ax.imshow(df_each.values,cmap=plt.cm.binary,aspect='auto')
        # all ticks...
        ax.set_xticks(np.arange(len(df_each.columns)))
        ax.set_yticks(np.arange(len(df_each.index)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(df_each.columns)
        ax.set_yticklabels(df_each.index)
        
        plt.setp(ax.get_xticklabels(),  ha="center")
        #add minor grid
        ax.set_xticks(np.arange(-.5, len(df_each.columns), 1), minor=True);
        ax.set_yticks(np.arange(-.5, len(df_each.index), 1), minor=True);
        ax.grid(which="minor",linestyle='-', linewidth=1, color='green', alpha=0.3)
        #add title
        ax.set(title=netid)
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='major',      # both major and minor ticks are affected
            bottom=True,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=True,
            length=0) # labels along the bottom edge are off
        fig.tight_layout()
        plt.savefig(outfile)
        plt.close()
    
    for netid in netids:
        plotForEach(netid)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a filename of a excel file, create pdf file showing times for each netid')
    parser.add_argument('file_excel', help = 'location of the excel file')
    f = parser.parse_args()
    plotAll(filename=f.file_excel)