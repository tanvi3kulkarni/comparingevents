# Import statements

from blimpy import Waterfall
import blimpy as bl
import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
import csv
import pandas as pd
import copy
from pathlib import Path
import sys

matplotlib.use('TkAgg')

# Creating LWA-1 and LWA-SV dataframes

paths = ['LWAmetadataCSVs/'+x for x in os.listdir('LWAmetadataCSVs')]
lwa1df = pd.read_csv(paths[0])
lwasvdf = pd.read_csv(paths[1])

# Mapping files across both dataframes

df1map = pd.DataFrame(lwa1df['datafile'])
df1map = df1map.rename(columns={'datafile':'lwa1file'})

dfsvmap = pd.DataFrame(lwasvdf['datafile'])
dfsvmap = dfsvmap.rename(columns={'datafile':'lwasvfile'})

mappingdf = pd.concat([df1map, dfsvmap], axis="columns")

def findotherfile(x):
    lwa1 = mappingdf['lwa1file']
    lwasv = mappingdf['lwasvfile']
    newlwa1 = lwa1.str.find(x, start=0, end=None)
    newlwasv = lwasv.str.find(x, start=0, end=None)
    if len(newlwa1[newlwa1 == 1]) == 1:
        lwa1index = newlwa1[newlwa1 == 1].index[0]
        lwasvfile = mappingdf['lwasvfile'].iloc[lwa1index]
        return ("lwa1 pandas row number: ", lwa1index), ("corresponding LWA-SV file: ", lwasvfile)
    elif len(newlwasv[newlwasv == 1]) == 1:
        lwasvindex = newlwasv[newlwasv == 1].index[0]
        lwa1file = mappingdf['lwa1file'].iloc[lwasvindex]
        return ("lwasv pandas row number: ", lwasvindex), ("corresponding LWA-1 file: ", lwa1file)
    else:
        return "file you inputted not found!"

# Inputting our .dat directory containing two files of specific row and tuning frequency

datadir = sys.argv[1]
pathlist1 = [datadir + '/' + x for x in os.listdir(datadir) if x.endswith('.dat')]

# Putting two files in order - LWA1 file followed by LWASV file

def order(pathlist):
    orderedpathlist = [0, 0]
    for x in pathlist:
        slash = x.find('/')
        hyphen = x.find('-')
        file = x[slash+1:hyphen]
        station, row = findotherfile(file)[0]
        if station.startswith('lwa1'):
            orderedpathlist[0] = x
        elif station.startswith('lwasv'):
            orderedpathlist[1] = x
        else:
            return "something's wrong here"
    return orderedpathlist

orderedpathlist = order(pathlist1)

# Opening files and loading contents as single string

openedfile1 = open(orderedpathlist[0], 'r')
stringfile1 = openedfile1.read()
openedfile1.close()

openedfile2 = open(orderedpathlist[1], 'r')
stringfile2 = openedfile2.read()
openedfile2.close()

# Functions to remove tabs, newlines, whitespaces, and everything except for header and data

def removeblankspace(filelistofstrings):
    while '' in filelistofstrings:
        filelistofstrings.remove('')
    return filelistofstrings

def removebeginning(listofstrings):
    if 'Top_Hit_#' in listofstrings:
        startpoint = listofstrings.index('Top_Hit_#')
        newlistofstrings = [listofstrings.pop(0) for x in range(0, startpoint)]
    return listofstrings

def removestuff(file):
    newfile = removebeginning(removeblankspace(file.rsplit(' ')))
    newfile.pop(12)
    return newfile

readylist1 = removestuff(stringfile1)
readylist2 = removestuff(stringfile2)

# Column names

headerlist1 = []
for x in readylist1[0:12]:
    if x.startswith('\t'):
        y = x[1:]
        headerlist1.append(y)
    else:
        headerlist1.append(x)

# Removing headers and sorting out data so that each element is a single data value

def datarows(readylist):

    rows = []
    startingrow = 12

    for x in range(startingrow, len(readylist)):

        if x == startingrow:
            nindex = readylist[startingrow].find('\n')
            tindex = readylist[startingrow].find('\t')
            row1_1 = readylist[startingrow][nindex+1:tindex]
            rows.append(row1_1)
        else:
            xstringlist = readylist[x].split('\t')
            for z in xstringlist:
                if z == '' or z == '\n':
                    continue
                elif z.startswith('\n'):
                    rows.append(z[1:])
                else:
                    rows.append(z)
    return rows

rows1 = datarows(readylist1)
rows2 = datarows(readylist2)

# Converting data to a pandas dataframe

def makedataframe(data_rows):

    dataframe = []
    singlerow = []
    for x in range(0, len(data_rows)):
        if x == 0:
            singlerow.append(data_rows[x])
        elif x%12 == 11:
            singlerow.append(data_rows[x])
            dataframe.append(copy.copy(singlerow))
            singlerow.clear()
        else:
            singlerow.append(data_rows[x])

    createddf = pd.DataFrame(dataframe, columns=headerlist1)
    createddf = createddf.apply(pd.to_numeric)
    return createddf

df1 = makedataframe(rows1)
df2 = makedataframe(rows2)

newdf1 = df1.loc[:, ['Top_Hit_#', 'Drift_Rate', 'SNR', 'Corrected_Frequency']]
newdf2 = df2.loc[:, ['Top_Hit_#', 'Drift_Rate', 'SNR', 'Corrected_Frequency']]

# Filtering rows from both dataframes to find *almost* matching frequencies and drift rates -> potential events

count = 0

svtophitlist = []
svfqlist = []
svdriftlist = []
tophitlist = []
frequencylist = []
driftlist = []
indexlist = []

for x in newdf1.index:
    frequency = newdf1['Corrected_Frequency'][x]
    fqhigh = frequency+0.002
    fqlow = frequency-0.002
    driftrate = newdf1['Drift_Rate'][x]
    tophitnumber = newdf1['Top_Hit_#'][x]

    y = newdf2['Corrected_Frequency'].between(fqlow, fqhigh, inclusive='both')
    boolseries = y.any(axis=0, skipna=True)
    if boolseries == True:

        mfqindex = y[y.values].index
        mfqindex2 = mfqindex.item()
        matchfq = newdf2.loc[mfqindex2]


        dr = matchfq["Drift_Rate"].item()

        summ = dr+driftrate

        if abs((driftrate - dr)/(driftrate)) <= 0.2:

            print('\n', 'potential event found!', '\n', 'lwa1 row, fq, and dr: ', tophitnumber, \
              frequency, driftrate,'\n', 'lwasv row and matching fq: ', matchfq)

            count += 1

            svtophit = newdf2['Top_Hit_#'].loc[mfqindex2].item()
            svfrequency = newdf2['Corrected_Frequency'].loc[mfqindex2].item()
            svdrift = newdf2['Drift_Rate'].loc[mfqindex2].item()

            svtophitlist.append(svtophit)
            svfqlist.append(svfrequency)
            svdriftlist.append(svdrift)

            tophitlist.append(tophitnumber)
            frequencylist.append(frequency)
            driftlist.append(driftrate)
            indexlist.append(x)

print('\n', "event count: ", count)

# Displaying matching events as two dataframes (LWA-1 followed by LWA-SV)

adf = newdf1.query("`Drift_Rate` in @driftlist and `Corrected_Frequency` in @frequencylist")
bdf = newdf2.query("`Drift_Rate` in @svdriftlist and `Corrected_Frequency` in @svfqlist")

print(adf)
print(bdf)

# Storing frequencies from events dataframes for future

adflist = []
bdflist = []
for fq in adf['Corrected_Frequency']:
    adflist.append(fq)
for fq in bdf['Corrected_Frequency']:
    bdflist.append(fq)

# Inputting our .h5 directory containing two files of specific row and tuning frequency

datadir2 = sys.argv[2]
h5paths = [datadir2 + '/' + x for x in os.listdir(datadir2) if x.endswith('.h5')]

h5pathlist  = order(h5paths)

# Creating waterfall plots for each event found and saving them as individual .png files to current directory

file_path = h5pathlist[0]
obs = Waterfall(file_path, max_load=11)

file_path1 = h5pathlist[1]
obs1 = Waterfall(file_path1, max_load=11)

filenamelist = []
for x in pathlist1:
    slash = x.find('/')
    hyphen = x.find('-')
    file = x[slash+1:hyphen]
    filenamelist.append(file)

for freq in range(0, len(adflist)):
    afstart = adflist[freq] - 0.5
    afstop = adflist[freq] + 0.5
    obs.plot_waterfall(logged=True, f_start=afstart, f_stop=afstop)
    plt.savefig(fname=str(filenamelist[0])+"event"+str(freq)+".png", bbox_inches='tight', pad_inches=1)
    plt.clf()
for freq in range(0, len(bdflist)):
    bfstart = bdflist[freq] - 0.5
    bfstop = bdflist[freq] + 0.5
    obs1.plot_waterfall(logged=True, f_start=bfstart, f_stop=bfstop)
    plt.savefig(fname=str(filenamelist[1])+"event"+str(freq)+".png", bbox_inches='tight', pad_inches=1)
    plt.clf()
