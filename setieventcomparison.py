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
from PIL import Image

# Creating LWA-1 and LWA-SV dataframes

paths = ['LWAmetadataCSVs/'+x for x in os.listdir('LWAmetadataCSVs')]
for file in paths:
    if 'lwa1' in file:
        lwa1file = file
    elif 'lwa-sv' in file:
        lwasvfile = file
lwa1df = pd.read_csv(lwa1file)
lwasvdf = pd.read_csv(lwasvfile)

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

# Function to find duration of observation before processing files

def findduration(singlestringfile):
    obsposition = int(singlestringfile.find('obs_length: '))
    newlineposition = int(singlestringfile[obsposition:].find('\n'))+obsposition
    duration = float(singlestringfile[obsposition+12:newlineposition])
    return duration

lwa1duration = findduration(stringfile1)
lwasvduration = findduration(stringfile2)

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
svsnrlist = []
tophitlist = []
frequencylist = []
driftlist = []
snrlist = []
indexlist = []

for x in newdf1.index:
    frequency = newdf1['Corrected_Frequency'][x]
    fqhigh = frequency+0.05
    fqlow = frequency-0.05
    driftrate = newdf1['Drift_Rate'][x]
    tophitnumber = newdf1['Top_Hit_#'][x]
    snr = newdf1['SNR'][x]

    y = newdf2['Corrected_Frequency'].between(fqlow, fqhigh, inclusive='both')
    boolseries = y.any(axis=0, skipna=True)
    if boolseries == True:

        mfqindex = y[y.values].index
        mfqindexlength = len(mfqindex)

        for value in mfqindex:
            matchfq = newdf2.loc[value]
            dr = matchfq["Drift_Rate"].item()

            summ = dr+driftrate

            if abs((driftrate - dr)/(driftrate)) <= 0.2:


                print('\n', 'potential event found!', '\n', 'lwa1 row, fq, and dr: ', tophitnumber, \
                  frequency, driftrate,'\n', 'lwasv row and matching fq: ', matchfq)

                count += 1

                svtophit = newdf2['Top_Hit_#'].loc[value].item()
                svfrequency = newdf2['Corrected_Frequency'].loc[value].item()
                svdrift = newdf2['Drift_Rate'].loc[value].item()
                svsnr = newdf2['SNR'].loc[value].item()

                svtophitlist.append(svtophit)
                svfqlist.append(svfrequency)
                svdriftlist.append(svdrift)
                svsnrlist.append(svsnr)

                tophitlist.append(tophitnumber)
                frequencylist.append(frequency)
                driftlist.append(driftrate)
                snrlist.append(snr)

                indexlist.append(x)

if count == 0:
    print('\n', "no events found...")
    print('\n', "exiting script now")
    sys.exit()
else:
    print('\n', "event count: ", count, '\n')


# Displaying matching events as two dataframes (LWA-1 followed by LWA-SV)

adf = pd.DataFrame({'Telescope': 'LWA-1','Top_Hit_#': tophitlist,'Corrected_Frequency': frequencylist, \
                        'Drift_Rate': driftlist,'SNR': snrlist})
bdf = pd.DataFrame({'Telescope': 'LWA-SV','Top_Hit_#': svtophitlist,'Corrected_Frequency': svfqlist, \
                         'Drift_Rate': svdriftlist,'SNR': svsnrlist})

print("For the following tables of events, row n in table A corresponds to row n in table B. Together, they represent a potential event!")
print("LWA-1 event table:", adf)
print("LWA-SV event table:", bdf)

# Storing frequencies from events dataframes for future

afreqlist = []
adrlist = []
bfreqlist = []
bdrlist = []
for fq in adf['Corrected_Frequency']:
    afreqlist.append(fq)
for dr in adf['Drift_Rate']
    adrlist.append(dr)
for fq in bdf['Corrected_Frequency']:
    bfreqlist.append(fq)
for dr in bdf['Drift_Rate']:
    bdrlist.append(dr)

# Inputting our .h5 directory containing two files of specific row and tuning frequency

datadir2 = sys.argv[2]
h5paths = [datadir2 + '/' + x for x in os.listdir(datadir2) if x.endswith('.h5')]

h5pathlist  = order(h5paths)

# Creating waterfall plots for each event found and saving them as .png files to current directory

file_path = h5pathlist[0]
obs = Waterfall(file_path, max_load=11)

file_path1 = h5pathlist[1]
obs1 = Waterfall(file_path1, max_load=11)

filenamelist = []
for x in h5pathlist:
    slash = x.find('/')
    hyphen = x.find('-')
    filetype = x.find('.')
    if 'cor' in x:
        cor = x.find('cor')
        file = x[slash+1:hyphen]+"-"+x[filetype-4:filetype]+x[cor:cor+3]
        filenamelist.append(file)
    else:
        file = x[slash+1:hyphen]+"-"+x[filetype-4:filetype]
        filenamelist.append(file)

turbo = datadir.find('turbo')
rowtuning = datadir[:turbo]

for x in range(0, len(afreqlist)):
    afstart = afreqlist[x] - 0.005
    afstop = afreqlist[x] + 0.005
    adrift = adrlist[x]/(1e6)

    bfstart = bfreqlist[x] - 0.005
    bfstop = bfreqlist[x] + 0.005
    bdrift = bdrlist[x]/(1e6)

    obs.plot_waterfall(logged=True, f_start=afstart, f_stop=afstop)
    plt.title("LWA-1 event"+str(x)+": freq="+str(afreqlist[x])+", dr="+str(adrlist[x]))
    plt.plot([afreqlist[x], afreqlist[x]+(adrift*lwa1duration)], [0, lwa1duration], ls = (0, (3, 10)), label = 'Dotted Line', c = '#ff5f1f')
    pngname1 = str(filenamelist[0])+"event"+str(x)+"-lwa1.png"
    plt.savefig(fname=pngname1, bbox_inches='tight', pad_inches=0.5)
    plt.clf()

    obs1.plot_waterfall(logged=True, f_start=bfstart, f_stop=bfstop)
    plt.title("LWA-SV event"+str(x)+": freq="+str(bfreqlist[x])+", dr="+str(bdrlist[x]))
    plt.plot([bfreqlist[x], bfreqlist[x]+(bdrift*lwasvduration)], [0, lwasvduration], ls = (0, (3, 10)), label = 'Dotted Line', c = '#ff5f1f')
    pngnamesv = str(filenamelist[1])+"event"+str(x)+"-lwasv.png"
    plt.savefig(fname=pngnamesv, bbox_inches='tight', pad_inches=0.5)
    plt.clf()

    png1 = Image.open(pngname1)
    pngsv = Image.open(pngnamesv)

    len1, ht1 = png1.size
    lensv, htsv = pngsv.size

    lentotal = len1 + lensv
    height = max(ht1, htsv)

    eventfig = Image.new("RGB", (lentotal, height))
    eventfig.paste(png1, (0, 0))
    eventfig.paste(pngsv, (len1, 0))
    eventfig.save(fp=rowtuning+"event"+str(x)+".png")

print('\n', "plots have been generated and saved as .png files!")
