import csv
import copy
import numpy as np
import math

lwa1tun1path = 'modifiedlwa1tun1.csv' #a
lwa1tun2path = 'modifiedlwa1tun2.csv'     #b
lwasvtun1path = 'modifiedlwasvtun1.csv'   #c
lwasvtun2path = 'modifiedlwasvtun2.csv'   #d

def process_csv(filename):
    example_file = open(filename, encoding="utf-8")
    example_reader = csv.reader(example_file)
    example_data = list(example_reader)
    example_file.close()
    return example_data

afile = process_csv(lwa1tun1path)
afileheader = afile[0]
afiledata = afile[1:]
asefd = process_csv('lwa1tun1sefd.csv')[1:]

bfile = process_csv(lwa1tun2path)
bfileheader = bfile[0]
bfiledata = bfile[1:]
bsefd = process_csv('lwa1tun2sefd.csv')[1:]

cfile = process_csv(lwasvtun1path)
cfileheader = cfile[0]
cfiledata = cfile[1:]
csefd = process_csv('lwasvtun1sefd.csv')[1:]

dfile = process_csv(lwasvtun2path)
dfileheader = dfile[0]
dfiledata = dfile[1:]
dsefd = process_csv('lwasvtun2sefd.csv')[1:]

sigma = 6
dstar = 10*(3.0857e16)
bandwidth = 2  # Hz
np = 2

def addeirp(filedata, fileheader, sefdfile, newcsvname):

    newfiledata = []
    newfileheader = []
    newfileheader.extend(fileheader)
    newfileheader.append('EIRP')

    count = 0
    for row in filedata:
        newrow = []
        newrow.extend(row)
        sefd = float(sefdfile[count][-1])*1000*(1e-26)
        duration = row[8]
        firstcolon = duration.find(':')
        minutes = float(duration[firstcolon+1:firstcolon+3])
        seconds = float(duration[-2:])
        tobs = minutes*60 + seconds
        eirp = sigma*4*math.pi*(dstar**2)*sefd*math.sqrt(bandwidth/(np*tobs))
        scieirp = '{:e}'
        scieirp = scieirp.format(eirp)
        newrow.append(scieirp)
        newfiledata.append(newrow)
        count += 1

    with open(newcsvname, 'w', newline='') as file:
        writeobject = csv.writer(file)
        writeobject.writerow(newfileheader)
        for row in newfiledata:
            writeobject.writerow(row)
    print("done creating ", newcsvname, "!")
    return

addeirp(afiledata, afileheader, asefd, 'lwa1tun1eirp.csv')
addeirp(bfiledata, bfileheader, bsefd, 'lwa1tun2eirp.csv')
addeirp(cfiledata, cfileheader, csefd, 'lwasvtun1eirp.csv')
addeirp(dfiledata, dfileheader, dsefd, 'lwasvtun2eirp.csv')
