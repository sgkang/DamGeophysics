# from SimPEG import DCIP as DC
from SimPEG.EM.Static import DC
import numpy as np

def readReservoirDC(fname):
    f = open(fname, 'r')
    data = f.readlines()
    temp = data[3].split()
    nelec, ndata, aspacing = int(temp[0]), int(temp[1]), float(temp[2])
    height_water = float(data[4+ndata+3].split()[0])
    height_dam = float(data[4+ndata+4].split()[0])
    ntx = nelec-2
    datalist = []
    for iline, line in enumerate(data[4:4+ndata]):
    #     line = line.replace(ignorevalue, 'nan')
        linelist = line.split()
        datalist.append(np.array(map(float, linelist)))
    DAT = np.vstack(datalist)
    datalistSRC = []
    srcList = []
#     for i in range(ntx-1):
    for i in range(ntx-1):
        txloc = np.array([i+2, i+1.])
        ind = (DAT[:,:2] == txloc).sum(axis=1) == 2.
        temp = DAT[ind,:]
        datalistSRC.append(temp)
        e = np.zeros_like(temp[:,2])
        rxtemp = DC.Rx.Dipole(np.c_[temp[:,2]*aspacing, e, e], np.c_[temp[:,3]*aspacing, e, e])
        srctemp = DC.Src.Dipole([rxtemp], np.r_[txloc[1]*aspacing, 0., 0.], np.r_[txloc[0]*aspacing, 0., 0.])
        srcList.append(srctemp)
    DAT_src = np.vstack(datalistSRC)
    survey = DC.Survey(srcList)
    survey.dobs = DAT_src[:,-1]
    survey.height_water = height_water
    survey.height_dam = height_dam
    return survey

def readReservoirDC_all(fname):
    f = open(fname, 'r')
    data = f.readlines()
    temp = data[3].split()
    nelec, ndata, aspacing = int(temp[0]), int(temp[1]), float(temp[2])
    element = data[4+ndata+3].split()[0]
    try:
        height_water = float(data[4+ndata+3].split()[0])
        isheight = True
    except ValueError:
        isheight = False

    if isheight:
        height_water = float(data[4+ndata+3].split()[0])
        height_dam = float(data[4+ndata+4].split()[0])
    else:
        height_water = np.nan
        height_dam = np.nan

    ntx = nelec-2
    datalist = []
    ID = []
    for iline, line in enumerate(data[4:4+ndata]):
    #     line = line.replace(ignorevalue, 'nan')
        linelist = line.split()
        datalist.append(np.array(map(float, linelist)))
        ID.append(linelist[0]+linelist[1]+linelist[2]+linelist[3])
    DAT = np.vstack(datalist)
    height = np.r_[height_water, height_dam]
    return DAT,  height, ID


def readReservoirDC_data(fname):
    f = open(fname, 'r')
    data = f.readlines()
    temp = data[3].split()
    nelec, ndata, aspacing = int(temp[0]), int(temp[1]), float(temp[2])
    element = data[4+ndata+3].split()[0]


    ntx = nelec-2
    datalist = []
    for iline, line in enumerate(data[4:4+ndata]):
    #     line = line.replace(ignorevalue, 'nan')
        linelist = line.split()
        datalist.append(np.array(map(float, linelist)))

    DAT = np.vstack(datalist)
    return DAT
