import os
import glob

def getFnames(dir, dtype="apr", minimumsize=7000.):
# def getFnames(dir, type="apr", minimumsize=7000.):

    fnames = []
#    os.chdir("../data/ChungCheonDC/")
    os.chdir("../data/IdongDC/")

    # print glob.glob("*.apr")
    for file in glob.glob("*.apr"):
        if os.path.getsize(file) > minimumsize:
            fnames.append(file)
    os.chdir("../../notebook")
    return fnames
