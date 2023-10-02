import os

import pandas as pd
from nucdens import access
import numpy as np
from pyperclip import copy


def main():
    # print all columns of the table with densities
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)

    # connect to database, choose working directory for downloading densities
    # workdir = os.environ["HOME"]+"/work/densitywork"
    workdir = os.environ["HOME"]+r"/OneDrive/"
    densdf = access.database(workdir=workdir)
    print("Proceeding with workdir="+workdir)
    myState = state(densdf, workdir)
    # print(myState)

    myState.Z, myState.N, myState.name = getNuc()
    myState.kind, myState.subfolder = getNumBodies()
    thetas = getTheta(myState)
    if len(thetas) > 1:
        print("Multiple Angles selected, using same energy for all angles")
        myState.theta = thetas[0]  # temporary
        myState.omega = getOmega(myState)
        for theta in thetas:
            myState.theta = theta
            # print(myState)
            # label is a dict with all the properties of the density
            myState.label = getlabel(myState)

            # print(densdf.pddf[selection][["Z", "N", "theta", "omega",
            #               "hashname"]])
            downloadFile(myState)
    else:
        myState.theta = thetas[0]
        myState.omega = getOmega(myState)
        myState.label = getlabel(myState)
        downloadFile(myState)

    os.remove(workdir+"densities_table.h5")


class state:
    def __init__(self, densdf, workdir):
        self.densdf = densdf
        self.workdir = workdir
        self.Z = None
        self.N = None
        self.name = None
        self.kind = None
        self.subfolder = None
        self.theta = None
        self.omega = None
        self.label = None

    def __str__(self):
        # this function is called when the "print" is called on it
        mDict = self.__dict__
        out = "{\n"
        for key, value in mDict.items():
            out += str(key)+":"+str(value)+"\n"
        out += "}"
        return out


def getNuc():
    inputStr = "\nInput which nucleus\n1. 3H\n2. 3He\n3. 4He\n4. 6Li\n"
    num = int(input(inputStr))
    if num == 1:
        Z = 1
        N = 2
        name = "3H"
    elif num == 2:  # 3He
        Z = 2
        N = 1
        name = "3He"
    elif num == 3:  # 4He
        Z = 2
        N = 2
        name = "4He"
    elif num == 4:  # 6Li
        Z = 3
        N = 3
        name = "6Li"
    else:
        raise ValueError("Wrong number entered")

    return Z, N, name


def getNumBodies():

    numBodies = int(
        input("\nEnter 1 or 2\n1. Onebody density\n2. Twobody density\n"))
    if numBodies == 1:
        kind = "one"
        subfolder = "1Ndensities"
    elif numBodies == 2:
        kind = "two"
        subfolder = "2Ndensities"
    else:
        raise ValueError("Wrong number entered")
    return kind, subfolder


def getTheta(mState):
    densdf = mState.densdf
    s = ((densdf.pddf.Z == mState.Z) & (
        densdf.pddf.kind == mState.kind) & (
        densdf.pddf.N == mState.N))

    printframe = densdf.pddf[s][["Z", "N", "theta"]
                                ].sort_values("theta")
    thetaValues = np.unique(printframe["theta"].to_numpy())
    print("\n")
    print("Possible Angles Are")
    print("------------------------")
    for i in range(len(thetaValues)):
        print(i+1, ". ", thetaValues[i])

    print("For multiple angles enter integer on left with a space 1 2 3")
    i = np.array(input("Input integer for the angle desired: ").split())
    i = i.astype(int)
    i = i-1
    theta = thetaValues[i]
    print("Selected theta=", theta, "\n")
    return theta


def getOmega(mState):
    densdf = mState.densdf
    s = ((
        densdf.pddf.Z == mState.Z) & (
        densdf.pddf.kind == mState.kind) & (
        densdf.pddf.theta == mState.theta) & (
        densdf.pddf.N == mState.N)
    )

    printframe = densdf.pddf[s][["Z", "N", "omega"]
                                ].sort_values("omega")

    omegaValues = np.unique(printframe["omega"].to_numpy())
    print("Possible Energy Values Are")
    print("------------------------")
    for i in range(len(omegaValues)):
        print(i+1, ". ", omegaValues[i], "MeV")
    i = int(input("Input integer for the energy desired: "))-1
    omega = omegaValues[i]
    return omega


def getlabel(mState):
    """
    Returns a dictionary that can be passed to densdf.get_file(**label)
    in order to actually do the download
    """
    densdf = mState.densdf
    selection = ((
        densdf.pddf.Z == mState.Z) & (
        densdf.pddf.N == mState.N) & (
        densdf.pddf.theta == mState.theta) & (
        densdf.pddf.kind == mState.kind) & (
        densdf.pddf.omega == mState.omega)
    )

    # tmp = densdf.pddf[selection]

    if len(densdf.pddf[selection]) == 1:
        label = densdf.pddf[selection].to_dict("records")[0]
    elif len(densdf.pddf[selection]) > 1:
        selectList = ["kind", "Z", "N", "omega", "theta",
                      "uniquefilename", "addtime",
                      "moditime", "LambdaNN", "tnforder"]
        print(densdf.pddf[selection][selectList])

        print("\nMore than one entry selected. Printed all columns")
        row = int(input("Enter row (number in left) of file to use: "))
        label = densdf.pddf.to_dict("records")[row]
    else:
        raise ValueError("No density found with these paramters")
    return label


def downloadFile(mState):
    name = mState.name
    angle = mState.theta
    omega = mState.omega
    kind = mState.kind
    subfolder = mState.subfolder
    densdf = mState.densdf
    label = mState.label
    workdir = mState.workdir

    densName0 = name + "-"+str(angle)+"=theta-"+str(omega)\
        + "=MeV-"+kind+"body.h5"

    pot = label["MODENN"]
    LambdaNN = label["LambdaNN"]
    densName1 = name + "-"+str(angle)+"=theta-"+str(omega)\
        + "=MeV-"+kind+"body-"\
        + str(pot)+"=poten-"\
        + str(LambdaNN)+"=LambdaNN"\
        + ".h5"

    densNames = [densName0, densName1]

    print("Download file with which name")
    print("------------------------------")
    for i in range(len(densNames)):
        print(i+1, ". ", densNames[i])
    i = int(input("Input integer for name desired: "))-1
    densName = densNames[i]

    path = "densities-"+name+r"/"+subfolder+r"/"
    if not os.path.exists(path):
        print("Directory named\n"+path+"\nDoes not exist, creating one now")
        os.makedirs(path)

    fullPath = path+densName
    if os.path.isfile(fullPath):
        yesno = input("File already exists, okay to overwrite? y/n: ").lower()
        if yesno == "n":
            os.remove("densities_table.h5")
            return 0
    # full print uniquefilename
    hashname, uniquename = densdf.get_file(**label)
    assert(os.path.isfile(hashname))
    fileLocation = workdir+fullPath
    os.replace(hashname, fileLocation)
    print("Downloaded file to: " + fileLocation)
    copy("'"+fileLocation+"'")
    print("Copied file location to clipboard")


if __name__ == "__main__":
    main()
