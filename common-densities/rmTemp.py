# -*- coding: utf-8 -*-

"""
@author: alexl
"""
import os
# import numpy as np


def main():
    filesHere = [f for f in os.listdir() if os.path.isfile(f)]
    filesInSub = [f for f in os.listdir(
        r"./density-modules/") if os.path.isfile(f)]
    intersection = list(set(filesHere) & set(filesInSub))
    intersection.sort()
    print(intersection)
    for file in intersection:
        os.remove(file)


if __name__ == "__main__":
    main()
