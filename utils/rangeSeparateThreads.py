import pandas as pd 

def singleThreads_Range(chromfile):
    cfrange = pd.read_csv(chromfile)
    Range_List = []
    for pos in range(cfrange.shape[0]):
        Range_List.append([int(cfrange["START"][pos]),int(cfrange["END"][pos])])
    return Range_List

def multiThreads_Range(chromfile,threads=2):
    file = open(chromfile,'r')
    lines = len(file.readlines()) - 1
    Total_List = []
    sep = lines // threads
    integer = lines % threads
    group = []
    for i in range(threads):
        if integer > 0:
            element = sep + 1
        elif integer <= 0:
            element = sep
        integer -= 1
        group.append(element)

    num_start = 0
    num_end = group[0] - 1

    cfrange = pd.read_csv(chromfile)

    for i in range(threads):
        Range_List = []
        for pos in range(cfrange.shape[0]):
            if pos >= num_start and pos <= num_end:
                Range_List.append([int(cfrange["START"][pos]),int(cfrange["END"][pos])])
            else:
                continue
        num_start = num_end + 1
        if i < threads - 1:
            num_end = num_end + group[i+1] 
        else:
            num_end = lines
        Total_List.append(Range_List)
    return Total_List
