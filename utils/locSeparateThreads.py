import pandas as pd 

def singleThreads_Location(chromfile):
    locate = pd.read_csv(chromfile)
    Chrom_List = []
    for pos in range(locate.shape[0]):
        Chrom_List.append([f'chr{locate["CHR"][pos]}',int(locate["STR"][pos]),int(locate["END"][pos])])
    return Chrom_List

def multiThreads_Location(chromfile,threads=2):
    file = open(chromfile,'r')
    lines = len(file.readlines()) - 1
    Wind_List = []
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

    locate = pd.read_csv(chromfile)

    for i in range(threads):
        Chrom_List = []
        for pos in range(locate.shape[0]):
            if pos >= num_start and pos <= num_end:
                Chrom_List.append([f'chr{locate["CHR"][pos]}',int(locate["STR"][pos]),int(locate["END"][pos])])
            else:
                continue
        num_start = num_end + 1
        if i < threads - 1:
            num_end = num_end + group[i+1] 
        else:
            num_end = lines
        Wind_List.append(Chrom_List)
    return Wind_List
