import pandas as pd 

def Threads_Location(chromfile,threads=2):
    Wind_List = []
    sep = 24 // threads
    integer = 24 % threads
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
            num_end = 24
        Wind_List.append(Chrom_List)
    return Wind_List
