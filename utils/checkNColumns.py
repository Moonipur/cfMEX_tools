import gzip

def check_n_Columns(Input):
    with gzip.open(Input,'rt') as file:
        data = file.readline().strip().split('\t')
        if len(data) == 5:
            return 'c5'
        elif len(data) == 6:
            return 'c6'
