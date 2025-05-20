from time import time

def Time_Stamp(start_time, last_time=0):
    if last_time == 0:
        run_time = time() - start_time
        print(f'Running time: {run_time:.2f} sec.')
    else:
        run_time = time() - last_time
        print(f'Running time: {run_time:.2f} sec.')
