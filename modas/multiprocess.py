import multiprocessing as mp
import subprocess
import os

def run(cmd):
    """
    run linux shell by subprocess
    """
    null_f = open(os.devnull, 'w')
    a = subprocess.call(cmd, shell=True, stdout=null_f, stderr=null_f)
    null_f.close()

    return a


def parallel(run, values, num_threads=1):
    """parallel run linux command line"""
    try:
        p = mp.Pool(processes=num_threads)
        results = [p.apply_async(run, cmd) for cmd in values]
        p.close()
        p.join()

    except KeyboardInterrupt as e:
        p.terminate()
        p.join()
        raise e
    return [result.get() for result in results]
