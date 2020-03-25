import os
import time


def time_elapsed(fun):
    """Function to be used as a decorator for measuring time elapsed."""
    def wrapper(*args, **kwargs):
        t0 = time.time()
        fun(*args, **kwargs)
        print("\nTime elapsed = %s" % (time.time() - t0))
    return wrapper


def elements(filelist, start, end):
    """Generator over the lines in filelist, only yielding lines from line niumber
    start up to and including end."""
    process = False
    n = 1
    with open(filelist) as fh:
        for line in fh:
            fname = line.strip()
            if n == start:
                process = True
            if n > end:
                return
            if process:
                yield (n, fname)
            n += 1


def print_element(n, fname):
    print("%s  %07d  %s" % (time.strftime("%Y%m%d:%H%M%S"), n, fname))


def ensure_directory(*fnames):
    """Ensure the directory part of all file names exists."""
    for fname in fnames:
        directory = os.path.split(fname)[0]
        if not os.path.exists(directory):
            os.makedirs(directory)


