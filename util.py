from __future__ import division, print_function, unicode_literals
from collections import OrderedDict
from subprocess import Popen, PIPE
import tempfile
import StringIO
import itertools
import re

import matplotlib as mpl
import numpy as np
import pandas


_NS_REGEX = re.compile(r'(\d+)', re.U)


def set_postmortem_hook():
    import sys, traceback, ipdb
    def _excepthook(exc_type, value, tb):
        traceback.print_exception(exc_type, value, tb)
        print()
        ipdb.pm()
    sys.excepthook = _excepthook


def natsort_key(s):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


def by_chrom(func, *tables, **kwargs):
    """
    Split one or more dataframes by chromosome.
    Apply function `func()` to each chromosome "part".
    Yield results.
    
    Input
    -----
    func: function to apply to split dataframes.
        The expected signature is `func(chrom, df1[, df2[, ...])`, 
        where `df1, df2, ...` are subsets of the input dataframes.
        The function can return anything.
        
    tables: sequence of BED-like `pandas.DataFrame`s.
        The first column of each dataframe must be chromosome labels.
    
    chroms: sequence of str, optional
        Select which chromosome subsets of the data to apply the function to.
        Defaults to all unique chromosome labels in the first dataframe input,
        in natural sorted order.
        
    ret_chrom: bool, optional (default: False)
        Yield "chromosome, value" pairs as output instead of only values.
        
    parallel: bool, optional (default: False)
        Use IPython to fork out execution on the individual parts.
    
    Returns
    -------
    Generator that yields the output of running `func` on each chromosome
        -or-
    Future object from IPython's `map_async`.
    
    """
    chroms = kwargs.setdefault('chroms', None)
    parallel = kwargs.setdefault('parallel', False)
    ret_chrom = kwargs.setdefault('ret_chrom', False)

    if chroms is None:
        chrom_field = tables[0].columns[0]
        chroms = natsorted(tables[0][chrom_field].unique())

    grouped_tables = [table.groupby(table.columns[0]) for table in tables]

    def iter_partials(): 
        for chrom in chroms:
            partials = []
            for gby in grouped_tables:
                try:
                    partials.append(gby.get_group(chrom))
                except KeyError:
                    partials.append(gby.head()[0:0])
            yield partials
    
    if parallel:
        from IPython.parallel import Client
        pool = Client()[:]
        imap = pool.map_async
    else:
        imap = itertools.imap

    if ret_chrom:
        def run_job(chrom, partials):
            return chrom, func(chrom, *partials)
    else:
        def run_job(chrom, partials):
            return func(chrom, *partials)
    
    return imap(run_job, chroms, iter_partials())


def chrom_sorted(df, sort_by=None, reset_index=True):
    if sort_by is None:
        return pandas.concat(
            by_chrom(lambda c,x:x, df),
            axis=0,
            ignore_index=reset_index
        )
    else:
        return pandas.concat(
            by_chrom(lambda c,x:x, df.sort(sort_by)),
            axis=0,
            ignore_index=reset_index
        )
    

# To use pandas DataFrames with bedtools...
def dftemp(df, **kwargs):
    """
    Write pandas DataFrame to a temporary csv file.
    Useful in a `with` block (file is deleted at context teardown).

    >>> with dftemp(my_df) as f:
            ... command that requires a csv file ...

    """
    fh = tempfile.NamedTemporaryFile()
    df.to_csv(fh, sep=b'\t', index=False, header=False, **kwargs)
    fh.flush()  # DON'T FORGET TO FLUSH!!!
    return fh


def run(cmd, outnames=None):
    """
    Run an external command in a subprocess, capturing stdout.
    We assume the command prints a CSV file and load it into a dataframe.

    """
    print(' '.join(cmd))
    
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise IOError("process failed: %d\n%s\n%s" % (p.returncode, out, err))
    
    f_out = StringIO.StringIO(out)
    return pandas.read_csv(f_out, sep=b'\t', names=outnames)

