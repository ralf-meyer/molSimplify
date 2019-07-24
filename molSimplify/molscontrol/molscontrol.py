"""
molscontrol.py
An automonous on-the-fly job control system for DFT geometry optimization aided by machine learning techniques.

Written by Chenru Duan, Kulik Group at MIT.
crduan@mit.edu, hjkulik@mit.edu
"""

import sys
from .io_tools import get_configure
from .dynamic_classifier import dft_control
import argparse


def main():
    """
    The main function of job controls.

    Parameters
    ----------
    pid: the process (DFT geometry optimization).

    Returns
    -------
    None
    """
    try:
        pid = sys.argv[1]
    except:
        pid = False
        print("NO PID to control. Should be in a test mode.")
    kwargs = get_configure()
    kwargs.update({"pid": pid})
    dftjob = dft_control(**kwargs)
    stop = False
    while not stop:
        stop = dftjob.update_and_predict()


if __name__ == "__main__":
    main()
