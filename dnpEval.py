from typing import List

from functions import *
from functions import NMRData

path = '/mnt/EPR-HADI/HADI/Data/ODNP/Tests/20180802_CWODNP_TempolTest4'
powerFile = ''
kwargs = {
    't1Calc': 'PCreal',  # 'PCreal' for real phase cycled channel or 'PCmagn' or 'real' or 'magn'
    'phase': 'first',  # all, none, first
    'ftWindow': 200,  # X domain for FT plot, int. and phase calculation
    'maxWin': 1000,  # Does not allow peaks out of this domain to be calculated in int.
    'lS': 70,  # Left shift points or 'auto'
    'lSt1': 70,  # LS for t1 experiments (this differs from DNP sometimes)
    'rS': 110,  # Right shidt points
    'rSt1': 75,
    'lB': 2,  # Line boradening [Hz]
    'offCor': True,  # offset correctioon
    'basCor': True,  # baseline correction
    'evalPath': 'evalHadiTest',
    'plotDpi': 250,  # plot file resolution
    'plotExts': [],  # remove all if you do not want plots to be saved
    'process': True,
    'debug': False,
    'dumpToCsv': True,
    'figSize': (13, 8),
    'powerFile': powerFile,
    't1SeriesEval': True,
    't1SeriesPolDeg': 1,  # Polynomial degree for T1 series fit (default = 1)
    'kSigmaCalc': True,
    }

if __name__ == '__main__':
    exps = return_exps(path, **kwargs)  # type: List[NMRData]
