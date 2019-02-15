from typing import List
import os
from functions import return_exps
from functions import NMRData

foldersInDir = []

# For evaluating several folders at once
path = '/home/hadi/Desktop/20190209_PhaseCyclingTest'
# Comment the nex line if you only want to evaluate one set of data. This is for evaluating several folders at the same time
# foldersInDir = [os.path.join(path,d) for d in os.listdir(path) if os.path.isdir(os.path.join(path,d))]

powerFile = ''
kwargs = {
    't1Calc': 'PCreal',  # 'PCreal' for real phase cycled channel or 'PCmagn' or 'real' or 'magn'
    'phase': 'first',  # all, none, first
    'ftWindow': 200,  # X domain for FT plot, int. and phase calculation
    'maxWin': 1000,  # Does not allow peaks out of this domain to be calculated in int.
    'lS': 60,  # Left shift points or 'auto'
    'lSt1': 70,  # LS for t1 experiments (this differs from DNP sometimes)
    'rS': 110,  # Right shidt points
    'rSt1': 75,
    'lB': 5,  # Line boradening [Hz]
    'offCor': True,  # offset correctioon
    'basCor': True,  # baseline correction
    'evalPath': 'evalTestPhaseCorrection',
    'plotDpi': 100,  # plot file resolution
    'plotExts': ['jpg'],  # remove all if you do not want plots to be saved
    'make3dPlots': True,
    'process': True,
    'debug': True,
    'dumpToCsv': True,
    'figSize': (13, 8),
    'powerFile': powerFile,
    't1SeriesEval': True,
    't1SeriesPolDeg': 1,  # Polynomial degree for T1 series fit (default = 1)
    't1ErrorTol': 1, # The tolerance for T1 experiment error, this is to prevent bad T1 series fits
    'kSigmaCalc': True,
    'enhCalc': True,
    }

if __name__ == '__main__':
    if foldersInDir:
        for path in foldersInDir:
            exps = return_exps(path, **kwargs)  # type: List[NMRData]
    else:
        exps = return_exps(path, **kwargs)  # type: List[NMRData]
