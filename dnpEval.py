from typing import List
import os
from functions import return_exps, process_exps
from functions import NMRData

foldersInDir = []

# For evaluating several folders at once
path = '/run/media/hadi/0ea14823-f28b-4ca3-ab38-748260d2f31e/Work/RUB/HADI/Data/ODNP/Tests/20180802_CWODNP_TempolTest4'
# Comment the nex line if you only want to evaluate one set of data. This is for evaluating several folders at the same time
# foldersInDir = [os.path.join(path,d) for d in os.listdir(path) if os.path.isdir(os.path.join(path,d))]

powerFile = ''
kwargs = {
    # 'PCreal' for real phase cycled channel or 'PCmagn' or 'real' or 'magn'
    't1Calc': 'PCreal',
    'phase': 'first',  # all, none, first
    'ftWindow': 200,  # X domain for FT plot, int. and phase calculation
    # Does not allow peaks out of this domain to be calculated in int.
    'maxWin': 1000,
    'lS': 60,  # Left shift points or 'auto'
    'lSt1': 70,  # LS for t1 experiments (this differs from DNP sometimes)
    'rS': 110,  # Right shidt points
    'rSt1': 75,
    'lB': 10,  # Line boradening [Hz]
    'offCor': True,  # offset correctioon
    'basCor': True,  # baseline correction
    'evalPath': 'evalTest',
    'plotDpi': 100,  # plot file resolution
    'plotExts': ['jpg'],  # remove all if you do not want plots to be saved
    'make3dPlots': True,
    'process': True,
    'debug': False,
    'dumpToCsv': True,
    'figSize': (13, 8),
    'powerFile': powerFile,
    't1SeriesEval': True,
    't1SeriesPolDeg': 1,  # Polynomial degree for T1 series fit (default = 1)
    # The tolerance for T1 experiment error, this is to prevent bad T1 series fits
    't1ErrorTol': 15,
    'kSigmaCalc': True,
    'enhCalc': True,
}

if __name__ == '__main__':
    if foldersInDir:
        for path in foldersInDir:
            exps = return_exps(path, **kwargs)  # type: List[NMRData]
            if kwargs['process']:
                results = process_exps(exps, **kwargs)
    else:
        exps = return_exps(path, **kwargs)  # type: List[NMRData]
        if kwargs['process']:
            results = process_exps(exps, **kwargs)
