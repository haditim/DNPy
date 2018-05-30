from functions import *

path = '/mnt/EPR-HADI/HADI/Data/ODNP/Tests/20180528_H2O_07TubeNew_NewLongFixedCoil'
kwargs = {
    't1Calc': 'PCreal',
    'phase': 'first',  # all, none, first
    'ftWindow': 200,  # X domain for FT plot, int. and phase calculation
    'maxWin': 1000,  # Does not allow peaks out of this domain to be calculated in int.
    'lS': 30,  # Left shift points or 'auto'
    'lSt1': 70,  # LS for t1 experiments (this differs from DNP sometimes)
    'rS': 110,  # Right shidt points
    'rSt1': 75,
    'lB': 2,  # Line boradening [Hz]
    'offCor': True,  # offset correctioon
    'basCor': True,  # baseline correction
    'evalPath': 'evalHadiTest',
    'plotDpi': 250,  # plot file resolution
    'plotExts': ['jpg'],  # remove all if you do not want plots to be saved
    'process': True,
    'debug': False,
    'dumpToCsv': True,
    'figSize': (13, 8),
    }
exps = return_exps(path, **kwargs)
