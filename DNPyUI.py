from PyQt5 import uic
from PyQt5.QtWidgets import *
import io
import sys
from functions import *
import threading
import time
from contextlib import redirect_stdout
import subprocess


baseUIClass, baseUIWidget = uic.loadUiType("ui/DNPyUI.ui")

class AppWindow(baseUIWidget, baseUIClass):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.setupUi(self)
        self.pathButton.clicked.connect(self.open_exp_path)
        self.toolButton.clicked.connect(self.open_exp_powers)
        self.buttonBox.accepted.connect(self.start)

    def start(self):
        path = self.path.text()
        powerFile = self.powerFile.text()
        phase = ''
        if self.apCombo.currentText() == "First experiment":
            phase = 'first'
        elif self.apCombo.currentText() == "All experimn":
            phase = 'all'
        elif self.apCombo.currentText() == "No auto phasing":
            phase = 'none'
        ftWindow = self.ftWindow.value()
        lB = self.lB.value()
        maxWin = self.maxWin.value()
        evalPath = self.evalPath.text()
        basCor = self.basCor.isChecked()
        offCor = self.offCor.isChecked()
        dumpToCsv = self.dumpToCsv.isChecked()
        kSigmaCalc = self.kSigmaCalc.isChecked()
        debug = self.Debug.isChecked()
        lS = self.lS.value()
        rS = self.rS.value()
        t1SeriesEval = self.t1SeriesEval.isChecked()
        t1Calc = ''
        if self.t1CalcCombo.currentText() == 'Phase cycled real':
            t1Calc = 'PCreal'
        if self.t1CalcCombo.currentText() == 'Phase cycled magnitude':
            t1Calc = 'PCmagn'
        if self.t1CalcCombo.currentText() == 'Real':
            t1Calc = 'real'
        if self.t1CalcCombo.currentText() == 'Magnitude':
            t1Calc = 'magn'
        lSt1 = self.lSt1.value()
        rSt1 = self.rSt1.value()
        t1SeriesPolDeg = self.t1SeriesPolDeg.value()
        makeFigs = self.makeFigs.isChecked()
        plotExts = []
        if makeFigs:
            for index, val in enumerate(['jpg','png','pdf','eps']):
                if getattr(self,val).isChecked():
                    plotExts.append(val)
        figWidth = self.figWidth.value()
        figHeight = self.figHeight.value()
        figDpi = self.figDpi.value()

        kwargs = {
            't1Calc': t1Calc,
            'phase': phase,
            'ftWindow': ftWindow,
            'maxWin': maxWin,
            'lS': lS,
            'lSt1': lSt1,
            'rS': rS,
            'rSt1': rSt1,
            'lB': lB,
            'offCor': offCor,
            'basCor': basCor,
            'evalPath': evalPath,
            'plotDpi': figDpi,
            'plotExts': plotExts,
            'process': True,
            'debug': debug,
            'dumpToCsv': dumpToCsv,
            'figSize': (figWidth, figHeight),
            'powerFile': powerFile,
            't1SeriesEval': t1SeriesEval,
            't1SeriesPolDeg': t1SeriesPolDeg,
            'kSigmaCalc': kSigmaCalc,
        }
        stdout = sys.stdout
        sys.stdout = io.StringIO()
        thread = threading.Thread(target=return_exps, args=(path,), kwargs=kwargs)
        thread.daemon = False  # Daemonize thread
        thread.start()  # Start the execution
        time.sleep(5)
        output = sys.stdout.getvalue()
        sys.stdout = stdout
        print(output)
        # stdout = sys.stdout
        # sys.stdout = io.StringIO()
        # exps = return_exps(path, **kwargs)
        # output = sys.stdout.getvalue()
        # sys.stdout = stdout
        # print(output)

    def open_exp_path(self):
        self.path.setText(str(QFileDialog.getExistingDirectory(self, "Select Experiment Directory")))
    def open_exp_powers(self):
        self.powerFile.setText(str(QFileDialog.getOpenFileName(self, "Select Powers CSV File","","csv files (*.csv)")[0]))



def main():
    app = QApplication(sys.argv)
    w = AppWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
