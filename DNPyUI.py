from PyQt5 import uic
from PyQt5 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
import sys
from functions import *
import subprocess


baseUIClass, baseUIWidget = uic.loadUiType("ui/DNPyUI.ui")


class AppWindow(baseUIWidget, baseUIClass):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.setupUi(self)
        self.pathButton.clicked.connect(self.open_exp_path)
        self.toolButton.clicked.connect(self.open_exp_powers)
        self.cancelButton.clicked.connect(self.close)
        self.startButton.clicked.connect(self.start)
        sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)



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
        self.generalGroup.setEnabled(False)
        self.dnpGroup.setEnabled(False)
        self.t1SeriesEval.setEnabled(False)
        self.makeFigs.setEnabled(False)
        self.resultsBrowser.setEnabled(True)
        self.startButton.setEnabled(False)
        self.cancelButton.setText("Abort")
        self.dnpThread = Thread(path, kwargs)
        self.dnpThread.dataReady.connect(self.onDataReady)
        self.dnpThread.start()

    def normalOutputWritten(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        self.resultsBrowser.append(text)



    def open_exp_path(self):
        self.path.setText(str(QFileDialog.getExistingDirectory(self, "Select Experiment Directory")))
    def open_exp_powers(self):
        self.powerFile.setText(str(QFileDialog.getOpenFileName(self, "Select Powers CSV File","","csv files (*.csv)")[0]))

    def onDataReady(self, exps):
        self.exps = exps
        self.cancelButton.setText("Close")
        sys.stdout = sys.__stdout__
        print(self.resultsBrowser.toPlainText())

class Thread(QtCore.QThread):
    dataReady = QtCore.pyqtSignal(list)
    def __init__(self, path, kwargs):
        super().__init__()
        self.path = path
        self.kwargs = kwargs

    def stop(self):
        self._flag = False

    @QtCore.pyqtSlot(list)
    def run(self):
        self._flag = True
        exps = return_exps(self.path, **self.kwargs)
        self.dataReady.emit(exps)



class EmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))

        
def main():
    app = QtWidgets.QApplication(sys.argv)
    w = AppWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
