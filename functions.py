# DNPy functions file
import datetime
import glob
import os
import traceback
import concurrent.futures
import time

import numpy as np
import scipy as sp
import matplotlib
# to prevent Tcl_AsyncDelete: async handler deleted by the wrong thread
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import fft
from scipy.fftpack import fftshift
from scipy.optimize import curve_fit
from scipy import stats
from scipy.interpolate import interp1d
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib import cm
import scipy.io as sio
style.use('seaborn-whitegrid')


class NMRData(object):
    """
    NMRData objects is used to import and store NMRData from several file formats.
    Usage:
    > NMRData = NMRData(path, type)
    supported type: TopSpin (Magritek, ntnmr, and spinSight support can be added by using \
    pyNMR package
    """

    def __init__(
            self,
            path,
            datatype,
            container=0,
            sizeTD1=0,
            ft_only=[],
            hiper_skip_footer=0,
            hiper_skip_header=3,
            endianess="<",
            autoPhase=True,
            ph=0,
            **kwargs
    ):
        """ This reads the data """
        process = kwargs.get('process', True)
        self.debug = kwargs.get('debug', False)
        self.lS = kwargs.get('lS', 'auto')
        self.rS = kwargs.get('rS', 100)
        self.lB = kwargs.get('lB', 2)
        self.lSt1 = kwargs.get('lSt1', self.lS)
        self.rSt1 = kwargs.get('rSt1', self.rS)
        self.basCor = kwargs.get('basCor', True)
        self.offCor = kwargs.get('offCor', True)
        self.autoPhase = autoPhase
        self.ph = ph
        self.ftWindow = kwargs.get('ftWindow', 200)
        self.maxWin = kwargs.get('maxWin', 1000)
        self.t1Calc = kwargs.get('t1Calc', 'PCmagn')
        self.t1ErrorTol = kwargs.get('t1ErrorTol', 0.5)
        if datatype == '':
            print("No Datatype - Setting it to ntnmr")
            datatype = "ntnmr"
        self.carrier = 0
        self.allFid = []
        self.allFid.append([])
        self.sizeTD1 = sizeTD1
        self.title = ['no title']
        self.vdList = []
        self.fidTimeHistory = {}
        self.parDictionary = {}
        self.expType = ''
        self.dbSet = 0
        self.powerDbm = 0
        self.powerMw = 0
        self.expTime = 0
        self.maxFreq = 0  # type: int
        endianness = kwargs.get('endianness', np.dtype('<i4'))
        if datatype == 'TopSpin':
            # The acqus file containts the spectral width SW_h and 2*SizeTD2 as ##$TD
            # The acqu2s file contains TD1 as ##$TD
            directory = os.path.abspath(path)
            acqusFile = open(os.path.join(directory, "acqus"), mode='r')
            self.expTime = os.path.getmtime(
                os.path.join(os.path.join(directory, "acqus")))  # time of the experiment
            # check if acqu2sfile exists, if yes, experiment is 2D!
            if os.path.isfile(os.path.join(directory, "acqu2s")):
                acqu2sFile = open(os.path.join(directory, "acqu2s"), mode='r')
                acqu2File = open(os.path.join(directory, "acqu2"), mode='r')
                self.is2D = True
            else:
                self.is2D = False
                self.sizeTD1 = 1
            # this could be crafted into a common routine which gives names of parameters
            # parameters and works the same for e.g., spinsight and topspin
            count = 0
            while True:
                # try:
                count += 1
                line = acqusFile.readline().strip()
                if "=" in line:
                    line = line.split("=")
                elif len(line) > 0:
                    line = line.split(" ")
                elif len(line) == 0 or count > 1000:
                    if self.debug:
                        print("Ended reading acqus file at line ", count)
                    break
                else:
                    next
                if line[0] == "##$SW_h":
                    # this line might be chopping the last digit off....
                    # self.sweepWidthTD2 = int(float(line[1][:-1]))
                    self.sweepWidthTD2 = int(float(line[1]))
                    if self.debug:
                        print("SweepWidthTD2: ", self.sweepWidthTD2)
                elif line[0] == "##$TD":
                    self.sizeTD2 = int(line[1]) / 2
                    if self.debug:
                        print("sizeTD2: ", self.sizeTD2)
                elif line[0] == "##$SFO1":
                    self.carrier = float(line[1]) * 1e6
                    if self.debug:
                        print("SFO1:", self.carrier)
                elif len(line) == 0:
                    break
                if len(line[0]) > 1:
                    if "@" in line[-1]:
                        # this line contains date, time, some unknown stuff and user, does not
                        # work with all bruker files, hence try only"
                        try:
                            self.parDictionary["date"] = line[1].strip()
                            self.parDictionary["time"] = line[2].strip()
                        except Exception as e:
                            if self.debug:
                                print(
                                    'error {} happened setting date of experiment'.format(e))
                            pass
                    elif line[0] == "##$D":
                        delays1 = acqusFile.readline().strip()
                        delays2 = acqusFile.readline().strip()
                        self.parDictionary["d"] = [float(d) for d in delays1.strip().split(" ")] + [
                            float(d) for d in delays2.strip().split(" ")]
                    elif line[0] == "##$L":
                        loopCounters = acqusFile.readline().strip()
                        self.parDictionary["l"] = [float(l)
                                                   for l in loopCounters.strip().split(" ")]
                    # We don't need to store all this
                    # else:
                    # self.parDictionary[line[0][2:].strip()] = line[1].strip()

            if self.is2D:
                count = 0
                while True:
                    # try:
                    count += 1
                    line = acqu2sFile.readline().strip()
                    if "=" in line:
                        line = line.split("=")
                    elif len(line) == 0 or count > 1000:
                        if self.debug:
                            print("Ended reading acqu2s file at line ", count)
                        break
                    else:
                        next
                    if line[0] == "##$TD" and self.sizeTD1 == 0:
                        self.sizeTD1 = int(line[1])
                        if self.debug:
                            print("sizeTD1: ", self.sizeTD1)
                    elif len(line) == 0:
                        break
                if os.path.isfile(directory + "/vdlist"):
                    if self.debug:
                        print("VD File exists!")
                    with open(directory + "/vdlist", 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            line = line.strip()
                            if line.endswith('m'):
                                self.vdList.append(float(line[:-1]) / 1000.)
                            else:
                                self.vdList.append(float(line))

            if self.debug:
                print("TD2: ", self.sizeTD2)
                print("TD1: ", self.sizeTD1)
                print("Carrier:", self.carrier)

            if self.is2D:
                self.f = open(os.path.join(path, "ser"), mode='rb')
            else:
                self.f = open(os.path.join(path, "fid"), mode='rb')

            dataString = np.frombuffer(self.f.read(), dtype=endianess + "i4")
            if self.debug:
                print("len(dataString) new: ", len(dataString))

            self.data = dataString
            self.sizeTD2 = len(self.data) / self.sizeTD1 / 2

            dwellTime = 1. / self.sweepWidthTD2
            self.fidTime = np.linspace(
                0, (self.sizeTD2 - 1) * dwellTime, num=self.sizeTD2)

            # here we create one array of complex numbers for each of the FIDs
            # i runs over all fids in a ser file, in case of a fid file i = 0
            # TD1 is number of FIDs, TD2 is number of datapoints in each FID
            for i in range(0, self.sizeTD1):
                realPart = self.data[int(
                    i * self.sizeTD2 * 2):int((i + 1) * self.sizeTD2 * 2):2]
                imagPart = sp.multiply(
                    self.data[int(i * self.sizeTD2 * 2 + 1):
                              int((i + 1) * self.sizeTD2 * 2 + 1):2], 1j)
                self.allFid[0].append(sp.add(realPart, imagPart))

            # here we read the experiment title (only the one stored in pdata/1):
            # could be made to read titles from all pdata folders (if needed)
            try:
                pathToTitle = os.path.join(directory, 'pdata', '1', 'title')
                titleFile = open(pathToTitle, mode='r')
                title = list(titleFile)
                self.title = ''.join([line.strip() for line in title])
                if 'DNP' in self.title or 'dnp' in self.title or 'baseline for integration' in self.title:
                    self.expType = 'dnp'
                elif 'T1' in self.title or '$T_1$' in self.title or 'T_{1,0}' in self.title or 't1' in self.title:
                    self.expType = 't1'
            except Exception as e:
                if self.debug:
                    print('{} error occured reading title file.'.format(e))
                else:
                    pass
            self.dwellTime = dwellTime
            self.expNum = os.path.basename(directory)
            self.directory = directory
        self.calculate_phase_cycles()
        # self.calculate_exp_integral
        # HADI: I want to keep track of time after processing
        self.fidTimeHistory['orig'] = self.fidTime
        self.power_calc(**kwargs)

    def power_calc(self, **kwargs):
        # Calculate powers
        if 'dBm' in self.title:
            self.powerDbm = float(self.title.split(" ")[-2])
            self.dbSet = None
            self.powerMw = 10.0 ** ((self.powerDbm) / 10.0)
        elif 'dB' in self.title:
            self.dbSet = float(self.title[:-2].split("set ", 1)[1])

        # TODO: I removed this which is for the Han's lab for now
        # elif 'dnpPowerMatFile' in kwargs and self.expType == 'dnp':
        #         if dnpCounter == 0:
        #             dnpFirstTime = self.expTime - \
        #                 float(dnpPowerMatFile[0][0])
        #         if kw['debug']:
        #             print('exp time: ', datetime.datetime.fromtimestamp(self.expTime),
        #                   ' time diff: ', datetime.datetime.fromtimestamp(
        #                       self.expTime - dnpFirstTime),
        #                   ' dnpCounter: ', dnpCounter,
        #                   ' dir: ', name)
        #         timeInd, time = find_nearest(
        #             dnpPowerMatFile[0], self.expTime - dnpFirstTime)
        #         self.powerDbm = float(dnpPowerMatFile[1][timeInd])
        #         self.powerMw = 10.0 ** ((self.powerDbm + 24.445) / 10.0)
        #     elif 't1PowerMatFile' in kwargs and self.expType == 't1':
        #         if t1Counter == 0:
        #             t1FirstTime = self.expTime - float(t1PowerMatFile[0][0])
        #         self.powerDbm = float(
        #             t1PowerMatFile[1][find_nearest(t1PowerMatFile[0], self.expTime - t1FirstTime)[0]])
        #         self.powerMw = 10.0 ** ((self.powerDbm + 24.445) / 10.0)

        else:
            try:
                self.dbSet = float(self.title.split()[-1])
                self.powerDbm = kwargs['powers_dict']
            except Exception:
                if kw['debug']:
                    print("Experiment %s does not have dB or dBm set in title.\
                        So I use 60dB attenuation" % self.title)
                self.dbSet = 60.
                self.powerDbm = -100.

    def offset_correction(self, fromPos, toPos, fraction=0.75, whichFid=0):
        self.check_to_pos(toPos)
        lenFid = len(self.allFid[fromPos][0])
        startOffset = int(fraction * lenFid)
        if self.debug:
            # print("fid: ", whichFid)
            print("startOffset: ", startOffset)
            print("len allFid: ", len(self.allFid))
            print("len allFid[0]: ", len(self.allFid[0]))
            print("===============================================================")
            print("len allFid[0][" + str(whichFid) + "]: ",
                  len(self.allFid[0][0]))
        oReal = np.mean(np.real(self.allFid[fromPos][whichFid][startOffset:]))
        stdReal = np.std(
            np.real(self.allFid[fromPos][whichFid][startOffset:]) - oReal)
        if self.debug:
            print("offsetReal: ", oReal)
            print("stdReal: ", stdReal)
        oImag = np.mean(np.imag(self.allFid[fromPos][whichFid][startOffset:]))
        stdImag = np.std(
            np.imag(self.allFid[fromPos][whichFid][startOffset:]) - oImag)
        if self.debug:
            print("offsetImag: ", oImag)
            print("stdImag: ", stdImag)
        self.allFid[toPos] = [np.real(fid) - oReal + 1j * (np.imag(fid) - oImag)
                              for fid in self.allFid[fromPos]]

    def line_broadening(self, fromPos, toPos, LB):
        """Applies line broadening of given width (in Hz) to the FID. Resulting spectrum
        (after fft is called) will be convolution of the original spectrum (fromPos)
        and Lorentzian with full-width-at-half-maximum equal to LB"""

        self.check_to_pos(toPos)
        self.allFid[toPos] = sp.multiply(
            self.allFid[fromPos][:], sp.exp(-self.fidTime * LB * np.pi))

    def fourier_transform(self, fromPos, toPos, only=[]):
        self.check_to_pos(toPos)
        if len(only) > 0:
            self.allFid[toPos] = np.array(
                [fftshift(fft(self.allFid[fromPos][fidIndex])) for fidIndex in only])
        else:
            self.allFid[toPos] = np.array(
                [fftshift(fft(fid)) for fid in self.allFid[fromPos]])

        self.frequency = np.linspace(
            -self.sweepWidthTD2 / 2, self.sweepWidthTD2 / 2, len(self.allFid[fromPos][0]))

    def baseline_correction(self, fromPos, toPos, fitRange, scale="Hz", order=1, applyLocally=False):
        """Polynomial base line correction.
        - fromPos: position where current spectra are stored
        - toPos: position where corrected spectra will be stored
        - fitRange: a list of lists, with each list element giving the limits of a baseline section.
        - order: polynomial order of the baseline correction. Use order 1 for linear correction
        - applyLocally: set to True for an only local baseline correction that extends only over
        the entire fitRange."""

        if len(self.allFid) > toPos:
            self.allFid[toPos] = []
        # if toPos does not exist, create it
        else:
            self.allFid.append([])
        # check that now the toPos exists:
        assert len(self.allFid) > toPos, "toPos too high"
        for k in range(self.sizeTD1):
            xVals = []
            yVals = []
            indices = []
            for pair in fitRange:
                i1, i2 = self.get_indices(pair, scale=scale)
                indices.extend([i1, i2])
                xVals.extend(self.frequency[i1:i2])
                yVals.extend(self.allFid[fromPos][k][i1:i2])
            z = np.polyfit(xVals, yVals, order)
            p = np.poly1d(z)
            if applyLocally:
                self.allFid[toPos].append(self.allFid[fromPos][k])
                self.allFid[fromPos][k][min(indices):max(indices)] -= p(
                    self.frequency[min(indices):max(indices)])
            else:
                self.allFid[toPos].append(
                    self.allFid[fromPos][k] - p(self.frequency))

    def baseline_correction_mean(self, fromPos, toPos, totalPoints):
        """HADI: baseline_correction does not work on our data. this is a simple mean correction."""

        self.check_to_pos(toPos)
        self.allFid[toPos] = [
            k - np.mean(k[len(k) - totalPoints:]) for k in self.allFid[fromPos]]

    def phase(self, fromPos, toPos, phase, applyIndex='all', degree=True):
        self.check_to_pos(toPos)
        self.fill_to_position(fromPos, toPos)
        if degree:
            phaseFactor = np.exp(-1j * float(phase) / 180. * np.pi)
        else:
            phaseFactor = np.exp(-1j * phase)
        if applyIndex == 'all':
            self.allFid[toPos] = [
                fid * phaseFactor for fid in self.allFid[fromPos]]
        else:
            applyIndex = [applyIndex] if not type(
                applyIndex).__name__ == 'list' else applyIndex
            for i in applyIndex:
                self.allFid[toPos][i] = self.allFid[fromPos][i] * phaseFactor

    def auto_phase(self, fromPos, toPos, index, start, stop, scale="Hz"):
        """This function should get fromPos and index pointing to a spectrum.
        It returns the phase for minimizing the integral over the imag. part
        in the spectral region of interest, in degrees. HADI: It puts the result of
        phasing to toPos"""

        i1, i2 = self.get_indices([start, stop], scale=scale)
        phiTest = np.linspace(0, 359, num=360)
        integrals = np.zeros(np.size(phiTest))
        for k in range(len(integrals)):
            integrals[k] = np.sum(np.imag(
                self.allFid[fromPos][index][i1:i2] * np.exp(-1j * float(phiTest[k]) / 180. * np.pi)))
        magnMin = np.argmin(integrals) + \
            (np.argmax(integrals) - np.argmin(integrals))//2
        self.phase(fromPos, toPos, phiTest[magnMin])
        return phiTest[magnMin]

    def auto_phase0(self, fromPos, toPos, phaseIndex, start, stop, applyIndex='all', scale="Hz"):
        """This function should get fromPos and phaseIndex pointing to a spectrum.
        It returns the phase for maximimizing the integral over the real part
        in the spectral region of interest, in degrees. 
        HADI: It puts the result of phasing to toPos with the applyIndex indices"""
        i1, i2 = self.get_indices([start, stop], scale=scale)
        phiTest = np.linspace(0, 359, num=360)
        integrals = np.zeros(np.size(phiTest))
        for k in range(len(integrals)):
            integrals[k] = sp.integrate.cumtrapz(np.real(
                self.allFid[fromPos][phaseIndex][i1:i2] * np.exp(-1j * float(phiTest[k]) / 180. * np.pi)))[-1]
        self.phase(fromPos, toPos, phiTest[np.argmax(
            integrals)], applyIndex=applyIndex)
        return phiTest[np.argmax(integrals)]

    def fwhm(self, fromPos, index, start, stop, scale="Hz"):
        i1, i2 = self.get_indices([start, stop], scale=scale)
        x = self.frequency[i1:i2]
        y = np.real(self.allFid[fromPos][index][i1:i2])
        linewidth = fwhm.fwhm(x, y)
        if self.debug:
            print(("Linewidth: {} Hz".format(linewidth)))
        return linewidth

    def phase_first_order(self, fromPos, toPos, phi1, degree=True,
                          noChangeOnResonance=False, pivot=0, scale="Hz"):
        """This should only be applied to Fourier transformed spectral data
        It will lead to a zero phase shift at the upper end of the spectrum and to a
        phase shift of phi1 at the lower end, linear interpolation inbetween.
        this is the spinsight convention.
        If a pivot is provided the phase will not be changed at the pivot,
        but the total change accross the entire spectrum will still amount to phi1.
        """

        self.check_to_pos(toPos)
        phaseValues = np.linspace(0, phi1, num=len(self.frequency))
        if noChangeOnResonance:
            pivot = 0
        elif pivot != 0:
            print("Using pivot for first order phase correction")
            index = self.get_index(pivot, scale=scale)
            phaseValues = phaseValues - phaseValues[index]
        if degree:
            phaseValues = phaseValues * np.pi / 180
        self.allFid[toPos] = [spectrum * np.exp(-1j * phaseValues)
                              for spectrum in self.allFid[fromPos]]

    def left_shift(self, fromPos, toPos, shiftPoints):
        self.check_to_pos(toPos)
        self.allFid[toPos] = [self.allFid[fromPos][k][shiftPoints:]
                              for k in range(len(self.allFid[fromPos]))]
        self.fidTimeHistory['bLeftShift'] = self.fidTime
        self.fidTime = self.fidTime[:len(self.fidTime) - shiftPoints]

    def right_shift(self, fromPos, toPos, shiftPoints):
        self.check_to_pos(toPos)
        self.allFid[toPos] = [self.allFid[fromPos][k][:len(self.allFid[fromPos][k]) - shiftPoints]
                              for k in range(len(self.allFid[fromPos]))]
        self.fidTimeHistory['bRightShift'] = self.fidTime
        self.fidTime = self.fidTime[:len(self.fidTime) - shiftPoints]

    def zero_filling(self, fromPos, toPos, totalPoints):
        self.check_to_pos(toPos)
        z = np.zeros(totalPoints)
        self.allFid[toPos] = [np.append(k, z) for k in self.allFid[fromPos]]
        self.fidTimeHistory['bZeroFilling'] = self.fidTime
        self.fidTime = np.linspace(0, (len(self.fidTime)
                                       - 1 + totalPoints) * self.dwellTime, num=len(self.fidTime)
                                   + totalPoints)

    def get_joined_partial_spectra(self, fromPos, start, stop, scale="Hz", returnX=False):
        spectra = []
        for index in range(self.sizeTD1):
            spectra.extend(self.get_partial_spectrum(
                fromPos, index, start, stop, scale=scale))
        if returnX:
            x = np.array(list(range(len(spectra)))) * \
                self.sizeTD1 / float(len(spectra)) + 0.5
            return x, spectra
        else:
            return spectra

    def get_partial_spectrum(self, fromPos, index, start, stop, scale="Hz"):
        i1, i2 = self.get_indices([start, stop], scale=scale)
        return np.real(self.allFid[fromPos][index][i1:i2])

    def integrate(self, fromPos, index, start, stop, scale="Hz"):
        """This function integrates the real part between start and stop. standard scale is Hz
        Arguments:
        - `fromPos`:
        - `index`: index of the spectrum
        - `startFreq`: lower limit of integration
        - `stopFreq`: upper limit of integration
        - `scale`: Hz or ppm
        - `part`: real or magnitude
        """

        i1, i2 = self.get_indices([start, stop], scale=scale)
        # get the step in x-variable
        if scale == "Hz":
            step = np.abs(self.frequency[1] - self.frequency[0])
        if scale == "ppm":
            step = np.abs(self.ppmScale[1] - self.ppmScale[0])
        else:
            step = 1
        integralImag = sp.integrate.cumtrapz(np.imag(
            self.allFid[fromPos][index][i1:i2]), x=self.frequency[i1:i2], dx=step, initial=0)
        integralReal = sp.integrate.cumtrapz(np.real(
            self.allFid[fromPos][index][i1:i2]), x=self.frequency[i1:i2], dx=step, initial=0)
        integralMagn = sp.integrate.cumtrapz(np.abs(
            self.allFid[fromPos][index][i1:i2]), x=self.frequency[i1:i2], dx=step, initial=0)
        return{
            'x-axis': self.frequency[i1:i2],
            'integralReal': integralReal,
            'integralImag': integralImag,
            'integralMagn': integralMagn,
        }

    def get_peak(self, fromPos, index, start, stop, negative=False, scale="Hz"):
        """This function returns peak intensities in a given range;
        it searches for negative peaks if negative = True"""

        i1, i2 = self.get_indices([start, stop], scale=scale)
        spec = self.allFid[fromPos][index]
        if not negative:
            maxVal = np.max(np.real(spec[i1:i2]))
        else:
            maxVal = -np.max(-np.real(spec[i1:i2]))
        return maxVal

    def get_center_frequency(self, fromPos, index, start, stop, scale="Hz"):
        ind = np.where(abs(self.allFid[fromPos][index]) == max(abs(
            self.allFid[fromPos][index][self.get_index_from_frequency(start):self.get_index_from_frequency(stop)])))[0][
            0]
        frequency = self.frequency[ind]
        return ind, frequency

    def get_index_from_frequency(self, freq):
        return np.argmin(abs(self.frequency - freq))

    def get_index(self, value, scale="Hz"):
        if scale == "Hz":
            return self.get_index_from_frequency(value)
        elif scale == "ppm":
            return self.get_index_from_ppm(value)

    def get_indices(self, interval, scale="Hz"):
        i1 = self.get_index(interval[0], scale=scale)
        i2 = self.get_index(interval[1], scale=scale)
        if i1 > i2:
            return i2, i1
        else:
            return i1, i2

    def check_to_pos(self, toPos):
        if len(self.allFid) <= toPos:
            self.allFid.append([])

    def fill_to_position(self, fromPos, toPos):
        if len(self.allFid[toPos]) < len(self.allFid[fromPos]):
            self.allFid[toPos] = np.zeros((np.shape(self.allFid[fromPos])[0], np.shape(
                self.allFid[fromPos])[1]), dtype=np.complex64)

    def export(self, pos, count, filename, scale="Hz", xlim=[], complexType="r", fmt="%.3f"):
        if scale == "Hz":
            xData = self.frequency
        elif scale == "ppm":
            xData = self.ppmScale
        elif scale == "Time":
            L = len(self.allFid[pos][count])
            xData = np.linspace(0, float(L) / self.sweepWidthTD2, L)
        yDataR = np.real(self.allFid[pos][count])
        yDataI = np.imag(self.allFid[pos][count])
        yDataM = np.abs(self.allFid[pos][count])
        if xlim == []:
            start = 0
            stop = -1
        else:
            if scale == "Hz":
                start = self.get_index_from_frequency(xlim[0])
                stop = self.get_index_from_frequency(xlim[1])
            elif scale == "ppm":
                start = self.get_index_from_ppm(xlim[0])
                stop = self.get_index_from_ppm(xlim[1])
        if complexType == "r":
            data = list(zip(xData[start:stop], yDataR[start:stop]))
        np.savetxt(filename, data, fmt=fmt, delimiter="\t")

    def auto_phase1(self, fromPos, index, start=-1e6, stop=1e6, derivative=1,
                    penalty=1e3, scale='Hz'):
        """Automatic phase correction (0 + 1 order) based on entropy
        minimization (Chen et al: J. Mag. Res. 158, 164-168 (2002)).
        Minimizes entropy of phased spectrum + a penalty function (which is
        equal to integral of intensity**2 in regions where intensity<0 multiplied
        by the "penalty" parameter given in autoPhase input).
        Returns phase correction coefs in radians in array [ph0, ph1]
        which can be used by method phase01 to apply the phase correction.
        Derivative should be set to 1-4, increasing penalty puts more
        emphasis on non-negative spectrum.
        By default the spectrum in range +/-1MHz around offset is considered,
        the interval can be set using the start and stop which can be
        in either 'Hz' or 'ppm' scale"""

        assert start < stop, "start should be smaller than stop"
        assert penalty > 0, "penalty shoud be possitive"
        assert type(
            derivative) is int, "derivative should be a (small possitive) integer"
        assert derivative > 0, "need derivative > 0"
        spectrum = np.array(self.allFid[fromPos][index])
        # normalize the spectrum:
        spectrum = spectrum / np.abs(spectrum).sum()
        # zero everything that is out of start-stop frequency window
        if scale == 'Hz':
            for i in range(len(spectrum)):
                if self.frequency[i] < start:
                    spectrum[i] = 0
                if self.frequency[i] > stop:
                    spectrum[i] = 0
        if scale == 'ppm':
            for i in range(len(spectrum)):
                if self.ppmScale[i] < start:
                    spectrum[i] = 0
                if self.ppmScale[i] > stop:
                    spectrum[i] = 0
        # record initial values of penalty and entropy:
        penalty_start = self.__penalty(spectrum, penalty)
        entropy_start = self.__entropy(spectrum, derivative)
        # find the phase correction that minimizes the objective function
        correction = [0, 0]
        res = sp.optimize.minimize(self.__tryPhase, correction,
                                   args=(spectrum, derivative, penalty,))
        if self.debug:
            spectrum = self.__phase01(spectrum, res.x)
            print('penalty change:', self.__penalty(
                spectrum, penalty) - penalty_start)
            print('entropy change:', self.__entropy(
                spectrum, derivative) - entropy_start)
        return res.x

    def phase01(self, fromPos, toPos, correction):
        """apply zero and first order phase correction to spectra at fromPos
        and write the result to toPos. correction angles are in radians and
        are stored in array correction = [ph0, ph1]. first order
        correction leads to no change at first point of spectrum and maximum
        (ph1) change at the last point.
        This function can be used to apply the phase correction returned by
        auto_phase1"""

        self.check_to_pos(toPos)
        # here we apply the correction
        self.allFid[toPos] = [self.__phase01(spectrum, correction)
                              for spectrum in self.allFid[fromPos]]

    def __phase01(self, spectrum, correction):
        """Returns a spectrum (np.array) to which a specified phase correction
        was applied. ph0 and ph1 are in rad"""

        ph0, ph1 = correction[0], correction[1]
        phaseValues = np.linspace(0, ph1, num=len(spectrum)) + ph0
        corrections = np.exp(1j * phaseValues)
        return np.array([spectrum[i] * corrections[i] for i in range(len(spectrum))])

    def __entropy(self, spectrum, m):
        """Calculates get m-th derivative of the real part of spectrum and
        returns entropy of its absolute value. """

        assert type(m) is int, 'm should be a (possitive) integer'
        assert m > 0, 'need m > 0'
        # get the real part of the spectrum
        spect = np.array(spectrum)
        spect = spect.real
        # calculate the m-th derivative of the real part of the spectrum
        spectrumDerivative = spect
        for i in range(m):
            spectrumDerivative = np.gradient(spectrumDerivative)
        # now get the entropy of the abslolute value of the m-th derivative:
        entropy = sp.stats.entropy(np.abs(spectrumDerivative))
        return entropy

    def __penalty(self, spectrum, gamma):
        """return penalty function for the spectrum - sum of squares of
        all negative points in normalized spectrum multiplied by gamma"""

        penalty = 0
        # normalize the real part of the spectrum:
        spect = spectrum.real / np.abs(spectrum.real).sum()
        # calculate the penalty for the normalized real part
        for point in spect:
            if point < 0:
                penalty += point ** 2
        return penalty * gamma

    def __tryPhase(self, correction, spectrum, m, gamma):
        """Apply the phase correction to the spectrum, evaluate
        entropy and penalty of the resulting spectrum and return
        their sum (aka objective function)"""

        phased = self.__phase01(spectrum, correction)
        objective = self.__entropy(phased, m) + self.__penalty(phased, gamma)
        return objective

    def normalize_intensity(self, fromPos, toPos, scaling=1.0):
        """Takes spectrum in fromPos, divides the intensities by number of scans
        and receiver gain and writes the resulting spectrum to toPos. Useful for
        comparing intensities of specta taken with different settings. Optional
        parameter scaling can be used to correct for additional effects (e.g.
        mass). Spectrum will be divided by this factor - higher scaling means
        lower resulting intensity."""
        self.check_to_pos(toPos)
        # delete whatever is in toPos:
        self.allFid[toPos] = []
        # scaling factor:
        factor = float(self.parDictionary["$NS"]) * \
            float(self.parDictionary["$RG"]) * scaling
        for spectrum in self.allFid[fromPos]:
            self.allFid[toPos].append([point / factor for point in spectrum])

    def get_index_from_ppm(self, value):
        return np.argmin(abs(self.ppmScale - ppm))

    def calculate_phase_cycles(self):
        try:
            self.phaseCycles = len(self.allFid[0]) / len(self.vdList)
            self.vdListLen = len(self.vdList)
            if int(self.phaseCycles) != self.phaseCycles:
                raise ValueError
        except ValueError:
            print('Exp. {}: The number of FIDs ({}) and vdList ({}) do not match. '
                  'There should have been a problem in your experiment'.format(self.title, len(self.allFid[0]),
                                                                               len(self.vdList)))
            pass
        except Exception as e:
            if self.debug:
                print('No phase cycling channel found ({}).'.format(e))
            self.phaseCycles = 1
            self.vdListLen = 1

    def process(self):
        """ Process routine for DNP enhancement, k_sigma and T1 data. """
        if self.expType == 't1':
            self.left_shift(0, 1, self.lSt1)
            # remove the bruker filter, run only once, otherwise the data.fidTime gets
            # messed up and the following plot fails:
            self.right_shift(1, 1, self.rSt1)  # removing right handside zeros
        else:
            self.left_shift(0, 1, self.lS)
            self.right_shift(1, 1, self.rS)
        if self.basCor:  # try to remove Offset if asked to
            self.baseline_correction_mean(1, 1, len(self.fidTime) // 10)
        if self.offCor:
            self.offset_correction(1, 1)
        self.zero_filling(1, 2, int((2. * len(self.fidTime))))  # zero filling
        """ Next we apply some line broadening - we multiply the FID by an exponentially
        decaying function. This is not strictly necessary, but can be handy.
        We take the FID at position 1 apply the line broadening (2Hz) to it and write
        the result to position 2. """
        self.line_broadening(2, 3, self.lB)
        self.fourier_transform(3, 4)  # FT
        """ automatic zero-order phase correction to get the phase which gives
        maximum real amplitude of spectrum at pos 3 (index 0) in frequency interval """

    def calculate_exp_integral(self):
        # real and magnitude integrals
        # time(from vdList), integral for each channel
        self.real = [[] for i in range(0, int(self.phaseCycles + 1))]
        self.magn = [[] for i in range(0, int(self.phaseCycles + 1))]
        # power, phasecycled real, phasecycles magn.
        self.phC = [[] for i in range(0, 3)]
        try:
            self.fwhm = fwhm(self.frequency, np.real(self.allFid[5][0]))
        except Exception:
            self.fwhm = None
        for i in range(len(self.allFid[0])):  # calculate the integral
            # HADI: hack for T1 integration problems
            if (self.expType == 't1' and i == 0) or self.expType == 'dnp':
                self.maxValInd, self.maxFreq = self.get_center_frequency(
                    5, i, min(self.frequency), max(self.frequency))
            if not -self.maxWin < self.maxFreq < self.maxWin:
                self.maxFreq = 0
            integration = self.integrate(
                5,
                i,
                self.maxFreq - self.ftWindow,
                self.maxFreq + self.ftWindow,
                scale="Hz",
            )
            calIntReal = integration['integralReal'][-1]
            calIntImag = integration['integralImag'][-1]
            calIntMagn = integration['integralMagn'][-1]
            if self.debug:
                print("Center frequency for this exp.: {}".format(self.maxFreq))
            if i % (self.phaseCycles) == 0:  # appending time values
                if self.vdListLen > 1:
                    self.real[0].append(self.vdList[int(i / self.phaseCycles)])
                    self.magn[0].append(self.vdList[int(i / self.phaseCycles)])
                    self.phC[0].append(self.vdList[int(i / self.phaseCycles)])
                else:
                    self.real[0].append(0)
                    self.magn[0].append(0)
                    self.phC[0].append(0)
                # putting zeros for phase cycled channel
                self.phC[1].append(0)
                self.phC[2].append(0)
            # Appending integral values
            self.real[int(i % self.phaseCycles) + 1].append(calIntReal)
            self.magn[int(i % self.phaseCycles) + 1].append(calIntMagn)
        # phase cycling
        self.phC = np.asarray(self.phC)
        self.phC[1] = np.sum(self.real[1:], axis=0)
        self.phC[2] = np.sum(self.magn[1:], axis=0)
        if self.expType == 't1':
            if self.t1Calc == 'real':  # real, magn, PCreal, PCmagn
                self.fitData = [self.real[0], self.real[1]]
            elif self.t1Calc == 'magn':
                self.fitData = [self.magn[0], self.magn[1]]
            elif self.t1Calc == 'PCreal':
                # set to whatever you want to be fitted and plotted and source of T1
                self.fitData = [self.phC[0], self.phC[1]]
            elif self.t1Calc == 'PCmagn':
                # set to whatever you want to be fitted and plotted and source of T1
                self.fitData = [self.phC[0], self.phC[2]]
            self.t1fit = fit_t1(
                self.fitData[0], self.fitData[1], t1ErrorTol=self.t1ErrorTol)


def fwhm(x, y):
    """Calulate the FWHM for a set of x and y values.
    The FWHM is returned in the same units as those of x."""

    maxVal = np.max(y)
    maxVal50 = 0.5 * maxVal

    # this is to detect if there are multiple values
    biggerCondition = [a > maxVal50 for a in y]

    changePoints = []
    xPoints = []

    for k in range(len(biggerCondition) - 1):
        if biggerCondition[k + 1] != biggerCondition[k]:
            changePoints.append(k)

    assert len(
        changePoints) == 2, "More than two crossings of the threshold found."

    for k in changePoints:
        # do a polyfit
        # with the points before and after the point where the change occurs.

        # note that here we are fitting the x values as a function of the y values.
        # then we can use the polynom to compute the value of x at the threshold, i.e. at maxVal50.

        yPolyFit = x[k - 1:k + 2]
        xPolyFit = y[k - 1:k + 2]

        z = np.polyfit(xPolyFit, yPolyFit, 2)
        p = np.poly1d(z)
        xThis = p(maxVal50)
        xPoints.append(xThis)

    if len(xPoints) == 2:
        linewidth = xPoints[1] - xPoints[0]
    else:
        linewidth = None
    return linewidth


def result_worker(dir):
    try:
        result = NMRData(dir, "TopSpin", **kw)
        print('Exp. {} added.'.format(os.path.basename(dir)))
        print(result)
        return result
    except Exception as e:
        print("Problem adding {}. The error is: {}".format(
            os.path.basename(dir), e))
        # if kw['debug']:
        print(traceback.format_exc())
        # pass


def exp_process_worker(exp):
    exp.process()

def exp_integrate_worker(exp):
    exp.calculate_exp_integral()


def return_powers_csv(path, power_file):
    att_power = [[], []]
    if os.path.isfile(power_file):
        openfile = open(power_file)
    elif os.path.isfile(path + '/' + power_file):
        openfile = open(path + '/' + power_file, 'r')
    elif os.path.isfile(path + '/' + power_file + '.csv'):
        openfile = open(path + '/' + power_file + '.csv', 'r')
    else:
        raise Exception(
            "Error in reading power file. Please correct your input.")
    lines = openfile.readlines()
    if len(lines) == 1:
        lines = lines[0].split('\r')
    if len(lines[0].split('\r')[0].split(',')) == 2 and 'time' in lines[0]:
        print('This code is not compatible with "time, power" power logging.'
              '\nYou should have "time, dBm" or "dB, dBm" csv file.\n'
              'I continue evaluation but probably the output will not be usable.')
    elif len(lines[0].split('\r')[0].split(',')) == 2:
        lines.pop(0)
        for line in lines:
            att, power = line.split('\r')[0].split(',')
            att_power[0].append(float(att))
            att_power[1].append(float(power))
    elif len(lines[0].split('\r')[0].split(',')) == 3:
        lines.pop(0)
        for line in lines:
            time, power, att = line.split('\r')[0].split(',')
            att_power[0].append(float(att))
            att_power[1].append(float(power))
    else:
        raise Exception(
            "Could not use any power file in the directory. Aborting!")
    att_power = np.asarray(att_power)
    power_sets = set(att_power[0, :])
    powers_dict = {}
    for item in power_sets:
        if len(att_power[1, :][att_power[0, :] == item]) > 10:
            powers_dict[item] = np.sum(att_power[1, :][att_power[0, :] == item][10:-10]) / len(
                att_power[1, :][att_power[0, :] == item][10:-10])
        elif len(att_power[1, :][att_power[0, :] == item]) != 0:
            powers_dict[item] = np.sum(
                att_power[1, :][att_power[0, :] == item]) / len(att_power[1, :][att_power[0, :] == item])
    return powers_dict


def return_powers_mat(path, mat_files, **kwargs):
    """ Trying to go compatible with data from Songi Han's lab
        I do not know why byte order in Han's lab is big endians.
        This is to take care of that (not a good practice though). """
    kwargs['endianness'] = np.dtype('>i4')
    t1PowerMatFile = sio.loadmat(
        os.path.join(path, [i for i in main_files if i.endswith('.mat') and 't1' in i][0]))
    dnpPowerMatFile = sio.loadmat(
        os.path.join(path, [i for i in main_files if i.endswith('.mat') and 't1' not in i][0]))
    t1PowerMatFile = np.asarray(
        (t1PowerMatFile['timelist'], t1PowerMatFile['powerlist']))
    dnpPowerMatFile = np.asarray(
        (dnpPowerMatFile['timelist'], dnpPowerMatFile['powerlist']))
    return t1PowerMatFile, dnpPowerMatFile


def phase_dataset(results):
    """This is supposed to phase the dataset accroding to the kw['phase'] variable"""
    if kw['phase'] == 'first':
        exp = [result for result in results if result.expType == 'dnp'][0]
        for i in range(int(exp.phaseCycles)):
            jIndex = len(exp.allFid[4])/exp.phaseCycles
            applyIndex = [int(j*exp.phaseCycles)+i for j in range(int(jIndex))]
            exp.ph = exp.auto_phase0(
                4, 5, applyIndex[-1], -exp.ftWindow, exp.ftWindow, applyIndex=applyIndex)
            phase = exp.ph
    for i, exp in enumerate(results):
        if kw['phase'] == 'all' or exp.expType == 't1':
            for i in range(int(exp.phaseCycles)):
                jIndex = len(exp.allFid[4])/exp.phaseCycles
                applyIndex = [int(j*exp.phaseCycles) +
                              i for j in range(int(jIndex))]
                exp.ph = exp.auto_phase0(
                    4, 5, applyIndex[-1], -exp.ftWindow, exp.ftWindow, applyIndex=applyIndex)
                phase = 0
        elif kw['phase'] == 'first':
            exp.phase(4, 5, phase)
        if kw['phase'] == 'none':
            phase = 0
            exp.phase(4, 5, phase)
    return phase


def check_kwargs(path, **kwargs):
    kwargs['debug'] = kwargs.get('debug', False)
    kwargs['dumpToCsv'] = kwargs.get('dumpToCsv', False)
    kwargs['powerFile'] = kwargs.get('powerFile', '')
    kwargs['phase'] = kwargs.get('phase', 'first')
    kwargs['plotExts'] = kwargs.get('plotExts', [])
    kwargs['process'] = kwargs.get('process', True)
    kwargs['t1SeriesEval'] = kwargs.get('t1SeriesEval', True)
    kwargs['kSigmaCalc'] = kwargs.get('kSigmaCalc', True)
    kwargs['enhCalc'] = kwargs.get('enhCalc', True)
    kwargs['dumpToCsv'] = kwargs.get('dumpToCsv', True)
    kwargs['evalPath'] = kwargs.get('evalPath', 'eval')
    if os.listdir(path):
        kwargs['filesInDir'] = os.listdir(path)
    else:
        raise Exception(
            "The folder you picked for the experiment ({}) is empty.".format(path))
    return kwargs


def return_exps(path, **kwargs):
    global kw
    kw = check_kwargs(path, **kwargs)
    # inits
    dirs = []
    # We do not use these names anymore, this is for Han's lab
    ignore_list = [304, 503, 700, 701]
    # Taking care of evaluation dir
    try:
        os.mkdir(os.path.join(path, kw['evalPath']))
    except Exception as e:
        if '17' not in str(e):
            print('{} occured when trying to create evalPath'.format(e))
        pass
    # Plots ind. experiments and calculate T1, enhancement, etc.
    if not os.path.isdir(path):
        print('Error: the folder does not exist.')
        return False
    for name in kw['filesInDir']:
        try:
            name = int(name)
            if name not in ignore_list:
                dirs.append(name)
        except Exception as e:
            if kw['debug']:
                print('{} not NMR experiment({}).'.format(name, e))
            pass

    if not dirs:
        raise Exception('The folder does not contain any NMR experiments.')
    # sort and append path to dirs
    dirs = [os.path.join(path, str(dir)) for dir in sorted(dirs)]
    # Calculate powers the old way
    mat_files = glob.glob(os.path.join(path, '')+'*.mat')
    if kw['powerFile']:
        kw['powers_dict'] = return_powers_csv(path, kw['powerFile'])
    elif mat_files:
        kw['t1PowerMatFile'], kw['dnpPowerMatFile'] = return_powers_mat(
            path, mat_files, **kwargs)

    # Preparing experiment data
    start = time.time()
    # Single process
    # results = [result_worker(dir) for dir in dirs]
    # Multi-thread. Unfortunately, ProcessPoolExecutor cannot be used here since we pass unpickelable data
    with concurrent.futures.ThreadPoolExecutor() as pool:
        results = pool.map(result_worker, dirs)
    end = time.time()
    results = [result for result in results if not result == None]
    print('{} experiments were added in {:.2} seconds'.format(
        len(results), end-start))
    return results


def process_exps(results, path):
    # TODO
    # phasing
    with concurrent.futures.ThreadPoolExecutor() as pool:
        pool.map(exp_process_worker, results)
    phase = phase_dataset(results)
    with concurrent.futures.ThreadPoolExecutor() as pool:
        pool.map(exp_integrate_worker, results)

    # DNP enhancement
    if kw['enhCalc']:
        print(r"Fitting enhancement")
        calculate_dnp_enh(results, phase)
    else:
        kw['kSigmaCalc'] = False
        kw['enhCalc'] = False
        print('No DNP exps. found, no enhancement/kSigma curve will be created.')
    if kw['t1SeriesEval'] or kw['kSigmaCalc']:
        print(r"Fitting T1")
        if not calculate_t1_series(results):
            print(
                'Did not find any T1 experiment. no T1 series fitting, no kSigma calculation')
            kw['t1SeriesEval'] = False
            kw['kSigmaCalc'] = False
    if kw['kSigmaCalc'] and kw['t1SeriesEval']:
        print(r"Fitting kSigma")
        # expNum, powerMw, powerDbm, intReal, normIntReal, intMagn, normIntMagn, forward
        kw['kSigmaFit'] = k_sigma_calc()

    if kw['dumpToCsv']:
        print('Saving CSV files...')
        dumpAllToCSV(results, path=path)

    if kw['plotExts']:
        print('Plotting evaluation figures...')
        make_figures(results, path=path)

    print('All done')
    return results

def plot_exps_worker(result):
    try:
        if result.expType == 'dnp':  # FID plots
            print('Plotting exp {} figures (DNP)'.format(str(int(result.expNum))))
            figure = plt.figure(figsize=kw['figSize'])
            plt.plot(result.fidTimeHistory['bLeftShift'], np.real(result.allFid[0][0]), label='real')
            plt.plot(result.fidTimeHistory['bLeftShift'], np.imag(result.allFid[0][0]), label='imag')
            plt.title('FID at {:+.1f} dBm and {:.2f} mW power'.format(result.powerDbm, result.powerMw))
            plt.xlabel('time (ms)')
            plt.ylabel('Signal (a.u.)')
            plt.tight_layout()
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FID.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            print('plot saved at ', os.path.join(result.directory, ('FID')))
            plt.close(figure)
            kw['dnp_raw_fid_ax'].plot(result.fidTimeHistory['bLeftShift'], np.real(result.allFid[0][0]), label=(
                '{:+.1f} dBm\t{:.2f} mW power'.format(result.powerDbm, result.powerMw)).expandtabs())  # original
            kw['dnp_corrected_ax'].plot(result.fidTimeHistory['bZeroFilling'], np.real(result.allFid[1][0]), label=(
                '{:+.1f} dBm\t{:.2f} mW power'.format(result.powerDbm,
                                                        result.powerMw)).expandtabs())  # after left and right shift and baselineCorrection
            kw['dnp_zero_filled_ax'].plot(result.fidTime, np.real(result.allFid[2][0]), label=(
                '{:+.1f} dBm\t{:.2f} mW power'.format(result.powerDbm,
                                                        result.powerMw)).expandtabs())  # after zero filling
            kw['dnp_exp_win_ax'].plot(result.fidTime, np.real(result.allFid[3][0]), label=(
                '{:+.1f} dBm\t{:.2f} mW power'.format(result.powerDbm,
                                                        result.powerMw)).expandtabs())  # after exponential windowing
            # FT plots
            figure = plt.figure(figsize=kw['figSize'])
            plt.plot(result.frequency, np.real(
                result.allFid[5][0]), label='real')
            plt.plot(result.frequency, np.imag(
                result.allFid[5][0]), label='imag.')
            plt.plot(result.frequency, np.abs(
                result.allFid[5][0]), label='abs.')
            plt.title('{:+.1f} dBm, {:.2f} mW power and {}$^\circ$ phase'.format(
                result.powerDbm, result.powerMw, result.ph))
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.xlim(result.maxFreq - result.ftWindow,
                        result.maxFreq + result.ftWindow)
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FT.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # FT integral plots
            figure = plt.figure(figsize=kw['figSize'])
            plt.plot(result.frequency, np.real(
                result.allFid[5][0]), 'g-', label='real data', )
            plt.plot(result.frequency, np.imag(
                result.allFid[5][0]), 'y-', label='imag data', )
            plt.plot(result.frequency, np.abs(
                result.allFid[5][0]), 'r-', label='magn data', )
            plt.title('{:+.1f} dBm, {:.2f} mW power and {}$^\circ$ phase'.format(
                result.powerDbm, result.powerMw, result.ph))
            integration = result.integrate(
                5, 0, result.maxFreq - result.ftWindow, result.maxFreq + result.ftWindow, scale="Hz")
            plt.plot(
                integration['x-axis'], integration['integralReal'], 'g--', label='real integral',)
            plt.plot(
                integration['x-axis'], integration['integralImag'], 'y--', label='imag integral',)
            plt.plot(
                integration['x-axis'], integration['integralMagn'], 'r--', label='magn integral',)
            plt.title('FT+integral, {:+.1f} dBm, {:.2f} mW power and {}$^\circ$ phase'.format(
                result.powerDbm, result.powerMw, result.ph))
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.xlim(result.maxFreq - result.ftWindow,
                        result.maxFreq + result.ftWindow)
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FT_integral.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            kw['dnp_ft_real_ax'].plot(result.frequency, np.real(result.allFid[5][0]), label=(
                '{:.1f} dBm\t{:.2f} mW power'.format(result.powerDbm, result.powerMw)).expandtabs())
            kw['dnp_ft_real_ax'].set_title('FT after phasing to %.0f degrees' % result.ph)
            kw['dnp_ft_real_ax'].grid(True)
            kw['dnp_ft_real_ax'].set_xlim(result.maxFreq - result.ftWindow,
                            result.maxFreq + result.ftWindow)
            kw['dnp_ft_magn_ax'].plot(result.frequency, np.abs(result.allFid[5][0]) * result.real[1][0] / np.abs(result.real[1][0]), label=(
                '{:.1f} dBm\t{:.2f} mW power'.format(result.powerDbm, result.powerMw)).expandtabs())
            kw['dnp_ft_magn_ax'].set_title('FT after phasing to %.0f degrees' % result.ph)
            kw['dnp_ft_magn_ax'].grid(True)
            kw['dnp_ft_magn_ax'].set_xlim(result.maxFreq - result.ftWindow,
                            result.maxFreq + result.ftWindow)
            # Data for DNP figs
            # centerFreq.append([((result.expTime-expStart)/60.),
            #                    result.expCenterFreq, result.expNum, result.powerMw])
        elif result.expType == 't1' and kw['t1SeriesEval']:
            print('Plotting exp {} figures (T1)'.format(str(int(result.expNum))))
            figure = plt.figure(figsize=kw['figSize'])
            plt.errorbar(result.fitData[0], result.fitData[1],
                            yerr=(result.fitData[1] - result.t1fit['evalY']),
                            fmt='+',
                            label='data: '+str(kw['t1Calc']), capthick=2, capsize=2)
            plt.plot(result.t1fit['xdata'], result.t1fit['ydata'],
                        label=result.t1fit['t1FitFormula'])
            plt.annotate(r'$T_1={{{:.4f}}}\pm{{{:.4f}}}$'.format(
                result.t1fit['t1'], result.t1fit['t1error']), xy=(0.4, 0.5), xycoords='axes fraction')
            plt.title('$T_1$ at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('T1.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # Real part figure
            figure = plt.figure(figsize=kw['figSize'])
            for j in range(1, len(result.real)):
                plt.plot(result.real[0], result.real[j], label='PC' + str(j))
            plt.title('Real integral of phase cycle channels of $T_1$ at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.plot(result.phC[0], result.phC[1], '--',
                        label='PhaseCycled (sum after individual phasing)')
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('T1_PC_real.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # Magnitude part figure
            figure = plt.figure(figsize=kw['figSize'])
            for j in range(1, len(result.magn)):
                plt.plot(result.magn[0], result.magn[j], label='PC' + str(j))
            plt.title('Magnitude integral of phase cycle channels of $T_1$ at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.plot(result.phC[0], result.phC[2], '--',
                        label='PhaseCycled (sum after individual phasing)')
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('T1_PC_magn.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # Figure for fourier transforms
            figure = plt.figure(figsize=kw['figSize'])
            for j, time in enumerate(result.allFid[5]):
                if j % (result.phaseCycles) == 0:  # only first scan
                    plt.plot(result.frequency, np.real(result.allFid[5][j]),
                                label=('{:.2f}s'.format(result.vdList[int(j / result.phaseCycles)])))
            plt.title('FT (real) at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.xlim(result.maxFreq - result.ftWindow,
                        result.maxFreq + result.ftWindow)
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FT.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # Figure for FIDs
            # real
            figure = plt.figure(figsize=kw['figSize'])
            for j, time in enumerate(result.allFid[5]):
                if j % (result.phaseCycles) == 0:  # only first scan
                    plt.plot(result.fidTime, np.real(result.allFid[2][j]),
                                label='{:.2f}s'.format(result.vdList[int(j / result.phaseCycles)]))
            plt.title('$T_1$ FIDs at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.xlim(0, max(result.fidTime) / 3.)
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FID_corrected.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
            # After exponential windowing
            figure = plt.figure(figsize=kw['figSize'])
            for j, time in enumerate(result.allFid[5]):
                if j % (result.phaseCycles) == 0:  # only first scan
                    plt.plot(result.fidTime, np.real(result.allFid[3][j]),
                                label='{:.2f}s'.format(result.vdList[int(j / result.phaseCycles)]))
            plt.title('$T_1$ FIDs at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.xlim(0, max(result.fidTime))
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FID_after_windowing.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)

            figure = plt.figure(figsize=kw['figSize'])
            # after corrections
            for j, time in enumerate(result.allFid[5]):
                if j % (result.phaseCycles) == 0:  # only first scan
                    plt.plot(result.fidTimeHistory['bLeftShift'], np.real(result.allFid[0][j]),
                                label='{:.2f}s'.format(result.vdList[int(j / result.phaseCycles)]))
            plt.title('$T_1$ FIDs at {:.2f} mW power, {}$^\circ$ phase'.format(
                result.powerMw, result.ph))
            plt.xlabel('Time (s)')
            plt.ylabel('Intensity (a.u.)')
            plt.tight_layout()
            plt.legend(loc='best', fancybox=True,
                        shadow=True, fontsize='x-small')
            [plt.savefig(os.path.join(result.directory, ('FID_raw.' + x)), dpi=kw['plotDpi'])
                for x in kw['plotExts']]
            plt.close(figure)
        return result.expNum
    except Exception as e:
        return e


def dumpAllToCSV(results, path):
    # dnpEnh: expNum, powerMw, powerDbm, intReal, normIntReal, intMagn, normIntMagn, forward
    np.savetxt(os.path.join(path, kw['evalPath'], 'enhancements.csv'), kw['dnpEnh'], delimiter='\t', fmt='%s',
               header=('expNum\tpowerMw\tpowerDbm\tintReal\tnormIntReal\tintMagn\tnormIntMagn\tforward'))
    # t1Series: expNum, powerMw, powerDbm, t1, t1error
    if kw['t1SeriesEval']:
        np.savetxt(os.path.join(path, kw['evalPath'], kw['t1series']+'.csv'), kw['t1Series'], delimiter='\t', fmt='%s',
                   header=('expNum\tpowerMw\tpowerDbm\tt1\tt1error'))
        if kw['kSigmaCalc']:
            np.savetxt(os.path.join(path, kw['evalPath'], 'ksigma.csv'),
                       np.asarray(
                           (kw['dnpEnh'][:, 1], kw['kSigmaFit']['kSigmaCor'], kw['kSigmaFit']['kSigmaUncor'])).transpose(),
                       delimiter='\t', fmt='%s', header=('powerMw\tkSigmaCor\tkSigmaUncor'))

def make_figs_axes():
    """Initializes figures for DNP and T1 experiments"""
    # Initializing figures
    # Raw FIDs
    kw['dnp_raw_fid'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_raw_fid_ax'] = kw['dnp_raw_fid'].add_subplot(111)
    kw['dnp_raw_fid_ax'].set_title('Raw NMR FIDs (real)')
    kw['dnp_raw_fid_ax'].set_xlabel('time (ms)')
    kw['dnp_raw_fid_ax'].set_ylabel('Signal (a.u.)')
    # After digital filter and offset/baseline correction
    kw['dnp_corrected'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_corrected_ax'] = kw['dnp_corrected'].add_subplot(111)
    kw['dnp_corrected_ax'].set_title(
        'NMR FIDs after digital filter removal and offset/baseline correction (real)')
    kw['dnp_corrected_ax'].set_xlabel('time (ms)')
    kw['dnp_corrected_ax'].set_ylabel('Signal (a.u.)')
    # With zero filling
    kw['dnp_zero_filled'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_zero_filled_ax'] = kw['dnp_zero_filled'].add_subplot(111)
    kw['dnp_zero_filled_ax'].set_title('NMR FIDs with zero filling (real)')
    kw['dnp_zero_filled_ax'].set_xlabel('time (ms)')
    kw['dnp_zero_filled_ax'].set_ylabel('Signal (a.u.)')
    # After exponential windowing
    kw['dnp_exp_win'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_exp_win_ax'] = kw['dnp_exp_win'].add_subplot(111)
    kw['dnp_exp_win_ax'].set_title('NMR FIDs after exponential windowing (real)')
    kw['dnp_exp_win_ax'].set_xlabel('time (ms)')
    kw['dnp_exp_win_ax'].set_ylabel('Signal (a.u.)')
    # FT real
    kw['dnp_ft_real'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_ft_real_ax'] = kw['dnp_ft_real'].add_subplot(111)
    kw['dnp_ft_real_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_ft_real_ax'].set_xlabel('Frequency offset (Hz)')
    kw['dnp_ft_real_ax'].set_ylabel('Intensity (a.u.)')
    # FT magnitude
    kw['dnp_ft_magn'] = plt.figure(figsize=kw['figSize'])
    kw['dnp_ft_magn_ax'] = kw['dnp_ft_magn'].add_subplot(111)
    kw['dnp_ft_magn_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_ft_magn_ax'].set_xlabel('Frequency offset (Hz)')
    kw['dnp_ft_magn_ax'].set_ylabel('Intensity (a.u.)')
    # Time vs centerFreq
    kw['time_center_freq'] = plt.figure(figsize=kw['figSize'])  
    kw['time_center_freq_ax'] = kw['time_center_freq'].add_subplot(111)
    kw['time_center_freq_ax'].set_xlabel('time (mim)')
    kw['time_center_freq_ax'].set_ylabel('Center FT frequency (Hz)')
    # DNP enhancement
    kw['dnp_enhancement'] = plt.figure(figsize=kw['figSize'])  # DNP enhancement
    kw['dnp_enhancement_ax'] = kw['dnp_enhancement'].add_subplot(111)
    kw['dnp_enhancement_ax'].set_xlabel('Power (mW)')
    kw['dnp_enhancement_ax'].set_ylabel('FT Integral, normalized to MW off')
    # 3D figures
    if kw['make3dPlots']:
        kw['dnp_corrected_3d'] = plt.figure(figsize=kw['figSize'])
        kw['dnp_corrected_3d_ax'] = kw['dnp_corrected_3d'].gca(projection='3d')
        kw['dnp_corrected_3d_ax'].set_title(
            'NMR FIDs after baseline/offset correction (magn.)')
        kw['dnp_corrected_3d_ax'].set_xlabel('time (ms)')
        kw['dnp_corrected_3d_ax'].set_zlabel('Signal (a.u.)')
        kw['dnp_corrected_3d_ax'].set_ylabel('Experiment index')
        kw['dnp_corrected_3d_ax'].grid(True)
        kw['dnp_ft_real_3d'] = plt.figure(figsize=kw['figSize'])
        kw['dnp_ft_real_3d_ax'] = kw['dnp_ft_real_3d'].gca(projection='3d')
        kw['dnp_ft_real_3d_ax'].set_title('NMR FTs after exponential windowing (real)')
        kw['dnp_ft_real_3d_ax'].set_xlabel('Frequency offset (Hz)')
        kw['dnp_ft_real_3d_ax'].set_zlabel('Intensity (a.u.)')
        kw['dnp_ft_real_3d_ax'].set_ylabel('Experiment index')
        kw['dnp_ft_real_3d_ax'].grid(True)
        kw['dnp_ft_magn_3d'] = plt.figure(figsize=kw['figSize'])
        kw['dnp_ft_magn_3d_ax'] = kw['dnp_ft_magn_3d'].gca(projection='3d')
        kw['dnp_ft_magn_3d_ax'].set_title('NMR FTs after exponential windowing (magnitude)')
        kw['dnp_ft_magn_3d_ax'].set_xlabel('Frequency offset (Hz)')
        kw['dnp_ft_magn_3d_ax'].set_zlabel('Intensity (a.u.)')
        kw['dnp_ft_magn_3d_ax'].set_ylabel('Experiment index')
        kw['dnp_ft_magn_3d_ax'].grid(True)


def make_figures(results, path=''):
    """ Makes figures for ODNP experiment """
    make_figs_axes()
    if kw['enhCalc']:
        # Generating figures for dnp folders and collection of data
        powers = [x.powerMw for x in results if x.expType == 'dnp']
        inter = interp1d([min(powers), max(powers)], [0, 1])
        colors = [cm.jet(inter(x)) for x in powers]
        # DNP enhancement
        kw['dnp_enhancement_ax'].plot(
            kw['dnpEnh'][kw['dnpEnh'][:, 7] == 1][:, 1], kw['dnpEnh'][kw['dnpEnh'][:, 7] == 1][:, 6], 'bo', marker="o", label='forward magn.')
        kw['dnp_enhancement_ax'].plot(
            kw['dnpEnh'][kw['dnpEnh'][:, 7] == 0][:, 1], kw['dnpEnh'][kw['dnpEnh'][:, 7] == 0][:, 6], 'ro', marker="o", label='backward magn')
        kw['dnp_enhancement_ax'].plot(
            kw['dnpEnh'][kw['dnpEnh'][:, 7] == 1][:, 1], kw['dnpEnh'][kw['dnpEnh'][:, 7] == 1][:, 4], 'bo', marker="x", label='forward real')
        kw['dnp_enhancement_ax'].plot(
            kw['dnpEnh'][kw['dnpEnh'][:, 7] == 0][:, 1], kw['dnpEnh'][kw['dnpEnh'][:, 7] == 0][:, 4], 'ro', marker="x", label='backward real')
        for i in range(0, len(kw['dnpEnh'][:, 0])):
            if kw['dnpEnh'][i, 7] == 1:
                kw['dnp_enhancement_ax'].annotate('exp {:d}'.format(int(float(kw['dnpEnh'][i, 0]))),
                             xy=(kw['dnpEnh'][i, 1], kw['dnpEnh'][i, 6]),
                             xytext=(kw['dnpEnh'][i, 1] + (max(kw['dnpEnh'][:, 1]) -
                                                     min(kw['dnpEnh'][:, 1])) / 40, kw['dnpEnh'][i, 6]),
                             va='center', ha='left', size=9, color='blue', alpha=0.6)
            elif kw['dnpEnh'][i, 7] == 0:
                kw['dnp_enhancement_ax'].annotate('exp {:d}'.format(int(float(kw['dnpEnh'][i, 0]))),
                             xy=(kw['dnpEnh'][i, 1], kw['dnpEnh'][i, 6]),
                             xytext=(kw['dnpEnh'][i, 1] - (max(kw['dnpEnh'][:, 1]) -
                                                     min(kw['dnpEnh'][:, 1])) / 40, kw['dnpEnh'][i, 6]),
                             va='center', ha='right', size=9, color='red', alpha=0.6)
        try:
            kw['dnp_enhancement_ax'].plot(kw['enhancementFit']['xdata'], kw['enhancementFit']['ydata'], 'b--',
                     label='(magn.) ' + kw['enhancementFit']['enhancementFormula'])
            kw['dnp_enhancement_ax'].plot(kw['enhancementFit']['xdata'], kw['enhancementFit']['ydataExp'], 'g--',
                     label='(magn.) ' + kw['enhancementFit']['enhancementFormulaExp'])
            kw['dnp_enhancement_ax'].annotate(kw['enhancementFit']['annotation'], xy=(
                0.4, 0.5), xycoords='axes fraction', color='blue')
            kw['dnp_enhancement_ax'].annotate(kw['enhancementFit']['annotationExp'], xy=(
                0.55, 0.5), xycoords='axes fraction', color='green')
        except:
            pass
        kw['dnp_enhancement_ax'].set_title('Normalized DNP enhancement')
        kw['dnp_enhancement_ax'].legend(loc='best', fancybox=True, shadow=True, fontsize='x-small')
        kw['dnp_enhancement_ax'].legend(loc='upper right', fancybox=True,
                   shadow=True, fontsize='x-small')
        kw['dnp_enhancement'].tight_layout()
        [kw['dnp_enhancement'].savefig(os.path.join(path, kw['evalPath'], ('normalized_ODNP_enhancement.' + x)), dpi=kw['plotDpi'])
         for x in kw['plotExts']]
        plt.close(kw['dnp_enhancement'])
        print("Enhancement figure saved in {}".format(
            os.path.join(path, kw['evalPath'])))
    if kw['t1SeriesEval']:
        # Main T1 figure
        figure = plt.figure(figsize=kw['figSize'])
        # Generated linear fit
        plt.errorbar(kw['t1Series'][:, 1], kw['t1Series'][:, 3], yerr=kw['t1Series'][:, 4],
                     fmt='+', capthick=2, capsize=2, label=r'$T_1$ experiments')
        plt.plot(kw['t1FitSeries']['xdata'], kw['t1FitSeries']['ydata'], '--k',
                 label=kw['t1FitSeries']['fitFormula'])
        for i in range(0, len(kw['t1Series'][:, 0])):
            plt.annotate('exp {:d}'.format(int(kw['t1Series'][i, 0])),
                         xy=(kw['t1Series'][i, 1], kw['t1Series'][i, 3]),
                         xytext=(kw['t1Series'][i, 1] + (max(kw['t1Series'][:, 1]) -
                                                   min(kw['t1Series'][:, 1])) / 40, kw['t1Series'][i, 3]),
                         va='center', ha='left', size=9, color='blue', alpha=0.6)
        plt.xlabel('Power (mW)')
        plt.ylabel('Time (s)')
        plt.title(r'$T_1$ series')
        plt.legend()
        figure.tight_layout()
        [plt.savefig(os.path.join(path, kw['evalPath'], ('T1_time_series.' + x)), dpi=kw['plotDpi'])
         for x in kw['plotExts']]
        plt.close(figure)
        print("T1 figure saved in {}".format(os.path.join(path, kw['evalPath'])))
        if kw['kSigmaCalc']:
            # kSigma figure
            figure = plt.figure(figsize=kw['figSize'])
            plt.plot(kw['dnpEnh'][:, 1], kw['kSigmaFit']['kSigmaCor'],
                     'o', c='green', label='cor')
            plt.plot(kw['dnpEnh'][:, 1], kw['kSigmaFit']['kSigmaUncor'],
                     'o', c='red', label='uncor')
            plt.plot(kw['kSigmaFit']['xdata'], kw['kSigmaFit']['ydataCor'],
                     '--', c='green', label=kw['kSigmaFit']['corFormula'])
            plt.plot(kw['kSigmaFit']['xdata'], kw['kSigmaFit']['ydataUncor'],
                     '--', c='red', label=kw['kSigmaFit']['uncorFormula'])
            plt.xlabel('Power (mW)')
            plt.ylabel(r'$k_{\sigma}s(P)C$  $(\frac{1}{s})$')
            plt.title(r'$k_{\sigma}$ calculations')
            plt.legend()
            figure.tight_layout()
            [plt.savefig(os.path.join(path, kw['evalPath'], ('kSigma.' + x)), dpi=kw['plotDpi'])
             for x in kw['plotExts']]
            plt.close(figure)
            print("kSigma figure saved in {}".format(
                os.path.join(path, kw['evalPath'])))

    kw['figs'] = {}
    # Single thread
    test = [plot_exps_worker(result) for result in results]
    # with concurrent.futures.ThreadPoolExecutor() as pool:
    #     test = list(pool.map(plot_exps_worker, results))
    # concurrent.futures.wait(test, timeout=None)
    # concurrent.futures.wait(test)
    print('The result of concurrent on figures is: ')
    print(list(test))
    
    # DNP figures
    print("saving figures in {}".format(os.path.join(path, kw['evalPath'])))
    # saving figures in evaluation directory
    kw['dnp_raw_fid_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_raw_fid'].tight_layout()
    [kw['dnp_raw_fid'].savefig(os.path.join(path, kw['evalPath'], ('01_FIDs_raw.' + x)),
                  dpi=kw['plotDpi']) for x in kw['plotExts']]
    plt.close(kw['dnp_raw_fid'])
    del kw['dnp_raw_fid']
    kw['dnp_corrected_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_corrected'].tight_layout()
    [kw['dnp_corrected'].savefig(os.path.join(path, kw['evalPath'], ('02_FIDs_after_LS_RS_baseline.' + x)),
                  dpi=kw['plotDpi']) for x in kw['plotExts']]
    plt.close(kw['dnp_corrected'])
    del kw['dnp_corrected']
    kw['dnp_zero_filled_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_zero_filled'].tight_layout()
    [kw['dnp_zero_filled'].savefig(os.path.join(path, kw['evalPath'], ('03_FIDs_after_zero_filling.' + x)), dpi=kw['plotDpi'])
     for x in kw['plotExts']]
    plt.close(kw['dnp_zero_filled'])
    del kw['dnp_zero_filled']
    kw['dnp_exp_win_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_exp_win'].tight_layout()
    [kw['dnp_exp_win'].savefig(os.path.join(path, kw['evalPath'], ('04_FIDs_after_exp_windowing.' + x)), dpi=kw['plotDpi'])
     for x in kw['plotExts']]
    plt.close(kw['dnp_exp_win'])
    del kw['dnp_exp_win']
    if kw['make3dPlots']:
        try:
            [kw['dnp_corrected_3d_ax'].plot3D(value.fidTimeHistory['bZeroFilling'], np.abs(value.allFid[1][0]), value.expNum, zdir='y',
                          zorder=int(-value.expNum), color=colors[i]) for i, value in
             enumerate(x for x in results if x.expType == 'dnp')]
            xmin, xmax = min([min(a.fidTimeHistory['bZeroFilling']) for a in results if a.expType == 'dnp']), max(
                [max(a.fidTimeHistory['bZeroFilling']) for a in results if a.expType == 'dnp'])
            zmin, zmax = min([min(np.abs(a.allFid[1][0])) for a in results if a.expType == 'dnp']), max(
                [max(np.abs(a.allFid[1][0])) for a in results if a.expType == 'dnp'])
            zs = [a.expNum for a in results if a.expType == 'dnp']
            labels = ['{:d} ({:.2f} W)'.format(int(a.expNum), a.powerMw / 1000)
                      for a in results if a.expType == 'dnp']
            plt.yticks(zs, labels,
                       rotation=270)
            for i, label in enumerate(kw['dnp_corrected_3d_ax'].get_yticklabels()):
                if label.get_text():
                    label.set_color(colors[i])
            kw['dnp_corrected_3d_ax'].set_xlim3d(xmin, xmax)
            kw['dnp_corrected_3d_ax'].set_zlim3d(zmin, zmax)
            kw['dnp_corrected_3d_ax'].set_ylim3d(min(zs), max(zs))
            kw['dnp_corrected_3d_ax'].yaxis.labelpad = 45
            kw['dnp_corrected_3d_ax'].view_init(elev=17., azim=-23.)
            kw['dnp_corrected_3d'].tight_layout()
            [kw['dnp_corrected_3d'].savefig(os.path.join(path, kw['evalPath'], ('01_FIDs_raw_3D.' + x)), dpi=kw['plotDpi'])
             for x in kw['plotExts']]
            plt.close(kw['dnp_corrected_3d'])
            xmin, xmax = min([a.maxFreq - a.ftWindow for a in results if a.expType == 'dnp']), max(
                [a.maxFreq + a.ftWindow for a in results if a.expType == 'dnp'])
            zmin, zmax = min([min(np.real(a.allFid[5][0])) for a in results if a.expType == 'dnp']), max(
                [max(np.real(a.allFid[5][0])) for a in results if a.expType == 'dnp'])
            verts = [list(zip(a.frequency[np.logical_and(a.frequency > xmin, a.frequency < xmax)],
                              np.real(a.allFid[5][0][np.logical_and(a.frequency > xmin, a.frequency < xmax)]))) for a in
                     results if a.expType == 'dnp']
            poly = PolyCollection(verts, facecolors=colors)
            poly.set_alpha(0.6)
            plt.yticks(zs, labels, rotation=270)
            for i, label in enumerate(kw['dnp_ft_real_3d_ax'].get_yticklabels()):
                if label.get_text():
                    label.set_color(colors[i])
            kw['dnp_ft_real_3d_ax'].view_init(elev=17., azim=-23.)
            kw['dnp_ft_real_3d_ax'].set_xlim3d(xmin, xmax)
            kw['dnp_ft_real_3d_ax'].set_zlim3d(zmin, zmax)
            kw['dnp_ft_real_3d_ax'].set_ylim3d(min(zs), max(zs))
            kw['dnp_ft_real_3d_ax'].add_collection3d(poly, zs=zs, zdir='y')
            kw['dnp_ft_real_3d_ax'].view_init(elev=17., azim=-23.)
            kw['dnp_ft_real_3d_ax'].yaxis.labelpad = 45
            kw['dnp_ft_real_3d_ax'].legend(loc='upper right', fancybox=True,
                         shadow=True, fontsize='x-small')
            kw['dnp_ft_real_3d'].tight_layout()
            [kw['dnp_ft_real_3d'].savefig(os.path.join(path, kw['evalPath'], ('05_FT_after_phasing_real_3d.' + x)), dpi=kw['plotDpi'])
             for x in kw['plotExts']]
            plt.close(kw['dnp_ft_real_3d'])
            del kw['dnp_ft_real_3d']
            verts = [list(zip(a.frequency[np.logical_and(a.frequency > xmin, a.frequency < xmax)], (
                np.abs(a.allFid[5][0][np.logical_and(a.frequency > xmin, a.frequency < xmax)]) * a.real[1][0] / np.abs(
                    a.real[1][0])))) for a in results if a.expType == 'dnp']
            zmin, zmax = min(
                [min(np.abs(a.allFid[5][0]) * a.real[1][0] / np.abs(a.real[1][0])) for a in results if
                 a.expType == 'dnp']), max(
                [max(np.abs(a.allFid[5][0]) * a.real[1][0] / np.abs(a.real[1][0])) for a in results if a.expType == 'dnp'])
            zs = [a.expNum for a in results if a.expType == 'dnp']
            poly = PolyCollection(verts, facecolors=colors)
            poly.set_alpha(0.6)
            plt.yticks(zs, labels, rotation=270)
            for i, label in enumerate(kw['dnp_ft_magn_3d_ax'].get_yticklabels()):
                if label.get_text():
                    label.set_color(colors[i])
            kw['dnp_ft_magn_3d_ax'].set_xlim3d(xmin, xmax)
            kw['dnp_ft_magn_3d_ax'].set_zlim3d(zmin, zmax)
            kw['dnp_ft_magn_3d_ax'].set_ylim3d(min(zs), max(zs))
            kw['dnp_ft_magn_3d_ax'].add_collection3d(poly, zs=zs, zdir='y')
            kw['dnp_ft_magn_3d_ax'].view_init(elev=17., azim=-23.)
            kw['dnp_ft_magn_3d_ax'].yaxis.labelpad = 45
            kw['dnp_ft_magn_3d_ax'].legend(loc='upper right', fancybox=True,
                         shadow=True, fontsize='x-small')
            # dnp_ft_magn_3d.colorbar(test)
            kw['dnp_ft_magn_3d'].tight_layout()
            [kw['dnp_ft_magn_3d'].savefig(os.path.join(path, kw['evalPath'], ('06_FT_after_phasing_magn_3d.' + x)), dpi=kw['plotDpi'])
             for x in kw['plotExts']]
            plt.close(kw['dnp_ft_magn_3d'])
            del kw['dnp_ft_magn_3d']
        except Exception as e:
            print('Error "{}" occured while making 3D plots'.format(e))
    kw['dnp_ft_real_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_ft_real'].tight_layout()
    [kw['dnp_ft_real'].savefig(os.path.join(path, kw['evalPath'], ('05_FT_after_phasing_real.' + x)), dpi=kw['plotDpi'])
     for x in kw['plotExts']]
    plt.close(kw['dnp_ft_real'])
    kw['dnp_ft_magn_ax'].legend(loc='upper right', fancybox=True,
               shadow=True, fontsize='x-small')
    kw['dnp_ft_magn'].tight_layout()
    [kw['dnp_ft_magn'].savefig(os.path.join(path, kw['evalPath'], ('06_FT_after_phasing_magn.' + x)), dpi=kw['plotDpi'])
     for x in kw['plotExts']]
    plt.close(kw['dnp_ft_magn'])
    del kw['dnp_ft_magn']
    plt.close(kw['time_center_freq'])
    del kw['time_center_freq']


def fit_t1(time, intensity, si00=-1e7, t10=.5, c0=-1e7, spaceNo=500, t1ErrorTol=0.5):  # T1 fitting
    def t1ir(x, si0, c, t1):  # Defines the T1 function.
        return si0 + (-c - si0) * np.exp(-x / t1)

    # calls the curve fitting routine using the function described above
    popt, pcov = curve_fit(t1ir, time, intensity, [si00, c0, t10])
    # def fun(x, t, y):# Defines minimization of T1 for least squares (It's not needed since curve_fit also uses least sq)
    #    return x[0]+(-x[1]-x[0])*np.exp(-t/x[2]) - y
    # res = least_squares(fun, x0 = popt, args = (time, intensity)) # Do least squares minimization
    xdata = np.linspace(min(time), max(time), spaceNo)
    ydata = t1ir(xdata, *popt)
    rmsd = np.sqrt(
        ((np.array(intensity) - np.array(t1ir(time, *popt))) ** 2).mean())
    if np.sqrt(pcov[2, 2]) > t1ErrorTol:
        raise Exception('T1 error should not exceed t1ErrorTol value. T1 error is {:0.3} and t1ErrorTol {:0.3}'.format(
            np.sqrt(pcov[2, 2]), float(t1ErrorTol)))
    return {
        'xdata': xdata,
        'ydata': ydata,
        'evalY': t1ir(time, *popt),
        't1': popt[2],
        't1FitFormula': r"$M(t)={0:.2e}+({1:+.2e}{0:+.2e})exp(\frac{{-t}}{{{2:.3f}}})$".format(*popt),
        't1error': np.sqrt(pcov[2, 2]),
        'rmsd': rmsd,
    }


def fit_enhancement(power, enhancement):
    def func(x, a, b, c, d):
        return a * np.exp(-b * x ** c) + d

    def func_exp(x, a, b, d):
        return a * np.exp(-b * x) + d

    poptExp, pcovExp = curve_fit(func_exp, power, enhancement, maxfev=50000)
    popt, pcov = curve_fit(func, power, enhancement, maxfev=50000, p0=(
        poptExp[0], poptExp[1], 1, poptExp[2]))
    print("Enhancement fit values are {} for normal and {} for exponential fit".format(
        popt, poptExp))
    eMax = func(np.inf, *popt)
    eP = 1 - (1 - eMax) / 2
    sd = 0  # square deviation
    eMaxExp = func_exp(np.inf, *poptExp)
    ePExp = 1 - (1 - eMaxExp) / 2
    sdExp = 0  # square deviation
    for i in range(len(power)):
        sd += np.square(enhancement[i] - func(power[i], *popt))
        sdExp += np.square(enhancement[i] - func_exp(power[i], *poptExp))
    rmsd = np.sqrt(sd / len(power))
    rmsdExp = np.sqrt(sdExp / len(power))
    xdata = np.linspace(min(power), max(power), 50000)
    ydata = func(xdata, *popt)
    ydataExp = func_exp(xdata, *poptExp)
    eHalf = min(abs(ydata - eP))
    eHalfExp = min(abs(ydataExp - ePExp))
    pHalf = np.asarray(list(zip(xdata, ydata)))[abs(ydata - eP) == eHalf][0][0]
    pHalfExp = np.asarray(list(zip(xdata, ydataExp)))[
        abs(ydataExp - ePExp) == eHalfExp][0][0]
    return {
        'rmsd': rmsd,
        'xdata': xdata,
        'ydata': ydata,
        'eHalf': eHalf,
        'pHalf': pHalf,
        'eMax': eMax,
        'eP': eP,
        'enhancementFormula': r"$E(P) = {:5.2f}exp({:+5.3f} P^{{{:5.2f}}}) {:+5.2f}$, RMSD = {:5.2f}".
        format(*popt, rmsd),
        'annotation': r"$S_{{rel}}(P) = \frac{{P / {:5.2f} }}{{1+P / {:5.2f} }}$"
                      "\n"
                      r"$E_{{max}} = {:5.2f}$".format(
                          pHalf, pHalf, func(np.inf, *popt)),
        'rmsdExp': rmsdExp,
        'ydataExp': ydataExp,
        'eHalfExp': eHalfExp,
        'pHalfExp': pHalfExp,
        'eMaxExp': eMaxExp,
        'ePExp': ePExp,
        'enhancementFormulaExp': r"$E(P) = {:5.2f}exp({:+5.3f} P) {:+5.2f}$, RMSD = {:5.2f}".
        format(*poptExp, rmsd),
        'annotationExp': r"$S_{{rel}}(P) = \frac{{P / {:5.2f} }}{{1+P / {:5.2f} }}$"
                         "\n"
                         r"$E_{{max}} = {:5.2f}$".format(
                             pHalfExp, pHalfExp, func_exp(np.inf, *poptExp)),
    }


def k_sigma_calc(spaceNo=500):
    kSigmaUncor = [((1 - kw['dnpEnh'][:, 6][i]) * .00152 / kw['t1FitSeries']['fit'](kw['dnpEnh'][:, 1][0]))
                   for i in range(len(kw['dnpEnh'][:, 1]))]
    kSigmaCor = [((1 - kw['dnpEnh'][:, 6][i]) * .00152 / kw['t1FitSeries']['fit'](kw['dnpEnh'][:, 1][i]))
                 for i in range(len(kw['dnpEnh'][:, 1]))]

    def func(x, a, b):
        return a * x / (b + x)

    poptCor, pcovCor = curve_fit(
        func, kw['dnpEnh'][:, 1], kSigmaCor, maxfev=10000)
    poptUncor, pcovUncor = curve_fit(
        func, kw['dnpEnh'][:, 1], kSigmaUncor, maxfev=10000)
    kSigmaSmaxCor = func(10e100, *poptCor)
    kSigmaSmaxUncor = func(10e100, *poptUncor)
    xdata = np.linspace(min(kw['dnpEnh'][:, 1]),
                        max(kw['dnpEnh'][:, 1]), spaceNo)
    return {
        'kSigmaCor': kSigmaCor,
        'kSigmaUncor': kSigmaUncor,
        'fitCor': lambda x: func(x, *poptCor),
        'fitUncor': lambda x: func(x, *poptUncor),
        'xdata': xdata,
        'ydataCor': func(xdata, *poptCor),
        'ydataUncor': func(xdata, *poptUncor),
        'kSigmaSmaxCor': kSigmaSmaxCor,
        'kSigmaSmaxUncor': kSigmaSmaxUncor,
        'corFormula': r'$k_{{\sigma}}s(P)C = \frac{{{:.2e}\times P}}{{{:.2f}+P}}$'
                      '\n'
                      r'$k_{{\sigma}}s_{{max}}C = {:,.2e}$ $(\frac{{1}}{{s}})$'.format(
                          *poptCor, kSigmaSmaxCor),
        'uncorFormula': r'$k_{{\sigma}}s(P)C = \frac{{{:.2e}\times P}}{{{:.2f}+P}}$'
                        '\n'
                        r'$k_{{\sigma}}s_{{max}}C = {:,.2e}$ $(\frac{{1}}{{s}})$'.format(
                            *poptUncor, kSigmaSmaxUncor),
    }


def fit_t1_series(power, t1, t1error=0, degree=1, spaceNo=500):
    coefs = np.polyfit(np.asarray(power), np.asarray(t1),
                       degree, w=1 / np.asarray(t1error))
    fit = np.poly1d(coefs)
    xdata = np.linspace(min(power), max(power), spaceNo)
    ydata = fit(xdata)
    fitFormula = ''
    for i, coeff in enumerate(coefs):
        if i < len(coefs) - 2:
            fitFormula += r'${:+.1e}P^{:d}$'.format(coeff, len(coefs) - i - 1)
        elif i < len(coefs) - 1:
            fitFormula += r'${:+.1e}P$'.format(coeff)
        else:
            fitFormula += r'${:+.2f}$'.format(coeff)
    # line = slope*power+intercept
    # rmsd = np.sqrt(((np.array(t1)-np.array(line))**2).mean())
    return {
        'fit': fit,
        'coefs': coefs,
        'xdata': xdata,
        'ydata': ydata,
        'fitFormula': fitFormula,
    }


def cc(arg):
    return mcolors.to_rgba(arg, alpha=0.6)


def find_nearest(array, values):
    indices = np.abs(np.subtract.outer(array, values)).argmin(0)
    return indices, array[indices]


def calculate_dnp_enh(results, phase):
    if len([res for res in results if res.expType == 'dnp']) < 2:
        return False
    powerMw = -1
    dnpEnh = []  # expNum, powerMw, powerDbm, intReal, normIntReal, intMagn, normIntMagn
    dnpCounter = 0
    enhMinPower = [results[i] for i in range(len(results)) if results[i].magn[1][0] == min(
        [result.magn[1][0] for result in results if result.expType == 'dnp'])][0].powerMw
    for i, result in enumerate(results):
        if result.expType == 'dnp':
            if dnpCounter == 0:
                normReal = result.real[1][0]
                normMagn = result.magn[1][0]
            # Work-around for phasing problem
            if phase == 'all' and result.powerMw > enhMinPower:
                dnpEnhLine = [int(result.expNum), result.powerMw, result.powerDbm, -result.real[1][0],
                              -result.real[1][0] / normReal, -result.magn[1][0], -result.magn[1][0] / normMagn]
            else:
                dnpEnhLine = [int(result.expNum), result.powerMw, result.powerDbm, result.real[1][0],
                              result.real[1][0] / normReal, result.magn[1][0], result.magn[1][0] / normMagn]

            if result.powerMw >= powerMw:
                # expNum, powerMw, powerDbm, intReal, normIntReal, intMagn, normIntMagn, forward
                dnpEnhLine.append(1)
            else:
                dnpEnhLine.append(0)
            dnpEnh.append(dnpEnhLine)
            dnpCounter += 1
            powerMw = result.powerMw
    dnpEnh = np.asanyarray(dnpEnh)
    try:
        enhancementFit = fit_enhancement(dnpEnh[:, 1], dnpEnh[:, 6])
    except Exception as e:
        print("Error {} happened when fitting enhancements".format(e))
        enhancementFit = False
    print(r"Fitting enhancement")
    kw['dnpEnh'] = dnpEnh
    kw['enhancementFit'] = enhancementFit
    return True


def calculate_t1_series(results):
    if len([res for res in results if res.expType == 't1']) < 1:
        return False
    t1Series = []
    t1Counter = 0
    for result in results:
        if result.expType == 't1':
            # expNum, powerMw, powerDbm, t1, t1error
            t1Series.append([result.expNum,
                             result.powerMw,
                             result.powerDbm,
                             result.t1fit['t1'],
                             result.t1fit['t1error']])
            t1Counter += 1
    t1Series = np.asarray(t1Series)
    print(t1Series)
    t1SeriesPolDeg = kw.get('t1SeriesPolDeg', 1)
    try:
        t1FitSeries = fit_t1_series(
            t1Series[:, 1], t1Series[:, 3], t1Series[:, 4], degree=t1SeriesPolDeg)
    except Exception as e:
        print("Error {} happened when fitting T1 series".format(e))
        t1FitSeries = False
    kw['t1Series'] = t1Series
    kw['t1FitSeries'] = t1FitSeries
    return True
