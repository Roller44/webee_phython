import numpy as np
import copy
import struct

import pdb 
import matplotlib.pylab as plt
# Debugging.

WITHCP = 80
#WITHCP = 64
WITHOUTCP = 64
SAMPLERATE = int(WITHCP / 4)

ALIGNMENT = [1,0,1,0] #1 means RIGHT and 0 means LEFT

def complexmultiply(a,b):
	# This function mimics the muliplication between two complex numbers.
	# Note that although the output should be a complex number, its imagine
	# part does not include "j".
	c = [0, 0]
	c[0] = a[0]*b[0] - a[1]*b[1]
	c[1] = a[1]*b[0] + a[0]*b[1]
	return c

def cycleshift(symbol, offset):
	length = len(symbol)
	tsymbol = [0 for i in range(length)]
	for i in range(length):
		tsymbol[i] = symbol[i]

	csymbol = 0
	if offset >= 0:
		# e.g. when offset = 2, we have [1, 2, 3, 4, 5] -> [3, 4, 5, 1, 2]
		csymbol = tsymbol[offset:]+tsymbol[:offset]
		
	if offset <= 0:
		# e.g. when offset = -2 , we have [1, 2, 3, 4, 5] -> [4, 5, 1, 2, 3]
		offset = length + offset
		csymbol = tsymbol[offset:]+tsymbol[:offset]

	return csymbol

def generateQAM(N, side):
	# Generate all N possible QAM points (at WiFi).
	# "side" is the maximum absolute values of inphase and quadrature.
	qam = [[0,0] for i in range(N)]
	m = int(np.sqrt(N))
	interval = float(side * 2) / (m - 1)
	xstart = -1 * side
	ystart = -1 * side
	for i in range(m):
		for j in range(m):
			x = xstart + i * interval
			y = ystart + j * interval
			qam[i*m + j][0] = float("%.5f" %x) # type: ignore
			qam[i*m + j][1] = float("%.5f" %y) # type: ignore
	return qam

def Edis(point1,point2):
	N = len(point1)
	temp = 0
	for i in range(N):
		temp = temp + (point1[i] - point2[i])**2
	dis = np.sqrt(temp)
	return dis

def approxQAM(qam, IQs):
	# This fnction approximate the WiFi sampling instanceso ver a ZigBee symbol
	# to standard WiFi QAM points.
	# qam: values of all useable standard WiFi QAM points.
	# IQs: ideal WiFi sampling instances over a ZigBee symbol.
	N = len(IQs)
	M = len(qam)
	approxIQs = []
	for i in range(N):
		x1 = IQs[i][0]
		y1 = IQs[i][1]
		if x1 == 0 and y1 == 0:
			approxIQs.append([0,0])
		else:
			# Find the wave_index of qam whose value is the most similar to a given ideal
			# WiFi sampling instance.
			mindis = 1000
			minindex = -1
			for j in range(M):
				x2 = qam[j][0]
				y2 = qam[j][1]
				dis = Edis([x1,y1],[x2,y2])
				if dis < mindis:
					mindis = dis
					minindex = j

			approxIQs.append(qam[minindex])

	return approxIQs
def get_zigbee_signal(symbol, ratio):
	# This function generates the waveform of a ZigBee symbol.
	# The output IQ values are WiFi sampling instances over the ZigBee symbol.

	# Chip map constuction.
	chip_map = np.zeros((16, 32))
	chip_map[0] = [1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,0]
	chip_map[1] = [1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0]
	chip_map[2] = [0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0]
	chip_map[3] = [0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1]
	chip_map[4] = [0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1]
	chip_map[5] = [0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0]
	chip_map[6] = [1,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1]
	chip_map[7] = [1,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1]
	chip_map[8] = [1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1]
	chip_map[9] = [1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1]
	chip_map[10] = [0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1]
	chip_map[11] = [0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0]
	chip_map[12] = [0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0]
	chip_map[13] = [0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1]
	chip_map[14] = [1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0]
	chip_map[15] = [1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0]
	# Symbol-to-chip_seq.
	chip_seq = chip_map[symbol]
	# Set time values of WiFi sampling instances.
	# Each half-sine wave sampled 20 times by WiFi.
	num_sample_per_wave = SAMPLERATE
	# Generate waves.
	wave = np.sin(np.arange(0, np.pi, np.pi/num_sample_per_wave)) 
	wave = np.tile(wave, 16)
	# Generate in-phase waves.
	I_chip_seq = chip_seq[::2]
	I_chip_seq = (I_chip_seq - 0.5) * 2
	I_chip_seq = np.repeat(I_chip_seq, num_sample_per_wave)
	I_wave = I_chip_seq * wave
	I_wave = np.concatenate((I_wave, np.zeros(int(num_sample_per_wave / 2)))) * ratio
	# Generate quadrant waves.
	Q_chip_seq = chip_seq[1::2]
	Q_chip_seq = (Q_chip_seq - 0.5) * 2
	Q_chip_seq = np.repeat(Q_chip_seq, num_sample_per_wave)
	Q_wave = Q_chip_seq * wave
	Q_wave = np.concatenate((np.zeros(int(num_sample_per_wave / 2)), Q_wave)) * ratio

	return np.vectorize(complex)(I_wave, Q_wave)

def generateIQ(symbol, ratio):
	# This function generates the waveform of a ZigBee symbol.
	# The output IQ values are WiFi sampling instances over the ZigBee symbol.

	# Chip map constuction.
	chips = np.zeros((16, 32))
	chips[0] = [1,1,0,1,  1,0,0,1,  1,1,0,0, 0,0,1,1,  0,1,0,1, 0,0,1,0, 0,0,1,0, 1,1,1,0]
	chips[1] = [1,1,1,0,  1,1,0,1,  1,0,0,1,  1,1,0,0, 0,0,1,1,  0,1,0,1, 0,0,1,0, 0,0,1,0]
	chips[2] = [0,0,1,0, 1,1,1,0,  1,1,0,1,  1,0,0,1,  1,1,0,0, 0,0,1,1,  0,1,0,1, 0,0,1,0]
	chips[3] = [0,0,1,0, 0,0,1,0, 1,1,1,0,  1,1,0,1,  1,0,0,1,  1,1,0,0, 0,0,1,1,  0,1,0,1]
	chips[4] = [0,1,0,1, 0,0,1,0, 0,0,1,0, 1,1,1,0,  1,1,0,1,  1,0,0,1,  1,1,0,0, 0,0,1,1]
	chips[5] = [0,0,1,1,  0,1,0,1, 0,0,1,0, 0,0,1,0, 1,1,1,0,  1,1,0,1,  1,0,0,1,  1,1,0,0] 
	chips[6] = [1,1,0,0, 0,0,1,1,  0,1,0,1, 0,0,1,0, 0,0,1,0, 1,1,1,0,  1,1,0,1,  1,0,0,1]
	chips[7] = [1,0,0,1,  1,1,0,0, 0,0,1,1,  0,1,0,1,  0,0,1,0, 0,0,1,0, 1,1,1,0, 1,1,0,1]
	chips[8] = [1,0,0,0,  1,1,0,0,  1,0,0,1, 0,1,1,0,  0,0,0,0, 0,1,1,1, 0,0,1,1, 1,0,1,1]
	chips[9] = [1,0,1,1, 1,0,0,0,  1,1,0,0,  1,0,0,1, 0,1,1,0,  0,0,0,0, 0,1,1,1, 0,0,1,1]
	chips[10] = [0,0,1,1, 1,0,1,1, 1,0,0,0,  1,1,0,0,  1,0,0,1, 0,1,1,0,  0,0,0,0, 0,1,1,1]
	chips[11] = [0,1,1,1, 0,0,1,1, 1,0,1,1, 1,0,0,0,  1,1,0,0,  1,0,0,1, 0,1,1,0,  0,0,0,0]
	chips[12] = [0,0,0,0, 0,1,1,1, 0,0,1,1, 1,0,1,1, 1,0,0,0,  1,1,0,0,  1,0,0,1, 0,1,1,0]
	chips[13] = [0,1,1,0,  0,0,0,0, 0,1,1,1, 0,0,1,1, 1,0,1,1, 1,0,0,0,  1,1,0,0,  1,0,0,1]
	chips[14] = [1,0,0,1, 0,1,1,0,  0,0,0,0, 0,1,1,1, 0,0,1,1, 1,0,1,1, 1,0,0,0,  1,1,0,0]
	chips[15] = [1,1,0,0, 1,0,0,1, 0,1,1,0,  0,0,0,0, 0,1,1,1, 0,0,1,1, 1,0,1,1, 1,0,0,0]

	# Symbol-to-chip sequence.
	chip = chips[symbol]
	chip_len = len(chip)
	IQs = []

	# Divide a chip sequence into an I-chip sequence and a Q-chip seuqnce.
	Ichips = []
	Qchips = []
	for i in range(chip_len):
		if i % 2 == 0:
			Ichips.append(chip[i])
		else:
			Qchips.append(chip[i])

	# Set start time of each half-sine wave.
	# Here, each chip_seq is encoded into the amplitude of a half-sine wave.
	Interval = []
	for i in range(chip_len + 1):
		Interval.append(i*np.pi / 2)

	# Set time values of WiFi sampling instances.
	# Each half-sine wave sampled 20 times by WiFi.
	SampleRate = SAMPLERATE
	X = []
	for i in range(16):
		for j in range(int(SampleRate)):
			X.append(i*np.pi + j*np.pi / SampleRate)

	# Set value of each WiFi sampling instance.
	for i in range(len(X)):
		# Figure out which chip_seq the i-th sampling instance belongs to.
		x = X[i]
		index = 0
		while x > Interval[index]:
			index = index + 1 
		index = index - 1
		Iindex = index / 2
		Qindex = (index - 1) / 2

		if Qindex < 0:
			Qindex = 0

		# Let amplitudes 1 and -1 represent chip_seq values 1 and -1, respectively.
		if Ichips[int(Iindex)] == 1:
			Ibit = 1
		else:
			Ibit = -1

		if Qchips[int(Qindex)] == 1:
			Qbit = 1
		else:
			Qbit = -1

		Iphase = Ibit*np.sin(x - Iindex * np.pi)
		Qphase = Qbit*np.sin(x - np.pi / 2 - Qindex * np.pi)

		if abs(Iphase) < 0.01:
			Iphase = 0
		if abs(Qphase) < 0.01:
			Qphase = 0

		Iphase = Iphase * ratio
		Qphase = Qphase * ratio
		IQs.append([Iphase,Qphase])

	return IQs


def shiftCons(approxfftIQs, xdelta, ydelta):
	length = len(approxfftIQs)
	newIQs = [[0, 0] for i in range(length)]
	for i in range(length):
		newIQs[i][0] = approxfftIQs[i][0] + xdelta
		newIQs[i][1] = approxfftIQs[i][1] + ydelta

	return newIQs

# def cyclicprefixer(IQs):
# 	M = WITHOUTCP
# 	N = WITHCP
# 	length = len(IQs)
# 	G = int(length / M)
# 	afterIQs = []
# 	for g in range(G):
# 		tempIQs = IQs[g*M:(g+1)*M]
		
# 	for j in range(N-M):
# 		cpI = tempIQs[2*M - N + j][0]
# 		cpQ = tempIQs[2*M - N + j][1]
# 		afterIQs.append([cpI,cpQ])
		
# 	for j in range(M):
# 		afterIQs.append([tempIQs[j][0], tempIQs[j][1]])
		
# 	return afterIQs

def get_freq_qam(symbol, ratio, subc):

	return 0

def generateIQSignal(symbol, ratio, subc):
	WholeSignals = []
	WholeIQs = generateIQ(symbol, ratio)
	length = len(WholeIQs)
	groupsize = WITHCP
	groupbeforecp = WITHOUTCP
	delta = WITHCP - WITHOUTCP
	Groups = int(length / groupsize)
	approxIY = []
	approxQY = []
	for g in range(Groups):
		complexIQs = []
		for i in range(delta, WITHCP):
			complexIQs.append(np.complex64(WholeIQs[g*groupsize + i][0] + WholeIQs[g*groupsize + i][1]*1j))

		complexIQs = np.array(complexIQs)
		tempsp = np.fft.fft(complexIQs)
		tempsp = np.fft.fftshift(tempsp)

		tempSignals = []
		for jj in range(len(tempsp)):
			tempSignals.append([tempsp[jj].real, tempsp[jj].imag])


		for jj in range(len(tempsp)):
			comangle = 1*2*jj*np.pi*delta*1.0/WITHOUTCP
			offset = [np.cos(comangle), np.sin(comangle)]
			tempSignals.append(complexmultiply([tempsp[jj].real, tempsp[jj].imag],offset))

		sp = []
		for jj in range(groupbeforecp):
			sp.append(np.complex64(tempSignals[jj][0]+tempSignals[jj][1]*1j))

		for i in range(len(sp)):
			if abs(i - len(sp)/2) <= subc:
				WholeSignals.append([sp[i].real,sp[i].imag])
			else:
				WholeSignals.append([0,0])

	return WholeSignals


def generateISignal(symbol, ratio, subc):
	WholeSignals = []
	WholeIQs = generateIQ(symbol, ratio)
	length = len(WholeIQs)
	groupsize = WITHCP
	groupbeforecp = WITHOUTCP
	delta = WITHCP - WITHOUTCP
	Groups = int(length / groupsize)
	approxIY = []
	approxQY = []
	for g in range(Groups):
		complexIQs = []
		if ALIGNMENT[g] == 0:
			for i in range(0, WITHCP - delta):
				complexIQs.append(np.complex64(WholeIQs[g*groupsize + i][0]))

		if ALIGNMENT[g] == 1:
			for i in range(delta, WITHCP):
				complexIQs.append(np.complex64(WholeIQs[g*groupsize + i][0]))

		complexIQs = np.array(complexIQs)
		tempsp = np.fft.fft(complexIQs)
		tempsp = np.fft.fftshift(tempsp) 
		# See: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#/media/File:Fourier_transform,_Fourier_series,_DTFT,_DFT.svg

		
		tempSignals = []
		if ALIGNMENT[g] == 0:
			# Multiply I and Q ZigBee signal with cos-wave and sin-wave, respectively.
			# Effect: Time domain inphase signal is cyclic left shifted 0.8us.
			for jj in range(len(tempsp)):
				comangle = 1*2*jj*np.pi*delta*1.0/WITHOUTCP
				offset = [np.cos(comangle), np.sin(comangle)]
				tempSignals.append(complexmultiply([tempsp[jj].real, tempsp[jj].imag], offset))

		if ALIGNMENT[g] == 1:
			for jj in range(len(tempsp)):
				tempSignals.append([tempsp[jj].real, tempsp[jj].imag])

		sp = []
		for jj in range(groupbeforecp):
			sp.append(np.complex64(tempSignals[jj][0]+tempSignals[jj][1]*1j))

		# Only allows "2*subc" subcarriers in the middle of the spectum to be non-zero.
		for i in range(len(sp)):
			if abs(i - len(sp)/2) <= subc:
				WholeSignals.append([sp[i].real,sp[i].imag])
			else:
				WholeSignals.append([0,0])

	# WoleSignals is a complex signal, instead of, suggested by the function name, a in-phase signal.
	# Donnot know why.
	return WholeSignals


def generateQSignal(symbol, ratio, subc):
	WholeSignals = []
	WholeIQs = generateIQ(symbol, ratio)
	length = len(WholeIQs)
	groupsize = WITHCP
	groupbeforecp = WITHOUTCP
	delta = WITHCP - WITHOUTCP
	Groups = int(length / groupsize)
	delta1 = int(WITHCP / 8 )
	# "deltal" is the duration of a half-chip_seq signal flipping. 
	# See the "solution" in section 5.2 of the original paper.

	LeftWholeIQs = copy.deepcopy(WholeIQs)
	RightWholeIQs = copy.deepcopy(WholeIQs)

	LeftWholeIQs = cycleshift(RightWholeIQs, -1*delta1)
	RightWholeIQs = cycleshift(RightWholeIQs, delta1)

	#Left Shift
	if symbol == 5 or symbol == 2 or symbol == 4 or symbol == 6 or symbol == 12 or symbol == 13 or symbol == 14:
		for i in range(delta1):
			LeftWholeIQs[i][1] = -1 * LeftWholeIQs[i][1] # type: ignore

    #Right Shift
	if symbol == 5 or symbol == 2 or symbol == 4 or symbol == 6 or symbol == 12 or symbol == 13 or symbol == 14:
		for i in range(delta1):
			RightWholeIQs[length - delta1 + i][1] = -1 * RightWholeIQs[length - delta1 + i][1] # type: ignore

	approxIY = []
	approxQY = []
	for g in range(Groups):
		complexIQs = []
		if ALIGNMENT[g] == 0:
			for i in range(0, WITHCP - delta):
				complexIQs.append(np.complex64(LeftWholeIQs[g*groupsize + i][1]*1j)) # type: ignore
		if ALIGNMENT[g] == 1:
			for i in range(delta, WITHCP):
				complexIQs.append(np.complex64(RightWholeIQs[g*groupsize + i][1]*1j)) # type: ignore

		complexIQs = np.array(complexIQs)
		tempsp = np.fft.fft(complexIQs)
		tempsp = np.fft.fftshift(tempsp)

		# Cyclic left shift quadrature signal for 0.5us + 0.8us if ALIGMENT is 0.
		# Cyclic right shift quadrature signal for 0.5us if ALIGMENT is 1.
		tempSignals = []
		for jj in range(len(tempsp)):
			comangle = 0
			if ALIGNMENT[g] == 0:
				comangle = 1*2*jj*np.pi*delta1*1.0/WITHOUTCP
			if ALIGNMENT[g] == 1:
				comangle = -1*2*jj*np.pi*delta1*1.0/WITHOUTCP

			offset = [np.cos(comangle), np.sin(comangle)]
			temp = complexmultiply([tempsp[jj].real, tempsp[jj].imag],offset)

			if ALIGNMENT[g] == 0:
				comangle = 1*2*jj*np.pi*delta*1.0/WITHOUTCP
			if ALIGNMENT[g] == 1:
				comangle = 0
			offset = [np.cos(comangle), np.sin(comangle)]
			tempSignals.append(complexmultiply(temp,offset))

		sp = []
		for jj in range(groupbeforecp):
			sp.append(np.complex64(tempSignals[jj][0]+tempSignals[jj][1]*1j))

		for i in range(len(sp)):
			if abs(i - len(sp)/2) <= subc:
				WholeSignals.append([sp[i].real,sp[i].imag])
			else:
				WholeSignals.append([0,0])
	return WholeSignals


def generateIFFTSignal(approxfftIQs):
	length = len(approxfftIQs)
	newIQs = []
	complexfftIQs = []
	groupsize = WITHOUTCP
	Groups = int(length / groupsize)
	for g in range(Groups):
		fftIQs = [[0, 0] for i in range(groupsize)]
		for i in range(groupsize):
			fftIQs[i][0] = approxfftIQs[g*groupsize + i][0]
			fftIQs[i][1] = approxfftIQs[g*groupsize + i][1]

		length = len(fftIQs)
		complexfftIQs = []
		for i in range(length):
			Isample = fftIQs[i][0]
			Qsample = fftIQs[i][1]
			complexfftIQs.append(np.complex64(Isample+Qsample*1j))

		#print complexfftIQs
		complexfftIQs = np.array(complexfftIQs)
		complexfftIQs = np.fft.fftshift(complexfftIQs)
		sp = np.fft.ifft(complexfftIQs)


		#spec = np.sqrt(sp.real*sp.real + sp.imag*sp.imag)
		for i in range(len(sp)):
			newIQs.append([sp[i].real,sp[i].imag])

	return newIQs




def compensateCFO(IQs, deltaf, Fs):
	N = len(IQs)
	newIQs = [[0,0] for i in range(N)]
	for i in range(N):
		comangle = float(2*np.pi*deltaf*i) / Fs
		offset = [np.cos(comangle), np.sin(comangle)]
		sample = [IQs[i][0], IQs[i][1]]
		com = complexmultiply(sample,offset)
		newIQs[i][0] = com[0]
		newIQs[i][1] = com[1]
	return newIQs

def generateMultipleSignal(count, chip_map, ratio):
	Groups = 5
	GroupSize = 48
	subc = 5
	starti = [0,13,25,38]
	channels = [-21,-7,7,21]
	WholeSignals = [[0,0] for i in range(Groups * GroupSize)]
	comdeltaf = [-0.4375,0.1875,0.8125,1.4375]
	#comdeltaf = [0,0,0,0]
	for cc in range(count):
		chip_seq = chip_map[cc]
		#subchannels = channels[cc]
		WholeIQs = generateIQ(chip_seq, ratio)
		deltaf = comdeltaf[cc]
		WholeIQs = compensateCFO(WholeIQs,deltaf,15)
		approxIY = []
		approxQY = []
		for g in range(Groups):
			IQs = [[0, 0] for i in range(GroupSize)]
			for i in range(GroupSize):
				IQs[i][0] = WholeIQs[g*GroupSize + i][0]
				IQs[i][1] = WholeIQs[g*GroupSize + i][1]

			length = len(IQs)
			complexIQs = []
			for i in range(length):
				Isample = IQs[i][0]
				Qsample = IQs[i][1]
				complexIQs.append(np.complex64(Isample+Qsample*1j))


			complexIQs = np.array(complexIQs)
			sp = np.fft.fft(complexIQs)
			#spec = np.sqrt(sp.real*sp.real + sp.imag*sp.imag)

			fftIQs = []
			for i in range(len(sp)):
				fftIQs.append([sp[i].real,sp[i].imag])


			qam = generateQAM(64,1.0801)
			qam.append([0,0])

			approxfftIQs = approxQAM(qam, fftIQs)

			for i in range(subc):
				loc = starti[cc] + subc + i
				WholeSignals[g*GroupSize + loc][0] = approxfftIQs[i+1][0]
				WholeSignals[g*GroupSize + loc][1] = approxfftIQs[i+1][1]
			for i in range(subc):
				loc = starti[cc] + i
				WholeSignals[g*GroupSize + loc][0] = approxfftIQs[GroupSize - subc + i - 1][0]
				WholeSignals[g*GroupSize + loc][1] = approxfftIQs[GroupSize - subc + i - 1][1]

	return WholeSignals