#calculate the source bits from the QAM points

import numpy as np
import copy
import random
import struct
from CRC16Kermit import CRC16Kermit
import numpy as np
# from GF import GF

# gf = GF(2)


def bitxor(a,b):
	if a == b:
		return 0
	else:
		return 1

def XOR(bits):
	# This function computes XOR of bits in a given sequence named "bits".
	N = len(bits)
	resbit = 0
	for i in range(N):
		resbit = bitxor(resbit,bits[i])

	return resbit

def getOtherSourceBit(index):
	return 0

def getOtherSourceBits(indexes, CERegisters):
	# This function apply the effect of coding of last symbol to the correct symbol.
	L = len(indexes)
	res = [0 for i in range(L)]
	for i in range(L):
		index = indexes[i]
		if index <= 0:
			# The ith coded bit (in "Y[i]") is effected by bits of the last symbol, of which
			# copies are left in CERegisters.
			index1 = int(-1*index)
			res[i] = CERegisters[index1]
		else:
			res[i] = getOtherSourceBit(index)

	return res

def generalizedEuclidianAlgorithm(a, b):
	if b > a:
		#print a, b
		return generalizedEuclidianAlgorithm(b,a)
	elif b == 0:
		return (1, 0)
	else:
        #print a,b
		(x, y) = generalizedEuclidianAlgorithm(b, a % b)
		return (y, x - (a / b) * y)

def inversemodp(a, p):
	a = a % p
	if (a == 0):
		#print "a is 0 mod p"
		return 0
	(x,y) = generalizedEuclidianAlgorithm(p, a % p)
	return y % p

# def identitymatrix(n):
# 	return [[long(x == y) for x in range(0, n)] for y in range(0, n)]

def inversematrix(A):
	# This function calculate the inverse of a given matrix, A, via Gaussian
	# elimination method.
	# https://en.wikipedia.org/wiki/Invertible_matrix#Gaussian_elimination
	N = len(A)

	# Generate a N*N identity matrix.
	B = [[0 for i in range(N)] for j in range(N)]
	for i in range(N):
		B[i][i] = 1

	for K in range (0, N):
		if (A[K][K] == 0):
			flag = 1
			i = K + 1
			while(flag):
				if (A[i][K]!=0):
					# Swap the K-th row and the i-th row.
					for L in range (0, N):
						s = A[K][L]
						A[K][L] = A[i][L]
						A[i][L] = s
						s = B[K][L]
						B[K][L] = B[i][L]
						B[i][L] = s
					flag = 0
				else:
					i = i + 1
					if i == N:
						print("NOT INVERSEABLE")

		for I in range (0, N):
			if (I!=K):
				if (A[I][K]):
					for M in range (0, N):
						# "^" means XOR operator
						# The reason of adopting XOR operator is that the 
						# matrix is in GF(2) field (according to the WEBee paper).
						A[I][M] = A[I][M] ^ A[K][M]
						B[I][M] = B[I][M] ^ B[K][M]

	return B

# def findinvert(a, p):
# 	if not p == 2:
# 		print "ERROR"
# 	else:
# 		return a

# def identitymatrix(n):
#     return [[long(x == y) for x in range(0, n)] for y in range(0, n)]

# def multiply_vector_scalar (vector, scalar, q):
#     kq = []
#     for i in range (0, len(vector)):
#         kq.append (vector[i] * scalar %q)
#     return kq

# def minus_vector_scalar(vector1, scalar, vector2, q):
#     kq = []
#     for i in range (0, len(vector1)):
#         kq.append ((vector1[i] - scalar * vector2[i]) %q)
#     return kq

# def inversematrix(matrix, q):
# 	n = len(matrix)
# 	A_invert = [[0 for i in range(n)] for j in range(n)]
# 	A =[]
# 	for j in range (0, n):
# 		temp = []
# 		for i in range (0, n):
# 			temp.append (matrix[j][i])
# 		A.append(temp)

# 	Ainv = identitymatrix(n)

# 	for i in range(0, n):
# 		factor = findinvert(A[i][i], q) #invert mod q
# 		A[i] = multiply_vector_scalar(A[i],factor,q)
# 		Ainv[i] = multiply_vector_scalar(Ainv[i],factor,q)
# 		for j in range(0, n):
# 			if (i != j):
# 				factor = A[j][i]
# 				A[j] = minus_vector_scalar(A[j], factor, A[i], q)
# 				Ainv[j] = minus_vector_scalar(Ainv[j], factor, Ainv[i], q)

# 	for i in range(n):
# 		for j in range(n):
# 			A_invert[i][j] = int(Ainv[i][j])
# 	return A_invert

def SolveXOREquations(A, Y):
	# This function multiplies a given matrix, A, to a given 
	# vector, Y, in GF(2) field.
	N = len(Y)
	res = [0 for i in range(N)]
	for i in range(N):
		temp = []
		for j in range(N):
			temp.append(A[i][j] & Y[j])
		res[i] = XOR(temp)
	return res

def MultiplyXORMatrix(A, B):
	N = len(A)
	res = [[0 for i in range(N)] for j in range(N)]
	for i in range(N):
		for j in range(N):
			temp = []
			for k in range(N):
				temp.append(A[i][k] & B[k][j])
			res[i][j] = XOR(temp)
	return res

#generate converlutional encoding mapping table
NCBPS = 288 # Number of coded bits per symbol

T1 = [0, -2, -3, -5, -6]
T2 = [0, -1, -2, -3, -6]

CETable = [[] for i in range(NCBPS)] # Elements in the table are indices of uncoded bits.

totalY = 0
for i in range(int(NCBPS/3 + 1)):
	tempY = [[0 for k in range(5)] for l in range(6)]
	for j in range(3):
		t = i * 3 + j
		for k in range(5):
			tempY[j*2+0][k] = t + T1[k] + 1
			tempY[j*2+1][k] = t + T2[k] + 1
	
	# The idea of the following code: puncturing in convolutional encoding with 3/4 rate.
	# Puncturing here: for every 6 bits coded with 1/2 rate, omit the 4-th and 5-th bits.
	if totalY < NCBPS:
		for k in range(5):
			CETable[totalY].append(tempY[0][k])

	totalY = totalY + 1

	if totalY < NCBPS:
		for k in range(5):
			CETable[totalY].append(tempY[1][k])

	totalY = totalY + 1

	if totalY < NCBPS:
		for k in range(5):
			CETable[totalY].append(tempY[2][k])

	totalY = totalY + 1

	if totalY < NCBPS:
		for k in range(5):
			CETable[totalY].append(tempY[5][k])

	totalY = totalY + 1

# Generate interleving mapping (table).
NCBPS = 288 # Number of coded bits per symbol
s = 3 # =max(NCBPS/2, 1)
# Indices of bits after interleaving. 
# Here, the k-th bit will be move to the position where the j-th bit at.
permutation = [0 for i in range(NCBPS)]
for k in range(NCBPS):
	i = int((NCBPS / 16) * (k % 16) + int(k / 16))
	j = s * int(i / s) + (i + NCBPS - int(((16 * i) / NCBPS))) % s
	permutation[k] = j # type: ignore
# Retrace the original indices of each permutated/interleaved bits. 
# For example, it will show that the j-th permutated bit is the original k-th bit.
reverseperm = [0 for i in range(NCBPS)]
for k in range(NCBPS):
	reverseperm[permutation[k]] = k

qaminteger = [-7,-5,-3,-1,1,3,5,7]
qamdict = {}
delta = float(1.0801 * 2) / 7
for i in range(8):
	key = round(-1.0801 + delta * i, 5)
	value = qaminteger[i]
	qamdict.update({key:value})
qamdict.update({0.0: 0})
qamtable = {-7:[0,0,0], -5:[0,0,1], -3:[0,1,1], -1:[0,1,0], 1:[1,1,0], 3:[1,1,1],5:[1,0,1],7:[1,0,0]}
qamtable.update({0: [2,2,2]})


for qamkk in range(100):
	QAMPoints = 48 # Number of QAM points, which is equal to number of WiFi data subcarriers.
	alloutputbits = [] # All WiFi bits that will be imputted to WiFi modulator to emulated a ZigBee packet.
	filename = "./data/WEBeeQAMs"+str(qamkk)+".txt"
	qamFile = open(filename, "r")
	index = 0
	outputbits = [] # WiFi bits per WiFi data symbol that will be used to emulated a ZigBee packet.
	# Demodulate each of 48 subcarriers
	for line in qamFile.readlines():
		strline = line.split(',')
		real = float(strline[0])
		img = float(strline[1])
		if (index == QAMPoints):
			alloutputbits.append(copy.deepcopy(outputbits))
			outputbits = []
			index = 0

		realindex = qamdict[real]
		imgindex = qamdict[img]
		for bit in qamtable[realindex]:
			outputbits.append(bit)
		for bit in qamtable[imgindex]:
			outputbits.append(bit)
		index = index + 1

	if (index == QAMPoints):
		alloutputbits.append(copy.deepcopy(outputbits))
		outputbits = []
		index = 0

	NWEBeeSybols = len(alloutputbits) # Number of WiFi data symbols.

	#build the bigraph
	# For each WiFi symbol, the index of bits that will be modulated as emulated 
	# ZigBee signal. Here the bits are assumed to be permutated.
	ReverseIndex = []
	symbols = len(alloutputbits[0])
	for i in range(symbols):
		if not alloutputbits[0][i] == 2:
			ReverseIndex.append(i)
	# The original index of each bit that will be modulated as emulated ZigBee
	# signal before it is permutated.
	Z = []
	for i in range(len(ReverseIndex)):
		Z.append(reverseperm[ReverseIndex[i]])

	# Correct Here
	# print ReverseIndex  

	#Begin to find a Matrix invertable

	# M = len(Z)
	# XS = []
	# G = [[] for i in range(M)]
	# for i in range(M):
	# 	for k in range(len(CETable[Z[i]])):
	# 		if CETable[Z[i]][k] > 0:
	# 			G[i].append(CETable[Z[i]][k])
	# 			if not CETable[Z[i]][k] in XS:
	# 				XS.append(CETable[Z[i]][k])

	# # XS correct

	# Rows = M
	# Left = []
	# Right = []
	# for i in range(Rows):
	# 	Left.append([])
	# 	Right.append(copy.deepcopy(CETable[Z[i]]))

	# for i in range(len(XS)):
	# 	x = XS[i]
	# 	for j in range(Rows):
	# 		if x in Right[j]:
	# 			Right[j].remove(x)
	# 			if not x in Left[j]:
	# 				Left[j].append(x)

	# XSSize = len(XS)
	# Xdict = {}
	# for i in range(XSSize):
	# 	Xdict.update({XS[i]:i})

	# MatrixRows = Rows
	# MatrixCols = XSSize

	# MMatrix = [[0 for i in range(MatrixCols)] for j in range(MatrixRows)]

	# for i in range(MatrixRows):
	# 	for key in Left[i]:
	# 		MMatrix[i][Xdict[key]] = 1


	# A = MMatrix
	# AM = np.matrix(A)
	# X = XS
	# M = gf.matrix_rank(AM)
	# N = len(A[0])
	# cols = N
	# for i in range(N-M):
	# 	flag = 1
	# 	while(flag):
	# 		j = random.randint(0,cols-1)
	# 		B = np.delete(AM, j, axis = 1)
	# 		if gf.matrix_rank(B) == M:
	# 			X.pop(j)
	# 			AM = B
	# 			cols = cols - 1
	# 			flag = 0
	# print X
	# print len(X)

	#Correct Here, We GET XS which can control the QAM points

	# XS = [85, 83, 82, 79, 97, 95, 94, 73, 71, 121, 120, 119, 115, 133, 131, 130, 127, 109, 108, 107, 106, 103, 157, 156, 151, 169, 167, 166, 145, 144, 142, 192, 191, 187, 205, 204, 202, 181, 180, 179, 26, 24, 21, 20, 2, 14, 12, 11, 9, 62, 60, 57, 35, 33, 50, 47, 45, 98, 92, 74, 69, 68, 86, 80, 129, 128, 105, 116, 170, 165, 141, 140, 158, 153, 206, 201, 200, 182, 176, 194, 189, 188, 3, 1, 15, 27, 25, 39, 49, 61, 75, 87, 99, 111, 123, 135, 147, 159, 195, 10, 28, 64, 58, 34, 88, 76, 124, 136, 160, 148, 184, 30, 6, 19, 43, 55, 54, 90, 138, 114, 126, 174, 213, 212, 5, 29, 41, 77, 89, 113, 173, 214]

	XS = [85, 84, 83, 82, 96, 73, 72, 70, 120, 119, 118, 115, 131, 130, 127, 109, 107, 106, 155, 154, 169, 166, 163, 144, 142, 139, 190, 187, 205, 181, 180, 178, 26, 23, 2, 14, 9, 59, 57, 56, 35, 32, 50, 48, 47, 45, 98, 93, 92, 69, 68, 81, 134, 128, 110, 105, 104, 122, 117, 116, 170, 146, 140, 158, 206, 201, 200, 176, 194, 189, 188, 3, 1, 13, 25, 39, 49, 63, 75, 87, 99, 135, 159, 171, 195, 16, 10, 28, 22, 4, 52, 46, 64, 58, 34, 88, 76, 124, 112, 160, 148, 196, 208, 184, 31, 30, 6, 18, 43, 54, 78, 138, 114, 174, 150, 162, 212, 210, 198, 17, 29, 65, 77, 89, 101, 125, 137, 149, 161, 173, 197, 209]

	M = len(Z)
	Rows = M
	Left = []
	Right = []
	# The i-th row of "Right" indicates the indices of uncoded bits that will be encoded as the i-th coded bit.
	for i in range(Rows):
		Left.append([])
		Right.append(copy.deepcopy(CETable[Z[i]]))
	

	for i in range(len(XS)):
		x = XS[i]
		for j in range(Rows):
			if x in Right[j]:
				Right[j].remove(x)
				if not x in Left[j]:
					Left[j].append(x)

	for i in range(Rows):
		Right[i] = [Right[i],[Z[i]]]

	XSSize = len(XS)
	Xdict = {}
	for i in range(XSSize):
		Xdict.update({XS[i]:i})

	MatrixRows = Rows
	MatrixCols = XSSize

	Matrix = [[0 for i in range(MatrixCols)] for j in range(MatrixRows)]

	for i in range(MatrixRows):
		for key in Left[i]:
			Matrix[i][Xdict[key]] = 1


	#print gf.matrix_rank(np.matrix(Matrix))
	tempMatrix = copy.deepcopy(Matrix)
	A_invert = inversematrix(tempMatrix)
	#print A_invert
	#print gf.matrix_rank(np.matrix(A_invert))

	# # res = MultiplyXORMatrix(Matrix, A_invert)
	# # N = len(res)
	# # for i in range(N):
	# # 	if not res[i][i] == 1:
	# # 		print 'NO INVERSE'


	NWEBeeSymbols = len(alloutputbits) # Number of WiFi data symbols.

	WEBeeBits = int(NCBPS * 3 / 4) # Number of bits before convolutional encoding with a coding rate of 3/4.

	CERegisters = [0,0,0,0,0,0,0]
	packetfile = "./data/WEBeeBits"+str(qamkk)+".txt"
	sourceFile = open(packetfile, "w")
	#print NWEBeeSymbols
	# Operations for each WiFi data symbol.
	for sym in range(NWEBeeSymbols):
		Y = [0 for i in range(M)]
		temp = []
		SourceBits = [0 for i in range(WEBeeBits)]
		for i in range(M):
			# The "correctbit" below is the ith bit that will be fed into mudulator to 
			# emulate ZigBee signals.
			correctbit = alloutputbits[sym][permutation[Z[i]]] 
			# The Y[i] below is the i-th coded bit.
			Y[i] = XOR(getOtherSourceBits(Right[i][0], CERegisters)+[correctbit])
		
		TY = [0 for i in range(WEBeeBits)] 

		# Reverse the convolutional encoding: find uncoded bits corresponding to coded bits 
		# that used to emulate ZigBee signals.
		CX = SolveXOREquations(A_invert, Y)
		# Y1 = SolveXOREquations(Matrix, CX)

		# for ii in range(len(Y)):
		# 	if not Y[ii] == Y1[ii]:
		# 		print "EEEEEEE"
		for i in range(WEBeeBits):
			if (i+1) in XS:
				cx = CX[Xdict[i+1]]
				SourceBits[i] = cx
			else:
				SourceBits[i] = getOtherSourceBit(i) # Set it to 0

		for i in range(len(CERegisters)):
			CERegisters[i] = SourceBits[WEBeeBits - i - 1]
		
		for b in range(len(SourceBits)):
			sourceFile.write(str(SourceBits[b]))
			if b < len(SourceBits) - 1:
				sourceFile.write(",")

		sourceFile.write("\n")

	sourceFile.close()
