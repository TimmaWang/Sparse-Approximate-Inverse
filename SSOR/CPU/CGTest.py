__author__ = 'Administrator'
import math
import datetime
import math
import scipy.io as sio

import scipy.sparse as ssp
from numpy import *
from scipy import *
from scipy.sparse import linalg


def SSOR(filename):


    fileinfo = sio.mminfo(filename)
    A = ssp.csr_matrix(sio.mmread(filename))
    print(fileinfo.__str__())

    b = random.rand(A.shape[0])



    starttime = datetime.datetime.now()

    x_direct, info_direct = linalg.cg(A, b)

    endtime_d = datetime.datetime.now()

    print("the directing time cost is %f" % ((endtime_d - starttime).microseconds/1000))
    print("the iterations is %d" % info_direct)


    D = ssp.diags(A.diagonal(), 0)


    # print("===================D")
    # print(D)


    w = 0.5

    D_like = (1/w)*D

    L = ssp.tril(A)


    D_inv = D_like
    for i in range(0, len(D_inv.data)):
        D_inv.data[i] = 1/D_inv.data[i]
    # print("===================D_inv")
    # print(D_inv)

    tmp = D_inv*L


    D_inv2 = D_inv
    for i in range(0, len(D_inv.data)):
        D_inv2.data[i] **= 0.5


    # one
    K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv)

    # two
    # K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv + D_inv2*L*D_inv*L*D_inv)


    K_T = K.T
    M = K_T*K

    a = 0.5

    print("when a=%f," % a)
    for row in range(0, len(M.indptr)-1):
            maxA = abs(A.data[A.indptr[row]])

            for indexA in range(A.indptr[row], A.indptr[row+1]):
                if abs(A.data[indexA]) > maxA:
                    maxA = abs(A.data[indexA])

            for indexM in range(M.indptr[row], M.indptr[row+1]):
                if abs(M.data[indexM]) <= (1-a)*maxA:
                    M.data[indexM] = 0

    A_like = A*M

    x, info = linalg.cg(A_like, b)


    starttime_p = datetime.datetime.now()

    x, info = linalg.cg(A_like, b)


    endtime = datetime.datetime.now()

    print("the SSOR time cost is %f" % (endtime - starttime_p).seconds)
    print("the iterations is %d" % info)


    # print("===================A*M")
    # print(A_like)

#filename = raw_input()

SSOR("685_bus.mtx")
SSOR("plbuckle.mtx")
SSOR("bcsstk10.mtx")
SSOR("msc00726.mtx")
SSOR("nasa2910.mtx")
SSOR("sts4098.mtx")
# SSOR("nasa4704.mtx")




