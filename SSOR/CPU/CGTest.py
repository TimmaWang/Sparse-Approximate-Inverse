__author__ = 'Administrator'
import math
import datetime
import scipy.io as sio
import scipy as sc
import numpy as np
import numpy.linalg as npli

import scipy.sparse as ssp

from scipy.sparse import linalg


def SSOR(filename):


    fileinfo = sio.mminfo(filename)
    A = ssp.csr_matrix(sio.mmread(filename))
    print(fileinfo.__str__())

    b = np.random.rand(A.shape[0])



    starttime = datetime.datetime.now()

    x_direct, info_direct = linalg.cg(A, b)

    endtime_d = datetime.datetime.now()

    print("the directing time cost is %f" % ((endtime_d - starttime).microseconds/1000))
    print("the iterations is %d" % info_direct)




    D = (A.diagonal())

    D = ssp.diags(A.diagonal(), 0)


    # print("===================D")
    # print(D)


    w = 0.5

    D_like = (1/w)*D

    L = ssp.tril(A)


    D_inv = ssp.csr_matrix(npli.inv(D_like.toarray()))
    # print("===================D_inv")
    # print(D_inv)

    tmp = D_inv*L


    D_inv2 = ssp.csr_matrix(D_inv.toarray()**0.5)


    # one
    K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv)

    # two
    # K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv + D_inv2*L*D_inv*L*D_inv)


    K_T = K.T
    M = K_T*K

    a = 0.5


    A_2 = A.toarray()
    M_2 = M.toarray()

    print("when a=%f," % a)
    for i in range(0, A_2.shape[0]):
        maxA = abs(A_2[i][0])
        for j in range(0, A_2.shape[0]):
            if abs(A_2[i][j]) > maxA:
                maxA = abs(A_2[i][j])

        if abs(A_2[i][j]) <= (1-a)*maxA:
            M_2[i][j] = 0

    A_like = A*ssp.csr_matrix(M_2)

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




