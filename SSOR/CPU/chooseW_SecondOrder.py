__author__ = 'Administrator'
import math
import scipy.io as sio
import scipy as sc
import numpy as np
import numpy.linalg as npli
import time
import scipy.sparse as ssp
from scipy import *

def SSOR(filename):

    A = ssp.csr_matrix(sio.mmread(filename))
    print("the row is %d" % A.shape[0])
    # print("===================A")
    # print(A)


    # cond_A = npli.cond(A.toarray())
    # print("the cond of A is %3f " % cond_A)


    D = (A.diagonal())

    D = ssp.diags(A.diagonal(), 0)


    # print("===================D")
    # print(D)


    w = 0.5
    while(w <= 1.5):
        spectral = 0

        print("the w is %f" % w)

        D_like = (1/w)*D

        L = ssp.tril(A)
        #
        # print("===================L")
        # print(L)



        D_inv = ssp.csr_matrix(npli.inv(D_like.toarray()))
        # print("===================D_inv")
        # print(D_inv)

        tmp = D_inv*L

        # print("===================tmp")
        # print(tmp)

        spectral = max(abs(sc.linalg.eigvals(tmp.toarray())))
        # spectral = sc.linalg.eigvals(tmp.toarray())

        print("the spectral is %3f" % spectral)
        # print("===================D_inv")
        # print(D_inv)

        D_inv2 = ssp.csr_matrix(D_inv.toarray()**0.5)

        # print("===================D_inv2")
        # print(D_inv2)

        K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv + D_inv2*L*D_inv*L*D_inv)

        # print("===================K")
        # print(ssp.csr_matrix(D_inv2*L*D_inv))

        K_T = K.T

        # M = K_T*K
        # print("===================M")
        # print(M)

        b = np.random.rand(A.shape[0])

        #need to change a #######################&&&&&&&&&
        a = 0.5

        A_2 = A.toarray()
        K_2 = K.toarray()

        print("when a=%f," % a)
        for i in range(0, A_2.shape[0]):
            maxA = abs(A_2[i][0])
            for j in range(0, A_2.shape[0]):
                if abs(A_2[i][j]) > maxA:
                    maxA = abs(A_2[i][j])

            if abs(A_2[i][j]) <= (1-a)*maxA:
                K_2[i][j] = 0

        K = ssp.csr_matrix(K_2)

        M = K.T*K

        M_2 = M.toarray()


        for i in range(0, A_2.shape[0]):
            for j in range(0, A_2.shape[0]):
                if K_2[i][j] == 0:
                    M_2[i][j] = 0

        A_like = A*ssp.csr_matrix(M_2)

        x, info = linalg.cg(A_like, b)

        print("the iteration times of A*M is %d " % info)

        w += 0.1


        # print("===================A*M")
        # print(A_like)



print("input file name:\n")

# SSOR("Trefethen_2000.mtx")
# SSOR("685_bus.mtx")
# SSOR("plbuckle.mtx")
# SSOR("bcsstm26.mtx")
# SSOR("plat1919.mtx")
# SSOR("Trefethen_700.mtx")
# SSOR("bcsstk10.mtx")
# SSOR("msc00726.mtx")
# SSOR("plbuckle.mtx")


########## more ##########
# SSOR("bcsstk13.mtx")
# SSOR("bcsstk15.mtx")
# SSOR("bcsstk24.mtx")
# SSOR("nasa2910.mtx")
# SSOR("nasa4707.mtx")
# SSOR("sts4098.mtx")