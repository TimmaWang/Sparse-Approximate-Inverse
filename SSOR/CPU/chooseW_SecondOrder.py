__author__ = 'Administrator'
import math
import scipy.io as sio

import scipy.sparse as ssp
from numpy import *
from scipy import *
from scipy.sparse import linalg

iterations = 0

def report(xk):
    global iterations
    iterations += 1
    if iterations > 10000:
        print("iteration times over 10000.")
    # print(iterations)

def SSOR(filename):

    A = ssp.csr_matrix(sio.mmread(filename))
    print("the row is %d" % A.shape[0])
    # print("===================A")
    # print(A)

    D = ssp.diags(A.diagonal(), 0)


    # print("===================D")
    # print(D)


    w = 0.5
    while(w <= 1.5):

        print("the w is %f" % w)

        D_like = (1/w)*D

        L = ssp.tril(A)
        #
        # print("===================L")
        # print(L)

        D_inv = D_like
        for i in range(0, len(D_inv.data)):
            D_inv.data[i] = 1/D_inv.data[i]
        # print("===================D_inv")
        # print(D_inv)

        tmp = D_inv*L

        # print("===================tmp")
        # print(tmp)


        D_inv2 = D_inv
        for i in range(0, len(D_inv.data)):
            D_inv2.data[i] **= 0.5

        # print("===================D_inv2")
        # print(D_inv2)

        K = math.sqrt(2-w)*(D_inv2 - D_inv2*L*D_inv + D_inv2*L*D_inv*L*D_inv)

        # print("===================K")
        # print(ssp.csr_matrix(D_inv2*L*D_inv))

        K_T = K.T

        # M = K_T*K
        # print("===================M")
        # print(M)

        b = random.rand(A.shape[0])

        #need to change a #######################&&&&&&&&&
        a = 0.5


        print("when a=%f," % a)

        for row in range(0, len(K.indptr)-1):
            maxA = abs(A.data[A.indptr[row]])

            for indexA in range(A.indptr[row], A.indptr[row+1]):
                if abs(A.data[indexA]) > maxA:
                    maxA = abs(A.data[indexA])

            for indexM in range(A.indptr[row], A.indptr[row+1]):
                if abs(A.data[indexM]) <= (1-a)*maxA:
                    K.data[indexM] = 0

        print("finish calculating K.")

        M = K.T*K


        for row in range(0, len(K.indptr)-1):
            for indexM in range(K.indptr[row], K.indptr[row+1]):
                if K.data[indexM] == 0:
                    M.data[indexM] = 0

        print("finish calculating M.")

        A_like = A*M

        print("finish calculating A_like.")

        global iterations

        x, info = linalg.bicgstab(A_like, b, callback=report)

        print("the iteration times of A_like is %d " % iterations)

        w += 0.1



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