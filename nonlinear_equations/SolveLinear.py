# SolveLinear.py
# Python module for PHY407
# Paul Kushner, 2015-09-26
# Modifications by Nicolas Grisouard, 2018-09-26
# This module contains useful routines for solving linear systems of equations.
# Based on gausselim.py from Newman
# from numpy import empty
# The following will be useful for partial pivoting
from numpy import empty, copy
import numpy as np


#


def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m + 1, N):
            mult = A[i, m]
            A[i, :] -= mult * A[m, :]
            v[i] -= mult * v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N - 1, -1, -1):
        x[m] = v[m]
        for i in range(m + 1, N):
            x[m] -= A[m, i] * x[i]
    return x


def PartialPivot(A_in, v_in):
    """Implement Gaussian Elimination with PartialPivot. This should be
    non-destructive for inputarrays, so we will copy A and v to temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):
        # partial pivoting at each round of m
        mcolumn = A[m:,
                  m]  # define all the mth column values for rows below and including m.
        maxrow = np.argmax(
            abs(mcolumn)) + m  # find the rows with max diagonal value

        if maxrow > m:  # if there are lower rows with values larger than the current row
            A[m, :], A[maxrow, :] = copy(A[maxrow, :]), copy(
                A[m, :])  # flip rows
            v[m], v[maxrow] = copy(v[maxrow]), copy(v[m])  # flip vectors
        # print(A)

        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m + 1, N):
            mult = A[i, m]
            A[i, :] -= mult * A[m, :]
            v[i] -= mult * v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N - 1, -1, -1):
        x[m] = v[m]
        for i in range(m + 1, N):
            x[m] -= A[m, i] * x[i]
    return x

# #test
# a = [1,4,5]
# a.index(max(a))
# b = np.array([[1,2,2],[2,4,4],[5,6,8]])
# g = b[1:,1]
# #find rows with max value
# maxr = np.amax(g)
# maxrow = np.where(b[:,1] == maxr)[0][0]
# b
# m=0
# b[m,:], b[maxrow,:] = copy(b[maxrow,:]),copy(b[m,:])
# a
# a[m], a[maxrow] = a[maxrow],a[m]
# a
# #test functions
# A = np.array([[2,1,4,1],[3,4,-1,-1],[1,-4,1,5],[2,-2,1,3]],float)
# v = np.array([-4,3,9,7],float)
# GaussElim(A,v)
# PartialPivot(A,v)

# #more tests
# b =  np.array([[2,1,4],[3,4,-1],[1,-4,1]],float)
# a =  np.array([-4,3,9],float)
# GaussElim(b,a)
# PartialPivot(b,a)

# #tests
# c = np.array([[0,1,4],[3,4,-1],[1,-4,1]],float)
# a =  np.array([-4,3,9],float)
# PartialPivot(c,a)
# np.linalg.solve(c,a)
