"""Skeleton file for your solution to the shortest-path kernel."""
import numpy as np
import math


def floyd_warshall(A):
    """Implement the Floyd--Warshall on an adjacency matrix A.

    Parameters
    ----------
    A : `np.array` of shape (n, n)
        Adjacency matrix of an input graph. If A[i, j] is `1`, an edge
        connects nodes `i` and `j`.

    Returns
    -------
    An `np.array` of shape (n, n), corresponding to the shortest-path
    matrix obtained from A.
    """
    n=len(A)
    S=np.matrix((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i,j]==0:
                    A[i,j]=999

    for k in range(n):
        for i in range(n):
            for j in range(n):
                if A[i,j] > A[i,k]+A[k,j]:
                    A[i,j]=A[i,k]+A[k,j]


    return A


def sp_kernel(S1, S2):
    """Calculate shortest-path kernel from two shortest-path matrices.

    Parameters
    ----------
    S1: `np.array` of shape (n, n)
        Shortest-path matrix of the first input graph.

    S2: `np.array` of shape (m, m)
        Shortest-path matrix of the second input graph.

    Returns
    -------
    A single `float`, corresponding to the kernel value of the two
    shortest-path matrices
    """
    n=min(len(S1),len(S2))
    k_sp=0
    for i in range(n):
        for j in range(n-i):
            if S1[i,j]==S2[i,j]:
                k_walk=1
            else:
                k_walk=0
            k_sp=k_sp+k_walk
    return k_sp
