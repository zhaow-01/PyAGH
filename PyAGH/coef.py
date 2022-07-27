###返回亲缘系数和近交系数
import numpy as np
import pandas as pd

def coefInbreeding(A):
    '''Calculate the inbreeding coefficients using kinship matrix.

    A: a list with two elements, A[0] is kinship matrix; A[1] is id series
    '''
    if not isinstance(A, list):
        print("ERROR: Parameter A should be a list!")
        return
    if len(A) != 2:
        print("ERROR: A should have 2 elements")
        return
    if not isinstance(A[1], pd.Series):
        print("ERROR: A should have 2 elements with numpy ndarray and pandas Series")
        return
    if not isinstance(A[0], np.ndarray):
        print("ERROR: A should have 2 elements with numpy ndarray and pandas Series")
        return
    if np.isnan(A[0]).any():
        print("ERROR: Nan in A")
        return
    if A[0].shape[0] != len(A[1]):
        print("ERROR: The dimension of A is not equal to the number of individual with genotyped")
        return
 
    coef = pd.DataFrame({'ID':list(A[1]),
                'F':list(np.round(np.diag(A[0])-1,6))
           })
    return coef

def getLowerTriangularIndices(n):
    return [list(x) for x in np.triu_indices(n)]
def coefKinship(A):
    '''Calculate the relationship coefficients using kinship matrix.

    A: a list with two elements, A[0] is kinship matrix; A[1] is id series
    '''
    if not isinstance(A, list):
        print("ERROR: Parameter A should be a list!")
        return
    if len(A) != 2:
        print("ERROR: A should have 2 elements")
        return
    if not isinstance(A[1], pd.Series):
        print("ERROR: A should have 2 elements with numpy ndarray and pandas Series")
        return
    if not isinstance(A[0], np.ndarray):
        print("ERROR: A should have 2 elements with numpy ndarray and pandas Series")
        return
    if np.isnan(A[0]).any():
        print("ERROR: Nan in A")
        return
    if A[0].shape[0] != len(A[1]):
        print("ERROR: The dimension of A is not equal to the number of individual with genotyped")
        return
    n = A[0].shape[0]
    mat = np.zeros((n,n))
    for i in range(0,n):
        for j in range(i,n):
            mat[i,j] = A[0][i,j]/(A[0][i,i]*A[0][j,j])**0.5

    mat = np.round(mat,6)
    name_index = getLowerTriangularIndices(n)
    inbreed = pd.DataFrame({'ID1':list(A[1][np.array(name_index[0])]),
        'ID2':list(A[1][np.array(name_index[1])]),
        'r':list(mat[np.triu_indices(n)]),
        
        })
    return inbreed
