import pandas as pd
import numpy as np
from numbers import Number
import sympy

def makeH(G,A,w):
    ##提供的参数对不对 列表里包括两个元素，第一个是np矩阵，第二个是series，两者的长度要一致 w是一个数字
    ##A[0],G[0]不能有Nan
    if not isinstance(G, list):
        print("ERROR: Parameter G should be a list!")
        return
    if len(G) != 2:
        print("ERROR: G should have 2 elements")
        return
    if not isinstance(G[1], pd.Series):
        print("ERROR: G should have 2 elements with numpy ndarray and pandas Series")
        return
    if not isinstance(G[0], np.ndarray):
        print("ERROR: G should have 2 elements with numpy ndarray and pandas Series")
        return
    if np.isnan(G[0]).any():
        print("ERROR: Nan in G")
        return
    if G[0].shape[0] != len(G[1]):
        print("ERROR: The dimension of G is not equal to the number of individual with genotyped")
        return
    ####A
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
        ## W
    if not isinstance(w, Number):
        print("ERROR: Parameter w should be a number")
        return
    
    if w<0 or w>1:
        print("ERROR: Parameter w should between 0 and 1")
        return
    
    ##要检查一下测序个体是不是都在系谱
    A[1] = A[1].astype(str)
    G[1] = G[1].astype(str)
    if not all(G[1].isin(A[1])):
        print("ERROR: not all individuals with genotyped are in A matrix")
        return

    geno_id_loc = A[1][A[1].isin(G[1])] 
    index_geno = geno_id_loc.iloc[list(map(geno_id_loc.tolist().index,G[1]))].index ###按照提供的测序个体的顺序提取的在A矩阵中的index
    index_nogeno = A[1][~A[1].isin(G[1])].index ##剩下的非测序个体的index
    ##提取子集
    A11 = A[0][index_nogeno,:][:,index_nogeno]
    A12 = A[0][index_nogeno,:][:,index_geno]
    A21 = A[0][index_geno,:][:,index_nogeno]
    A22 = A[0][index_geno,:][:,index_geno]
    iA22 = np.linalg.inv(A22)
    ave_dia_G = np.trace(G[0])/G[0].shape[0]
    ave_offdia_G = (np.sum(G[0])-np.trace(G[0]))/(G[0].shape[0]*G[0].shape[0]-G[0].shape[0])
    ave_dia_A22 = np.trace(A22)/A22.shape[0]
    ave_offdia_A22 = (np.sum(A22)-np.trace(A22))/(A22.shape[0]*A22.shape[0]-A22.shape[0])
    a = sympy.Symbol("a")
    b = sympy.Symbol("b")
    ab = sympy.solve([ave_dia_G*b+a-ave_dia_A22,ave_offdia_G*b+a-ave_offdia_A22],[a,b])

    G_star = ab[a]+ab[b]*G[0]

    G_w = (1-w) * G_star + w * A22
    G_w =np.around(G_w.astype(float),6)

    H11 = A11+A12.dot(iA22).dot(G_w-A22).dot(iA22).dot(A21)
    H12 = A12.dot(iA22).dot(G_w)
    H21 = G_w.dot(iA22).dot(A21)
    H22 = G_w
    H = np.round(np.hstack((np.vstack((H11,H21)),np.vstack((H12,H22)))),2)
    return [H,pd.concat([A[1][index_nogeno], A[1][index_geno]]).reset_index(drop=True)]
    