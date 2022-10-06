##显性效应怎么标准化？

import numpy as np
#import math
import multiprocessing
import _thread
import numba as nb 

''' @nb.njit(parallel=True)
def numba_(a, shf):
    b = np.empty_like(a)
    rows_num = a.shape[0]
    cols_num = a.shape[1]
    for i in nb.prange(rows_num):
        b[i, shf:] = a[i, :cols_num - shf]
        b[i, :shf] = a[i, cols_num - shf:]
    return b '''

class Sum: #again, this class is from ParallelPython's example code (I modified for an array and added comments)
    def __init__(self,geno):
        self.value = np.zeros((geno.shape[0],geno.shape[0])) #this is the initialization of the sum
        self.lock = _thread.allocate_lock()
        self.count = 0

    def add(self,value):
        self.count += 1
        self.lock.acquire() #lock so sum is correct if two processes return at same time
        self.value += value #the actual summation
        self.lock.release()
''' def makematrix(geno,start,end,method):
    if method == 'dd':
        geno[geno==2] = 0  ### W矩阵
        zz11 = np.zeros((geno.shape[0],geno.shape[0]))
        for i in range(start,end):
            zz00 = geno * numba_(geno, i) 
            zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
            print(i)
    if method =='aa':
        zz11 = np.zeros((geno.shape[0],geno.shape[0]))
        for i in range(start,end):
            zz00 = geno * numba_(geno, i) 
            zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
            print(i)
    if method =='ad':
        Z = geno  ###Z矩阵
        geno[geno==2] = 0  ### W矩阵
        zz11 = np.zeros((geno.shape[0],geno.shape[0]))
        for i in range(start,end):
            Z[:,-i:] = geno[:,-i:] 
            W = numba_(geno, i)
            W[:,-i:] = Z[:,:i]
            zz00 = Z * W
            zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
            print(i)
    if method =='da':
        Z = geno  ###Z矩阵
        geno[geno==2] = 0  ### W矩阵
        zz11 = np.zeros((geno.shape[0],geno.shape[0]))
        for i in range(start,end):
            geno[:,-i:] = Z[:,-i:] 
            Z_t = numba_(Z, i)
            Z_t[:,-i:] = geno[:,:i]
            zz00 = geno * Z_t
            zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
            print(i)
    return zz11 '''
import numba as nb 
@nb.njit(parallel=True)
def one_matrix(geno,start,end,step):
    zz11 = np.zeros((geno.shape[0],geno.shape[0]))
    for i in range(start,end,step):
        zz00 = geno[:,:-i] * geno[:,i:]
        zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
        print(i)
    return zz11

@nb.njit(parallel=True)
def ad_matrix(geno,Z,start,end,step):
    zz11 = np.zeros((geno.shape[0],geno.shape[0]))
    for i in range(start,end,step):
        zz00 = Z[:,:-i] * geno[:,i:]
        zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
        print(i)
    return zz11
    
@nb.njit(parallel=True)
def da_matrix(geno,Z,start,end,step):

    zz11 = np.zeros((geno.shape[0],geno.shape[0]))
    for i in range(start,end,step):
        zz00 =geno[:,:-i] * Z[:,i:]
        zz11 += (zz00.dot(zz00.T)).astype(np.float64)   ##转置相乘float32就够，但是下一步加在一起，可能几十万或者几百万会影响精度
        print(i)
    return zz11

def makematrix(geno,start,end,method,step):
    if method == 'dd':
        geno[geno==2] = 0  ### W矩阵
        zz11 = one_matrix(geno,start,end,step)
    if method =='aa':
        zz11 = one_matrix(geno,start,end,step)
    if method =='ad':
        Z = geno  ###Z矩阵
        geno[geno==2] = 0  ### W矩阵
        zz11 = ad_matrix(geno,Z,start,end,step)
    if method =='da':
        Z = geno  ###Z矩阵
        geno[geno==2] = 0  ### W矩阵
        zz11 = da_matrix(geno,Z,start,end,step)
    return zz11


def makeG_inter(geno,method,multi=1):
    '''Calculate the epistatic kinship matrix using genotype.

    geno: a numpy matrix of genotype which code as 0,1 and 2. Rows are individuals and columns are SNPs.
    method: str value of 'dd','aa','ad' or 'da'.
    multi: int value. Default value is 1. This function uses multi-threaded calculation by default, 
        if your computer has more than one cpu, you can set the value of multi equal to the number of cpu.
    '''
    method_list =['dd','aa','ad','da']
    if method not in method_list:
        print("ERROR: Parameter method should be in %s" %method_list)
        return
    ###输入的geno,行是个体，列是snp，编码为012
    if not isinstance(geno, np.ndarray):
        print("ERROR: geno data should be numpy ndarray type")
        return
    if np.isnan(geno).any():
        print("ERROR: Nan in geno")
        return
    if not isinstance(multi, int):
        print("ERROR: Parameter multi should be int type!")
        return
    if multi<1:
        print('Error: multi must more than 1.')
        return
    geno = geno.astype(np.float32)
    pool = multiprocessing.Pool(processes=multi)
    sumArr = Sum(geno) #create an instance of callback class and zero the sum
    #num = math.ceil(geno.shape[1]/2)-1
    #index = np.linspace(1,num+1,multi+1,dtype=int)
    for x in range(multi):
        start= x+1 
        end = geno.shape[1]

        singlepoolresult = pool.apply_async(makematrix,(geno,start,end,method,multi),callback=sumArr.add)

    pool.close()
    pool.join() #waits for all the processes to finish
    '''###加一个检测，如果是偶数，最后多加一次平分的
    ###还要在除以对角线均值
    if (geno.shape[1] % 2) == 0: ##偶数
        if method == 'aa':
            geno1,geno2 = np.hsplit(geno,2)
            zz00 = geno1 *  geno2
            zz11 = (zz00.dot(zz00.T)).astype(np.float64) 
            sumArr.value += zz11
        if method == 'dd':
            geno[geno==2] = 0
            geno1,geno2 = np.hsplit(geno,2)
            zz00 = geno1 *  geno2
            zz11 = (zz00.dot(zz00.T)).astype(np.float64) 
            sumArr.value += zz11
        if method == 'ad':
            Z = geno  ###Z矩阵
            geno[geno==2] = 0 
            mid = int(geno.shape[1]/2)
            zz00 = Z[:,:mid] * geno[:,mid:]
            zz11 = (zz00.dot(zz00.T)).astype(np.float64) 
            sumArr.value += zz11
        if method == 'da':
            Z = geno  ###Z矩阵
            geno[geno==2] = 0 
            mid = int(geno.shape[1]/2)
            zz00 = geno[:,:mid] * Z[:,mid:]
            zz11 = (zz00.dot(zz00.T)).astype(np.float64) 
            sumArr.value += zz11'''
    M = sumArr.value/sumArr.value.diagonal().mean()
    return M