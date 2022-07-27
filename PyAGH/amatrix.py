import pandas as pd
import ctypes 
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve_triangular
import gc
import FUNC

from multiprocessing import Pool
def makeA(data_ord:pd.DataFrame,Sparse=False):
    '''Calculate the additive kinship matrix using pedigree.

    data_ord: a dataframe of pedigree after sort.
    Sparse: bool value. Default value is False. Using sparse matrix can save memory but consume more time.
    '''
    if not isinstance(data_ord, pd.DataFrame): ###必须是data.frame
        print("Please provide data with dataframe type!")
        return
    if data_ord.shape[0] == 0:
        print("data is null!")
        return
    if data_ord.shape[1] != 3:
        print("Data must have three columns, please verify")
        return
    if data_ord.iloc[:,0].shape[0] != data_ord.iloc[:,1].shape[0]:
        print("ID has to be of the same length than sire and dam")
        return
    if data_ord.iloc[:,2].shape[0] != data_ord.iloc[:,1].shape[0]:
        print("sire and dam have to be of the same length")
        return
    ##如果存在na，使用sort_ped函数
    if any(data_ord.isnull().any()):
        print("There is Nan in data, first use the 'sort_ped' function")
        return
    data_ord = data_ord.astype(str)
    ###如果存在重复的行，使用sort——ped函数
    if any(data_ord.duplicated()):
        print("Duplicate in data, first use the 'sort_ped' function")
        return

    if any(data_ord.iloc[:,0].duplicated()): ##检查存在相同个体但父母却不相同的错误,即仅在个体列存在重复的情况
        print("some individuals appear more than once in the pedigree")
        return
    ###检查是否所有父母都为0
    if all(data_ord.iloc[:,1] =="0") and all(data_ord.iloc[:,2]=="0"):
        print("All dams and sires are missing")
        return
    ##检查个体同时出现在父亲列和母亲列的错误
    if any(data_ord[data_ord.iloc[:,1] != "0"].iloc[:,1].isin(data_ord[data_ord.iloc[:,1] != "0"].iloc[:,2])):
        print("Dams appearing as Sires")
        return
    ###检查一个个体本身是不是自己的父母
    if any(data_ord.iloc[:,0] == data_ord.iloc[:,1]) or any(data_ord.iloc[:,0] == data_ord.iloc[:,2]):
        print("Individual appearing as its own Sire or Dam")
        return
    N = len(data_ord.iloc[:,0])
    nPed = pd.DataFrame( {"id":pd.Series(range(N)),
                "sire":pd.Categorical(data_ord.iloc[:,1],categories=data_ord.iloc[:,0]).codes,
                "dam": pd.Categorical(data_ord.iloc[:,2],categories=data_ord.iloc[:,0]).codes
    })
    if any(nPed["id"] < nPed["sire"]) or any(nPed["id"] < nPed["dam"]):
        print("Offspring appearing before their Sires or Dams: first use the 'sort_ped' function")
        return
    nPed[nPed == -1] = N 
    f = ([0.0]*N)+[-1]
    dam = nPed["dam"]
    sire = nPed["sire"]
    dii = [0.0]*N

    dam_c =(ctypes.c_int * N) (*dam)
    sire_c = (ctypes.c_int * N) (*sire)
    f_c = (ctypes.c_double * (N+1)) (*f)
    dii_c = (ctypes.c_double * N) (*dii)
    FUNC.fcoeff(dam_c,
                    sire_c,  
                    f_c,       
                    dii_c, 
                    ctypes.c_int(N),
                    ctypes.c_int(1)
    )
    loc = pd.DataFrame({
            "x": pd.concat([nPed["id"],nPed["id"],pd.Series(range(N))]),
            "y": pd.concat([nPed["sire"],nPed["dam"],pd.Series(range(N))]),
            "data": np.append([-0.5] * (len(nPed["id"])*2),[1.0] *N)
            })
    loc=loc[loc["y"] != N]
    Tinv = sparse.coo_matrix((loc["data"], (loc["y"], loc["x"])), shape=(N, N), dtype=np.float64).tocsr()

    di =  sparse.coo_matrix((np.sqrt(1/np.array(dii_c)), (np.array(range(N)), np.array(range(N)))), shape=(N,N), dtype=np.float64).tocsr()
    #X = Tinv.dot(di).tocsc() ##upper
    X = sparse.csr_matrix(Tinv.dot(di),dtype=np.float32)

    Ip = sparse.identity(N,dtype='float32', format='dia').toarray() ## 最快的方法 7.7秒3W int8

    tu = spsolve_triangular(X,Ip,False,overwrite_A=True, overwrite_b=True)
    if Sparse:
        tu = sparse.csr_matrix(tu,dtype=np.float32)
    del Ip, X
    gc.collect()
    #tu = spsolve_triangular(X,Ip,False)
    #tu = sparse.linalg.inv(X).tocsr()
    A = tu.T.dot(tu)
    return [A, data_ord.iloc[:,0]]

def makeD(data_ord,multi=1):
    '''Calculate the additive kinship matrix using pedigree.

    data_ord: a dataframe of pedigree after sort.
    multi: int value. Default value is 1. This function uses multi-threaded calculation by default, 
    if your computer has more than one cpu, you can set the value of multi equal to the number of cpu.
    '''
    if not isinstance(multi, int):
        print("ERROR: Parameter multi should be int type!")
        return
    if multi<1:
        print('Error: multi must more than 1.')
        return
    A = makeA(data_ord,Sparse=True)
    N = A[0].shape[0]
    c =sparse.triu(A[0],format='csc')
    dA = c.diagonal()
    iAP = c.indices
    pAP = c.indptr
    xAP = c.data
    dij =  [0.0]*len(xAP)
    Di = [0]*len(iAP)
    Dp = [0]*N
    nPed = pd.DataFrame( {"id":pd.Series(range(N)),
                "sire":pd.Categorical(data_ord.iloc[:,1],categories=data_ord.iloc[:,0]).codes,
                "dam": pd.Categorical(data_ord.iloc[:,2],categories=data_ord.iloc[:,0]).codes
    })
    nPed[nPed == -1] = N 
    dam = nPed["dam"]
    sire = nPed["sire"]
    dam_c =(ctypes.c_int * N) (*dam)
    sire_c = (ctypes.c_int * N) (*sire)
    iAP_c = (ctypes.c_int * len(iAP)) (*iAP)
    pAP_c = (ctypes.c_int * len(pAP)) (*pAP)
    xAP_c = (ctypes.c_double * len(xAP)) (*xAP/2)
    
    if multi ==1:
        dij_c = (ctypes.c_double * len(xAP)) (*dij)
        Di_c = (ctypes.c_int * len(iAP)) (*Di)
        Dp_c = (ctypes.c_int * N) (*Dp)
        cnt_c = ctypes.c_int(0)
        FUNC.cald(dam_c,
            sire_c,  
            iAP_c,     
            pAP_c,	         
            xAP_c,
            ctypes.c_int(N),
            dij_c,
            Di_c,
            Dp_c,
            cnt_c
                )
        D = sparse.csr_matrix((dij_c[:cnt_c.value], Di_c[:cnt_c.value], np.append(Dp_c,cnt_c)), shape=(17, 17)).toarray()
        row, col = np.diag_indices_from(D) 
        D[row, col] = np.array(2-dA)
        return [D, data_ord.iloc[:,0]]
    else:
        global pyagh_multi_d_function
        def pyagh_multi_d_function(sub_lA):  ###listA的一部分
            lA_r = sub_lA.shape[0]
            indk_c = (ctypes.c_int * lA_r) (*sub_lA.iloc[:,0])
            indj_c = (ctypes.c_int * lA_r) (*sub_lA.iloc[:,1])
            dij =  [0.0]*lA_r
            dij_c = (ctypes.c_double * lA_r) (*dij)
            FUNC.dijp(
                dam_c,
                sire_c,
                ctypes.c_int(lA_r),
                indk_c,
                indj_c,
                iAP_c,
                pAP_c,
                xAP_c,
                dij_c
            )
            return list(dij_c),list(indk_c),list(indj_c)
        listA = pd.DataFrame({
                'Row':np.repeat(np.arange(0,len(pAP)-1),np.diff(pAP)),
                'Column':iAP
            })
        data = np.array_split(listA, multi)
        dij_list, indk_list,indj_list = [], [],[]
        with Pool() as p:
            for dij,indk,indj in p.map(pyagh_multi_d_function,data):   # or   imap()
                dij_list.extend(dij)
                indk_list.extend(indk)
                indj_list.extend(indj)
        D = sparse.coo_matrix((dij_list,(indk_list,indj_list)),shape=(N,N)).toarray()
        row, col = np.diag_indices_from(D) 
        D[row, col] = np.array(2-dA)
        return [D, data_ord.iloc[:,0]]