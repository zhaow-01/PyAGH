import pandas as pd
import ctypes 
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve_triangular
import PyAGH.FCOEFF
def makeA(data_ord):
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
    PyAGH.FCOEFF.fcoeff(dam_c,
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
    X = Tinv.dot(di).tocsr() ##upper

    Ip=np.identity(N)

    tu = spsolve_triangular(X,Ip,False)
    A = (tu.T).dot(tu)
    return [A, data_ord.iloc[:,0]]