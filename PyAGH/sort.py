import pandas as pd
class ped:
    def sort_ped(self,data):
        if not isinstance(data, pd.DataFrame): ###必须是data.frame
            print("Please provide data with dataframe type!")
            return
        if data.shape[1] != 3:
            print("Data must have three columns, please verify")
            return
        if data.iloc[:,0].shape[0] != data.iloc[:,1].shape[0]:
            print("ID has to be of the same length than sire and dam")
            return
        if data.iloc[:,2].shape[0] != data.iloc[:,1].shape[0]:
            print("sire and dam have to be of the same length")
            return
        ###输入的数据中缺失值为NA和0，如果存在NA，转换为0
        data.fillna(0)
        ###去除完全重复的行 
        data = data[~data.duplicated()]

        data = data.astype(str)

        if any(data.iloc[:,0].duplicated()): ##检查存在相同个体但父母却不相同的错误,即仅在个体列存在重复的情况
            print("some individuals appear more than once in the pedigree")
            return
        ###检查是否所有父母都为0
        if all(data.iloc[:,1] =="0") and all(data.iloc[:,2]=="0"):
            print("All dams and sires are missing")
            return
        ##检查个体同时出现在父亲列和母亲列的错误
        if any(data[data.iloc[:,1] != "0"].iloc[:,1].isin(data[data.iloc[:,1] != "0"].iloc[:,2])):
            print("Dams appearing as Sires")
            return
        ###检查一个个体本身是不是自己的父母
        if any(data.iloc[:,0] == data.iloc[:,1]) or any(data.iloc[:,0] == data.iloc[:,2]):
            print("Individual appearing as its own Sire or Dam")
            return
        
        tmp = pd.concat([data.iloc[:,1],(data.iloc[:,2])]).drop_duplicates().astype(str)
        tmp.drop(tmp[tmp == "0"].index,inplace=True) #缺失值为0
        index = ~tmp.isin(data.iloc[:,0].astype(str))
        missingP =pd.Series(dtype=str)
        if any(index):
            missingP = tmp[index]
        lable01 = pd.concat([missingP,data.iloc[:,0]]).reset_index(drop=True)  
        sire01 = pd.concat([pd.Series(["0"]*len(missingP)),data.iloc[:,1]]).reset_index(drop=True)
        dam01 = pd.concat([pd.Series(["0"]*len(missingP)),data.iloc[:,2]]).reset_index(drop=True)
        sire = pd.Categorical(sire01,categories=lable01).codes
        dam = pd.Categorical(dam01,categories=lable01).codes
        nped = len(lable01)
        lable = pd.Series(range(nped))
        sire = pd.Series(sire)
        dam = pd.Series(dam)
        
        self.__pede = pd.DataFrame({"id":lable,
                            "sire":sire,
                            "dam":dam,
                            "gene_max": None,
                            "gene_min":None
                            })
        noParents=(self.__pede["sire"] == -1) & (self.__pede["dam"] == -1)
        self.__pede.loc[noParents,"gene_max"] = 0
        self.__pede.loc[noParents,"gene_min"] = 0
        self.__pede = self.__pede.to_dict()
        
        self.__flag = True

        for i in range(nped):
            if self.__flag:
                if self.__pede["gene_max"][i]==None:
                    
                    self.__count = 0
                    self.__getGenAncestors(i)
            else:
                break
        if self.__flag:    
            ans = pd.DataFrame({"id":lable01,
                        "sire":sire01,
                        "dam":dam01,
                        "gene_max":pd.Series(self.__pede["gene_max"]),
                        "gene_min":pd.Series(self.__pede["gene_min"])
                        })
            ans.sort_values(["gene_max","gene_min"],inplace=True,ascending=[True,True])
            ans.reset_index(inplace=True,drop=True)
            return ans
        else:
            print("ERROR: infinite pedigree loop involving individual")
            print("The first error individual is %s" % lable01[i-1])
    def __getGenAncestors(self,i):
        
        self.__count += 1
        if self.__count < 100:
            sire = self.__pede["sire"][i]
            dam = self.__pede["dam"][i]
            genP1=0
            genP2=0
            if sire != -1:
                tmpgenP1 = self.__pede["gene_max"][sire]
                if tmpgenP1 == None:
                    self.__getGenAncestors(sire)
                    genP1 = 1+ genP1
                else:
                    genP1 = 1+tmpgenP1
            if dam != -1:
                tmpgenP2 = self.__pede["gene_max"][dam]
                if tmpgenP2 == None:
                    self.__getGenAncestors(dam)
                    genP2 = 1+ genP2
                else:
                    genP2 = 1+tmpgenP2
            self.__pede["gene_max"][i] = max(genP1,genP2)
            self.__pede["gene_min"][i] = min(genP1,genP2)  
        else:
            #print("infinite pedigree loop involving individual")
            self.__flag = False
            return
def sortPed(d):
    a = ped()
    return a.sort_ped(d)
