import numpy as np
import pandas as pd
import polars as pl
import gc

## genoFile_list:两个群体合并在一起的geno文件，格式必须是traw,使用一个txt文件整合文件名
## method=0 sum{[(xij - 2pi)*(xik - 2pi)] / [2pi(1-pi)]}/N as described in Yang et al. 2010 Nat Genet. 
## method=1 sum[(xij - 2pi)(xik - 2pi)] / sum[2pi(1-pi)].
## method=2 G_new if method ==2,分别提供两个群体的数量，n1: 1号群里个体数量 n2:2号群体数量
## method=3 G_chen

def makeG(genoFile_list,method,n1=0,n2=0): 
    method_list =[0,1,2,3]
    if not isinstance(method, int):
        print("ERROR: Parameter method should be a int")
        return
    if method not in method_list:
        print("ERROR: Parameter method should be in %s" %method_list)
        return
    try:
        genofile = pd.read_table(genoFile_list,header=None)  ##读取包含每条染色体文件名的文件，失败报错
        for i in range(genofile.shape[0]):
            if genofile.iloc[i,0].endswith('.traw'):  ##判断是不是traw格式,并不知道文件存不存在
                pass
            else:
                print("ERROR: genotype file type should be .traw of PLINK")
                return
    except Exception as reason:
        print(str(reason))
    try:
        temp = open(genofile.iloc[0,0],"r")  
        line1 = temp.readline().split() 
        ###提取出个体id，最后和G矩阵一起返回                  
        geno_id = line1[6:]                    

        N = len(line1)-6 ##个体数 ##只读第一行标题 N在方法2和3中有需要
        
        G = np.zeros(N)  ###G在方法0 1 中需要

    except Exception as reason:
        print(str(reason))
    if method == 2:
        ##判断一下n1+n2=N
        if not isinstance(n1, int) or not isinstance(n2, int):
            print("ERROR: Parameter n1 and n2 should be a int")
            return
        if  n1 <= 0  or n2 <=0:
            print("ERROR: Parameter n1 and n2 should be much than 0")
            return 
        
        if n1+n2 != N:
            print("ERROR: Parameter n1 + n2 is not equal to the number of indivadual. ")
            return 

        G11=np.zeros((n1,n1))
        G12=np.zeros((n1,n2))
        G22=np.zeros((n2,n2))
        sum2pq_xx=0
        sum2pq_tb=0
        for i in range(genofile.shape[0]):
            chr_file = genofile.iloc[i,0]
            try:
                all_gen = pl.read_csv(chr_file,sep="\t",null_values="NA")
            except Exception as reason:
                print(reason)
                break
            ###这里的Z矩阵是和文献相比转置的，行是snp，列是id
            xx_gen = all_gen[:,6:(n1+6)]
            tb_gen = all_gen[:,(n1+6):]
            ###检查是否有空值，再加一个检查是不是某一个snp全是空值，全是空值的话，没法补
            if xx_gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in first population, and imputing markers with mean value.")
                xx_gen = xx_gen.fill_null(xx_gen.mean(axis=1))
                if xx_gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals in first population, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    xx_gen = xx_gen.fill_null(0)   
            if tb_gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in second population, and imputing markers with mean value.")
                tb_gen = tb_gen.fill_null(tb_gen.mean(axis=1))
                if tb_gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals in second population, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    tb_gen = tb_gen.fill_null(0)   

            p_xx = xx_gen.mean(axis=1)/2
            q_xx = 1- p_xx
            temp_2pq_xx = (2 * p_xx * q_xx).sum() #每一个snp的2pq

            p_tb = tb_gen.mean(axis=1)/2
            q_tb = 1- p_tb
            temp_2pq_tb = (2 * p_tb * q_tb).sum()
        
            z_xx = (xx_gen - 2*p_xx).to_numpy()
            z_tb = (tb_gen - 2*p_tb).to_numpy()
            G22_temp = z_tb.T.dot(z_tb)
            G12_temp = z_xx.T.dot(z_tb)
            G11_temp = z_xx.T.dot(z_xx)
            sum2pq_xx = sum2pq_xx + temp_2pq_xx
            sum2pq_tb = sum2pq_tb + temp_2pq_tb
            
            G11 = G11 + G11_temp
            G12 = G12 + G12_temp
            G22 = G22 + G22_temp
            #print(i)
            del all_gen, xx_gen, tb_gen, z_tb, z_xx
            gc.collect()
        G11_new = G11/sum2pq_xx
        G12_new = G12/(sum2pq_tb ** 0.5 * sum2pq_xx ** 0.5)
        G21_new = G12_new.T
        G22_new = G22/sum2pq_tb
        gmat = np.vstack((np.hstack((G11_new,G12_new)),np.hstack((G21_new,G22_new))))
        #np.savetxt("G_new_wgs.txt",gmat,fmt="%.10f")
        return [gmat,geno_id]
    if method == 0:
        M = 0 
        for i in range(genofile.shape[0]):
            chr_file = genofile.iloc[i,0]
            try:
                all_gen = pl.read_csv(chr_file,sep="\t",null_values="NA")
            except Exception as reason:
                print(reason)
                break
            gen = all_gen[:,6:]
            if gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in data, and imputing markers with mean value.")
                gen = gen.fill_null(gen.mean(axis=1))
                if gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    gen = gen.fill_null(0)  
            p = gen.mean(axis=1)/2
            q = 1-p
            dii =1/ (2 * p * q)
            Z = (gen - 2 * p).to_numpy()
            ZD = Z.T * dii.to_numpy()
            G_temp = ZD.dot(Z)
            G = G + G_temp
            M = M + gen.shape[0]
            del all_gen,gen, dii, ZD, Z, G_temp
            gc.collect()
        G = G/M
        return [G,geno_id]
    if method == 1:
        sum2pq = 0
        for i in range(genofile.shape[0]):
            chr_file = genofile.iloc[i,0]
            try:
                all_gen = pl.read_csv(chr_file,sep="\t",null_values="NA")
            except Exception as reason:
                print(reason)
                break
            gen = all_gen[:,6:]
            if gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in data, and imputing markers with mean value.")
                gen = gen.fill_null(gen.mean(axis=1))
                if gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    gen = gen.fill_null(0)  
            p = gen.mean(axis=1)/2
            q = 1-p
            temp_2pq = (2 * p * q).sum()
            Z = (gen - 2*p).to_numpy()
            G_temp = Z.T.dot(Z)
            G = G + G_temp
            sum2pq = sum2pq + temp_2pq
            del all_gen, gen, Z, G_temp
            gc.collect()
        G = G/sum2pq
        return [G,geno_id]
    if method == 3:
        ##判断n1 和n2,  n1+n2应该等于N
        if not isinstance(n1, int) or not isinstance(n2, int):
            print("ERROR: Parameter n1 and n2 should be a int")
            return
        if  n1 <= 0  or n2 <=0:
            print("ERROR: Parameter n1 and n2 should be much than 0")
            return 
        if n1+n2 != N:
            print("ERROR: Parameter n1 + n2 is not equal to the number of indivadual. ")
            return 
        G11=np.zeros((n1,n1))
        G12=np.zeros((n1,n2))
        G22=np.zeros((n2,n2))
        sum2pq_xx=0
        sum2pq_tb=0
        xxbytb = 0
        for i in range(genofile.shape[0]):
            chr_file = genofile.iloc[i,0]
            try:
                all_gen = pl.read_csv(chr_file,sep="\t",null_values="NA")
            except Exception as reason:
                print(reason)
                break
            ###这里的Z矩阵是和文献相比转置的，行是snp，列是id
            xx_gen = all_gen[:,6:(n1+6)]
            tb_gen = all_gen[:,(n1+6):]
            if xx_gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in first population, and imputing markers with mean value.")
                xx_gen = xx_gen.fill_null(xx_gen.mean(axis=1))
                if xx_gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals in first population, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    xx_gen = xx_gen.fill_null(0)   
            if tb_gen.null_count().sum(axis=1)[0] > 0:
                print("Warning: there are mising values in second population, and imputing markers with mean value.")
                tb_gen = tb_gen.fill_null(tb_gen.mean(axis=1))
                if tb_gen.null_count().sum(axis=1)[0] > 0:
                    print("Warning: Some SNPS are missing in all individuals in second population, and imputing markers with 0.")
                    print("Warning: It is recommended that impute or filter the data first.")
                    tb_gen = tb_gen.fill_null(0)
            p_xx = xx_gen.mean(axis=1)/2
            q_xx = 1- p_xx
            xx_gen_list = 2 * p_xx * q_xx #每一个snp的2pq

            p_tb = tb_gen.mean(axis=1)/2
            q_tb = 1- p_tb
            tb_gen_list = 2 * p_tb * q_tb #每一个snp的2pq

            temp_2pq_xx = xx_gen_list.sum()  ##2pq的和
            temp_2pq_tb = tb_gen_list.sum()
            temp_xxbytb = ((xx_gen_list * tb_gen_list) ** 0.5).sum()  ##
            
            z_xx = (xx_gen - 2 * p_xx).to_numpy()
            z_tb = (tb_gen - 2 * p_tb).to_numpy()
            
            G11_temp = z_xx.T.dot(z_xx)
            G12_temp = z_xx.T.dot(z_tb)
            G22_temp = z_tb.T.dot(z_tb)
            
            sum2pq_xx = sum2pq_xx + temp_2pq_xx
            sum2pq_tb = sum2pq_tb + temp_2pq_tb
            xxbytb = xxbytb + temp_xxbytb
            G11 = G11 + G11_temp
            G12 = G12 + G12_temp
            G22 = G22 + G22_temp
            #print(i)
            del all_gen, xx_gen, tb_gen, z_tb, z_xx
            gc.collect()
        G11_chen = G11/sum2pq_xx
        G12_chen = G12/xxbytb
        G21_chen = G12_chen.T
        G22_chen = G22/sum2pq_tb
        gmat = np.vstack((np.hstack((G11_chen,G12_chen)),np.hstack((G21_chen,G22_chen))))
        #np.savetxt("G_new_wgs.txt",gmat,fmt="%.10f")
        return [gmat,geno_id]