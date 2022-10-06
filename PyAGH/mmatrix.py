#Ross, E. M., Moate, P. J., Marett, L. C., Cocks, B. G., & Hayes, B. J. (2013). Metagenomic predictions: From microbiome to complex health and environmental phenotypes in humans and cattle. 
import numpy as np
def makeM(otu):
    '''Calculate the kinship matrix using OTU.
    
    otu: OTU data in numpy ndarray type. Rows are individuals and columns are OTUs.
    '''
    ###检测没有空值
    if not isinstance(otu, np.ndarray):
        print("ERROR: OTU matrix should be numpy ndarray")
        return
    if np.isnan(otu).any():
        print("ERROR: Nan in OTU")
        return
    otu_log = np.log(otu+0.001)
    sd_logSj  = np.std(otu_log,axis=0)
    mean_logSj = np.mean(otu_log,axis=0)
    X =(otu_log-mean_logSj)/sd_logSj
    O = (1/otu.shape[1])*X.dot(X.T)
    return O

def makeT(trans):
    '''Calculate the kinship matrix using transcriptome data.
    
    trans: transcriptome data in numpy ndarray type. Rows are individuals and columns are transcriptome.
    '''
    ###检测没有空值
    if not isinstance(trans, np.ndarray):
        print("ERROR: trans matrix should be numpy ndarray")
        return
    if np.isnan(trans).any():
        print("ERROR: Nan in trans")
        return

    sd_trans  = np.std(trans,axis=0)
    mean_trans = np.mean(trans,axis=0)
    X =(trans-mean_trans)/sd_trans
    T = (1/trans.shape[1])*X.dot(X.T)
    return T