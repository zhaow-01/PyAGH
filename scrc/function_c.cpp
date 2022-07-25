#include<iostream>
#include<string>
#include<pybind11/pybind11.h>
#include <algorithm>
namespace py = pybind11;
using namespace std;

extern "C"{  

void fcoeff(
	int *dam,
	int *sire,
	double *f,
	double *dii,
	int *n,
	int *fmiss
){
 
  int     j, k, h, cnt, sj, dj;
  double  ai;
  double  *AN = new double[2*n[0]];
  double  *li = new double[n[0]];

  for(k=0; k<n[0]; k++){
     li[k]=0.0;               // set l to zero
  }
  for(k=0; k<n[0]; k++){
     AN[k]=-1;               // set AN to zero
  }

  for(k=0; k<n[0]; k++){  // iterate through each row of l 
    dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);
    
    if(fmiss[0] == 1){  // only do below if f coefficients NOT supplied by user
      if((k > 0) && (dam[k] == dam[k-1]) && (sire[k] == sire[k-1])){
        f[k] += f[k-1];
      } 
      else {
        li[k] = 1.0;                   // set l_ii to one
        ai=0.0;                        // set a_ii to zero
        j=k;
        cnt=0;
        while(j>=0){
          sj=sire[j];
          dj=dam[j];

          if(sj!= n[0]){
            AN[cnt] = sj;
            li[sj] += 0.5*li[j];
            cnt++;
          }

          if(dj!= n[0]){
            AN[cnt] = dj;
            li[dj] += 0.5*li[j];
            cnt++;
          }

          ai += li[j]*li[j]*dii[j];
          j=-1;

          for(h=0; h<cnt; h++){   // find eldest individual
            if(AN[h]>j){
              j = AN[h];
            }
          }
          for(h=0; h<cnt; h++){   // delete duplicates
            if(AN[h]==j){
              AN[h] -= n[0];
            }
          }
        }  // end of while
        f[k] = ai-1.0;
        for(h=0; h<=k; h++){
          li[h]  = 0.0;            // reset l to zero except l_ii =1
        }

      } // end else for checking if k has same parents as k-1
    }  // end if f missing
  } // end of for
  delete[] AN;
  delete[] li;
}
}  

extern "C"{  

void cald(
	int *dam,
	int *sire,
  int *iAP,     
	int *pAP,	         
	double *xAP,
	int *nAP,
	double *dij,
	int *Di,
  int *Dp,
	int *cnt
){         

  int     lb, step, it, k, j, m, kDam, kSire, jDam, jSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;

  for(k = 0; k < nAP[0]; k++){  // iterate through each column of "A"
    Dp[k] = cnt[0];
    kDam = dam[k];
    kSire = sire[k];
    if((kDam != -999) && (kSire != -999)){
      for(j = pAP[k]; j < pAP[k+1]; j++){ //iterate through all rows of k column
         if(k != iAP[j]){ 
           jDam = dam[iAP[j]];
           jSire = sire[iAP[j]];
           if((jDam != -999) && (jSire != -999)){
             rmfp = 0.0;
             rmmp = 0.0;
             rfmp = 0.0;
             rffp = 0.0;

             m = pAP[max(kDam, jSire)];
             lb = pAP[max(kDam, jSire)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kDam, jSire)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


             m = pAP[max(kDam, jDam)];
             lb = pAP[max(kDam, jDam)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kDam, jDam)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


             m = pAP[max(kSire, jDam)];
             lb = pAP[max(kSire, jDam)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kSire, jDam)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


             m = pAP[max(kSire, jSire)];
             lb = pAP[max(kSire, jSire)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kSire, jSire)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


             dij_tmp = (rmmp*rffp) + (rmfp*rfmp);
             if(dij_tmp != 0.0){
                dij[cnt[0]] = dij_tmp;   
                Di[cnt[0]] = iAP[j];
                cnt[0]++;
             }
           }   
         }

      }   
    }   
  }   
         

            
}
}
extern "C"{  

void dijp(
	int *dam,
	int *sire,
	int *lAr,
	int *indk,
    int *indj,
	int *iAP,     
	int *pAP,	         
	double *xAP,
	double *dij
){         

  int     lb, step, it, k, m, kDam, kSire, jDam, jSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;

  for(k = 0; k < lAr[0]; k++){ 
    kDam = dam[indk[k]];
    kSire = sire[indk[k]];
    if((kDam != -999) && (kSire != -999)){
       if(indk[k] != indj[k]){ 
         jDam = dam[indj[k]];
         jSire = sire[indj[k]];
         if((jDam != -999) && (jSire != -999)){
           rmfp = 0.0;
           rmmp = 0.0;
           rfmp = 0.0;
           rffp = 0.0;


           m = pAP[max(kDam, jSire)];
           lb = pAP[max(kDam, jSire)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kDam, jSire)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


           m = pAP[max(kDam, jDam)];
           lb = pAP[max(kDam, jDam)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kDam, jDam)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


           m = pAP[max(kSire, jDam)];
           lb = pAP[max(kSire, jDam)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kSire, jDam)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


           m = pAP[max(kSire, jSire)];
           lb = pAP[max(kSire, jSire)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kSire, jSire)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


           dij_tmp = (rmmp*rffp) + (rmfp*rfmp);
           if(dij_tmp != 0.0){
              dij[k] = dij_tmp;   
           }   
         }
       }   
    }   
  }   
                    
}
}


PYBIND11_MODULE(FUNC, m) {
  m.def("fcoeff",[](py::buffer dam, py::buffer sire, py::buffer f, py::buffer dii, py::buffer n, py::buffer fmiss) {\
  py::buffer_info dam_info = dam.request();
  py::buffer_info sire_info = sire.request();
  py::buffer_info f_info = f.request();
  py::buffer_info dii_info = dii.request();
  py::buffer_info n_info = n.request();
  py::buffer_info fmiss_info = fmiss.request();
  fcoeff(static_cast<int*>(dam_info.ptr), static_cast<int*>(sire_info.ptr),static_cast<double*>(f_info.ptr),static_cast<double*>(dii_info.ptr),static_cast<int*>(n_info.ptr),static_cast<int*>(fmiss_info.ptr));
  });
  m.def("cald",[](py::buffer dam, py::buffer sire, py::buffer iAP, py::buffer pAP, py::buffer xAP, py::buffer nAP, py::buffer dij, py::buffer Di, py::buffer Dp, py::buffer cnt) {
  py::buffer_info dam_info = dam.request();
  py::buffer_info sire_info = sire.request();
  py::buffer_info iAP_info = iAP.request();
  py::buffer_info pAP_info = pAP.request();
  py::buffer_info xAP_info = xAP.request();
  py::buffer_info nAP_info = nAP.request();
  py::buffer_info dij_info = dij.request();
  py::buffer_info Di_info = Di.request();
  py::buffer_info Dp_info = Dp.request();
  py::buffer_info cnt_info = cnt.request();
  cald(static_cast<int*>(dam_info.ptr), static_cast<int*>(sire_info.ptr),static_cast<int*>(iAP_info.ptr), static_cast<int*>(pAP_info.ptr),static_cast<double*>(xAP_info.ptr),static_cast<int*>(nAP_info.ptr),static_cast<double*>(dij_info.ptr),static_cast<int*>(Di_info.ptr),static_cast<int*>(Dp_info.ptr),static_cast<int*>(cnt_info.ptr));
  });
  m.def("dijp",[](py::buffer dam, py::buffer sire, py::buffer lAr, py::buffer indk, py::buffer indj, py::buffer iAP, py::buffer pAP, py::buffer xAP, py::buffer dij) {
  py::buffer_info dam_info = dam.request();
  py::buffer_info sire_info = sire.request();
  py::buffer_info lAr_info = lAr.request();
  py::buffer_info indk_info = indk.request();
  py::buffer_info indj_info = indj.request();
  py::buffer_info iAP_info = iAP.request();
  py::buffer_info pAP_info = pAP.request();
  py::buffer_info xAP_info = xAP.request();
  py::buffer_info dij_info = dij.request();

  dijp(static_cast<int*>(dam_info.ptr), static_cast<int*>(sire_info.ptr),static_cast<int*>(lAr_info.ptr), static_cast<int*>(indk_info.ptr),static_cast<int*>(indj_info.ptr),static_cast<int*>(iAP_info.ptr),static_cast<int*>(pAP_info.ptr),static_cast<double*>(xAP_info.ptr),static_cast<double*>(dij_info.ptr));
  });
}


