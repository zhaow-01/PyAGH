#include<iostream>
#include<string>
#include<pybind11/pybind11.h>
namespace py = pybind11;
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

PYBIND11_MODULE(FCOEFF, m) {
  m.def("fcoeff",[](py::buffer dam, py::buffer sire, py::buffer f, py::buffer dii, py::buffer n, py::buffer fmiss) {\
  py::buffer_info dam_info = dam.request();
  py::buffer_info sire_info = sire.request();
  py::buffer_info f_info = f.request();
  py::buffer_info dii_info = dii.request();
  py::buffer_info n_info = n.request();
  py::buffer_info fmiss_info = fmiss.request();
  fcoeff(static_cast<int*>(dam_info.ptr), static_cast<int*>(sire_info.ptr),static_cast<double*>(f_info.ptr),static_cast<double*>(dii_info.ptr),static_cast<int*>(n_info.ptr),static_cast<int*>(fmiss_info.ptr));
  });

}


