#include<iostream>
#include<string>

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
     li[k]=0.0;             
  }
  for(k=0; k<n[0]; k++){
     AN[k]=-1;            
  }

  for(k=0; k<n[0]; k++){  
    dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);
    
    if(fmiss[0] == 1){ 
      if((k > 0) && (dam[k] == dam[k-1]) && (sire[k] == sire[k-1])){
        f[k] += f[k-1];
      } 
      else {
        li[k] = 1.0;               
        ai=0.0;                     
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

          for(h=0; h<cnt; h++){   
            if(AN[h]>j){
              j = AN[h];
            }
          }
          for(h=0; h<cnt; h++){  
            if(AN[h]==j){
              AN[h] -= n[0];
            }
          }
        }  
        f[k] = ai-1.0;
        for(h=0; h<=k; h++){
          li[h]  = 0.0;           
        }

      } 
    }  
  } 
  delete[] AN;
  delete[] li;
}
}  

