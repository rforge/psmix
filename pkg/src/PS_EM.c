#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>

/*
  Mixture Modeling of Population Structure with EM Algorithm
  1) Admixture Modeling: individual and marker level
  2) Missing value handling: missing completely at random
*/

void ps_admix1(double *PIik, double *Zimak, double *Gkmj, int *K, int *N, int *M, int *Mn, int *MnCS,
	       int *Geno, int *itMax, double *eps)
{
  //   K = number of parent populations
  //   N = number of people
  //   M = number of markers
  //   Mn = number of allelles at each marker
  //   MnCS = c(0, cumsum(Mn))
  //   Geno = Genotype information
  //   itMax = maximum number of iterations
  //   Zimak = missing population indicator probability for each allele
  //   eps = convergence criterion
  
  int i, k, m, a, imak, ik, km, j, it, MnT = MnCS[*M];
  double tmp1, tmp2, delta;
  double *PIik0, *Gkmj0;

  PIik0 = (double *) R_alloc( (*N)*(*K), sizeof(double) );
  Gkmj0 = (double *) R_alloc( (*K)*MnT, sizeof(double) );

  for(it=0; it < *itMax; it++){
    memcpy( PIik0, PIik, (*N)*(*K)*sizeof(double) );
    memcpy( Gkmj0, Gkmj, (*K)*MnT*sizeof(double) );
    
    /* E-step: Estimate Zimak */
    // memset( Zimak, 0, (*N)*(*M)*2*(*K)*sizeof(double) );
    for(i=0; i < *N; i++){
      for(m=0; m < *M; m++){
	for(a=0; a < 2; a++){
	  j = Geno[2*(*N)*m + 2*i+a]-1;
	  if(j >= 0){
	    tmp2 = 0;
	    for(k=0; k < *K; k++){
	      imak = 2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a;
	      Zimak[imak] = PIik[(*N)*k+i]*Gkmj[(*K)*(MnCS[m]+j) + k];
	      tmp2 += Zimak[imak];
	    }
	    for(k=0; k < *K; k++){
	      Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a] /= tmp2;
	    }
	  }
	}
      }
    }

    /* M-step: Estimate PIik */
    memset( PIik, 0, (*N)*(*K)*sizeof(double) );
    for(i=0; i < *N; i++){
      tmp2 = 0;
      for(k=0; k < *K; k++){
	ik = (*N)*k + i;
	for(m=0; m < *M; m++){
	  for(a=0; a < 2; a++){
	    PIik[ik] += Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a];
	  }
	}
	tmp2 += PIik[ik];
      }
      for(k=0; k < *K; k++){
	PIik[(*N)*k+i] /= tmp2;
      }
    }
    
    /* M-step: Estimate Gkmj */
    memset( Gkmj, 0, (*K)*MnT*sizeof(double) );
    for(k=0; k < *K; k++){
      for(m=0; m < *M; m++){
	km = (*K)*MnCS[m] + k;
	tmp2 = 0;
	for(i=0; i < *N; i++){
	  for(a=0; a< 2; a++){
	    j = Geno[2*(*N)*m + 2*i+a]-1;
	    if(j >= 0){
	      tmp1 = Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a];
	      Gkmj[km + (*K)*j] += tmp1;
	      tmp2 += tmp1;
	    }
	  }
	}
	for(j=0; j < Mn[m]; j++){
	  Gkmj[km + (*K)*j] /= tmp2;
	}
      }
    }
    
    /* Check for convergence */
    delta = 0;
    for(j=0; j < (*N)*(*K); j++){
      tmp1 = PIik[j] - PIik0[j];
      delta += tmp1*tmp1;
    }
    for(j=0; j < (*K)*MnT; j++){
      tmp1 = Gkmj[j] - Gkmj0[j];
      delta += tmp1*tmp1;
    }
    delta = sqrt(delta);
    if(delta <= *eps){
      break;
    }

  }
  
  if(it < *itMax){
    Rprintf("Converge after %d steps! Distance = %lf. \n", it, delta);
  } else{
    Rprintf("No convergence after %d steps! Distance = %lf. \n", it, delta);
  }
  
}


void ps_admix2(double *PIimk, double *Zimak, double *Gkmj, int *K, int *N, int *M, int *Mn, int *MnCS,
	       int *Geno, int *itMax, double *eps)
{
  //   K = number of parent populations
  //   N = number of people
  //   M = number of markers
  //   Mn = number of allelles at each marker
  //   MnCS = c(0, cumsum(Mn))
  //   Geno = Genotype information
  //   itMax = maximum number of iterations
  //   Zimak = missing population indicator probability for each allele
  //   eps = convergence criterion
  
  int i, k, m, a, imak, imk, km, j, it, MnT = MnCS[*M];
  double tmp1, tmp2, delta;
  double *PIimk0, *Gkmj0;

  PIimk0 = (double *) R_alloc( (*N)*(*M)*(*K), sizeof(double) );
  Gkmj0 = (double *) R_alloc( (*K)*MnT, sizeof(double) );

  for(it=0; it < *itMax; it++){
    memcpy( PIimk0, PIimk, (*N)*(*M)*(*K)*sizeof(double) );
    memcpy( Gkmj0, Gkmj, (*K)*MnT*sizeof(double) );
    
    /* E-step: Estimate Zimak */
    // memset( Zimak, 0, (*N)*(*M)*2*(*K)*sizeof(double) );
    for(i=0; i < *N; i++){
      for(m=0; m < *M; m++){
	for(a=0; a < 2; a++){
	  j = Geno[2*(*N)*m + 2*i+a]-1;
	  if(j >= 0){
	    tmp2 = 0;
	    for(k=0; k < *K; k++){
	      imak = 2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a;
	      Zimak[imak] = PIimk[(*N)*(*M)*k + (*N)*m + i]*Gkmj[(*K)*(MnCS[m]+j) + k];
	      tmp2 += Zimak[imak];
	    }
	    for(k=0; k < *K; k++){
	      Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a] /= tmp2;
	    }
	  }
	}
      }
    }

    /* M-step: Estimate PIimk */
    memset( PIimk, 0, (*N)*(*M)*(*K)*sizeof(double) );
    for(i=0; i < *N; i++){
      for(m=0; m < *M; m++){
	tmp2 = 0;
	for(k=0; k < *K; k++){
	  imk = (*N)*(*M)*k + (*N)*m + i;
	  for(a=0; a < 2; a++){
	    PIimk[imk] += Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a];
	  }
	  tmp2 += PIimk[imk];
	}
	if(tmp2 > 0){
	  for(k=0; k < *K; k++){
	    PIimk[(*N)*(*M)*k+(*N)*m+i] /= tmp2;
	  }
	}
      }
    }
    
    /* M-step: Estimate Gkmj */
    memset( Gkmj, 0, (*K)*MnT*sizeof(double) );
    for(k=0; k < *K; k++){
      for(m=0; m < *M; m++){
	km = (*K)*MnCS[m] + k;
	tmp2 = 0;
	for(i=0; i < *N; i++){
	  for(a=0; a< 2; a++){
	    j = Geno[2*(*N)*m + 2*i+a]-1;
	    if(j >= 0){
	      tmp1 = Zimak[2*(*N)*(*M)*k + 2*(*N)*m + 2*i+a];
	      Gkmj[km + (*K)*j] += tmp1;
	      tmp2 += tmp1;
	    }
	  }
	}
	if(tmp2 > 0){
	  for(j=0; j < Mn[m]; j++){
	    Gkmj[km + (*K)*j] /= tmp2;
	  }
	}
      }
    }
    
    /* Check for convergence */
    delta = 0;
    for(j=0; j < (*N)*(*M)*(*K); j++){
      tmp1 = PIimk[j] - PIimk0[j];
      delta += tmp1*tmp1;
    }
    for(j=0; j < (*K)*MnT; j++){
      tmp1 = Gkmj[j] - Gkmj0[j];
      delta += tmp1*tmp1;
    }
    delta = sqrt(delta);
    if(delta <= *eps){
      break;
    }

  }
  
  if(it < *itMax){
    Rprintf("Converge after %d steps! Distance = %lf. \n", it, delta);
  } else{
    Rprintf("No convergence after %d steps! Distance = %lf. \n", it, delta);
  }
  
}

