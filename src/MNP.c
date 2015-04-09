/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai and David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cMNPgibbs(int *piNDim, int *piNCov, int *piNSamp, int *piNGen, 
	       double *b0,    /* prior mean for beta */
	       double *pdA0, int *piNu0, double *pdS, double *pdX, 
	       int *y,        /* response variable: -1 for missing */
	       double *pdbeta, double *pdSigma, int *piImp, 
	       int *invcdf,   /* use inverse cdf for TruncNorm? */
	       int *piBurnin, /* the number of burnin */
	       int *piKeep,
	       int *itrace,
	       int *verbose,  /* 1 if extra print is needed */ 
	       int *piMoP,    /* 1 if Multinomial ordered Probit */
	       int *latent,   /* 1 if W is stored */
	       double *pdStore){
  
  /* paramerters from R */
  int n_samp = *piNSamp; /* sample size */
  int n_gen = *piNGen;   /* number of gibbs draws */
  int n_cov = *piNCov;   /* The number of covariates */
  int n_dim = *piNDim;   /* The number of indentifiable dimension (p-1) */
  int nu0 = *piNu0;      /* nu0: degree of freedom, S: scale matrix */
  int imp = *piImp;      /* 0: improper prior for beta - algorithms 1 & 2 
			    1: proper prior for beta */
  int keep = 1;          /* counter for keep */
  int progress = 1;      /* counter for progress report */
  double *beta;	         /* The model parameters */
  double **Sigma;        /* The unidentifiable variance matrix */
  double **X;            /* The covariates and outcome var */
  double **A0;           /* prior precision for beta */
  double **S;            /* The prior parameters for Sigma */

  /* model parameters */
  int **Y;                   /* The response variable */
  double **W;                /* The missing data! */
  double cmean;  	     /* The conditional mean for Wij */
  double cvar;               /* The conditional variance for Wij */
  double maxw=0.0, minw=0.0; /* used for sampling W */
  int *maxy;
  double **Xbeta;
  double ***PerSig;          /* Permuted Sigma for sampling W */
  int kindx, lindx;          /* Indexes for the permutation of Sigma */
  double *mbeta;             /* Posterior mean of beta */
  double **Vbeta;   	     /* The var of beta, given Yaug */
  double **Xstar;            /* L*X where L is the Chol factor of SigInv */
  double **SigInv;           /* The Inverse of Sigma */
  double *epsilon;           /* The error term */
  double **R;                /* Sum of squares matrix e'e */
  double **SS;	             /* Sum of squares for sweep operator */
  double alpha2;             /* alpha^2: Inv-chisq df=nu0 */
  double ss;                 /* scale parameter for alpha^2 */
  int i, j, k, l, main_loop; /* used for loops */

  /* temporay storages */
  int itemp, itempMax, itempMin, *ivtemp, itempS = 0, itempP=ftrunc((double) n_gen/10); 
  double *vtemp;
  double **mtemp,**mtemp1,**mtemp2;

  /** get random seed **/
  GetRNGstate();

  /** defining vectors and matricies **/
  Y = intMatrix(n_samp, n_dim+1);
  W = doubleMatrix(n_samp, n_dim+1);
  X = doubleMatrix(n_samp*n_dim+n_cov, n_cov+1);
  Xbeta = doubleMatrix(n_samp, n_dim);
  SS = doubleMatrix(n_cov+1, n_cov+1);
  Sigma = doubleMatrix(n_dim, n_dim);
  SigInv = doubleMatrix(n_dim, n_dim);
  epsilon = doubleArray(n_samp*n_dim);
  R = doubleMatrix(n_dim, n_dim);
  beta = doubleArray(n_cov);
  mbeta = doubleArray(n_cov);
  Vbeta = doubleMatrix(n_cov, n_cov);
  A0 = doubleMatrix(n_cov, n_cov);
  S = doubleMatrix(n_dim, n_dim);
  vtemp = doubleArray(n_dim-1);
  maxy = intArray(n_samp);
  ivtemp = intArray(n_dim+1);
  Xstar = doubleMatrix(n_samp*n_dim+n_cov, n_cov+1);
  mtemp = doubleMatrix(n_cov, n_cov);
  mtemp1 = doubleMatrix(n_dim, n_dim);
  mtemp2 = doubleMatrix(n_dim, n_dim);
  PerSig = doubleMatrix3D(n_dim, n_dim, n_dim);

  /** Packing Y, X, A0, S, beta, Sigma  **/
  if(*piMoP) {
    itemp = 0;
    for (j = 0; j <= n_dim; j++) 
      for (i = 0; i < n_samp; i++) Y[i][j] = y[itemp++];
    for (i = 0; i < n_samp; i++) { /* recording the max of Y[i] */
      for (j=0; j<= n_dim; j++)
	ivtemp[j]=Y[i][j];
      R_isort(ivtemp, n_dim+1);
      maxy[i]=ivtemp[n_dim];
    }
  }

  itemp = 0;
  for (k = 0; k < n_cov; k++) 
    for (j = 0; j < n_dim; j++) 
      for (i = 0; i < n_samp; i++) 
	X[i*n_dim+j][k] = pdX[itemp++];
      
  /* PdoubleMatrix(X, n_dim*3, n_cov); */
  itemp = 0;
  for (k = 0; k < n_cov; k++) 
    for (j = 0; j < n_cov; j++)	A0[j][k] = pdA0[itemp++];

  itemp = 0;
  for (k = 0; k < n_dim; k++) 
    for (j = 0; j < n_dim; j++) S[j][k] = pdS[itemp++];

  itemp = 0;
  for (j=0; j<n_cov;j++) 
    beta[j]=pdbeta[itemp++];

  itemp = 0;
  for (k = 0; k < n_dim; k++) 
    for (j = 0; j < n_dim; j++) 
      Sigma[j][k] = pdSigma[itemp++];
  dinv(Sigma,n_dim,SigInv); 

  /** add prior information as additional data points **/
  if(imp){
    for(j=0;j<n_cov;j++)
      for(k=0;k<=n_cov;k++)
	Xstar[n_samp*n_dim+j][k]=0;
  }
  else{
    dcholdc(A0,n_cov,mtemp); /* Cholesky decomposition */
    for(j=0;j<n_cov;j++){
      Xstar[n_samp*n_dim+j][n_cov]=b0[j];
      for(k=0;k<n_cov;k++){ 
	Xstar[n_samp*n_dim+j][n_cov]+=mtemp[j][k]*b0[k]; 
	Xstar[n_samp*n_dim+j][k]=mtemp[j][k];
      }
    }
  }

  /** starting values for W **/
  for(i=0;i<n_samp;i++) {
    for(j=0;j<n_dim;j++){
      Xbeta[i][j]=0;
      for (k=0;k<n_cov;k++) Xbeta[i][j] += X[i*n_dim+j][k]*beta[k];
      if(*piMoP) /* you DO need to worry about ordering for this */
	W[i][j] = (double)(Y[i][j]-Y[i][n_dim])/(double) n_dim;
      else /* you DON'T need to worry about ordering for this */
	W[i][j] = Xbeta[i][j] + norm_rand();
    }
    W[i][n_dim]=0.0;
  }

  /*** GIBBS SAMPLER! ***/
  for(main_loop=1; main_loop<=n_gen; main_loop++){
    /* for Algorithm 1 (Scheme 1), draw alpha2 first */
    ss=0;
    for(j=0;j<n_dim;j++)      
      for(k=0;k<n_dim;k++) mtemp1[j][k]=0;
    for(i=0;i<n_dim;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++) mtemp1[j][k]+=S[j][i]*SigInv[i][k];
    for(j=0;j<n_dim;j++) ss+=mtemp1[j][j];
    alpha2=ss/(double)rchisq((double)nu0*n_dim);

    /** permutation of Sigma **/
    for(j=0;j<n_dim;j++){
      kindx = 0;
      for(k=0;k<n_dim;k++){ 
	lindx=0;
	for(l=0;l<n_dim;l++){ 
	  if(j!=k)
	    if(j!=l)
	      PerSig[j][k+kindx][l+lindx]=Sigma[k][l];
	    else{
	      lindx=-1;
	      PerSig[j][k+kindx][n_dim-1]=Sigma[k][l];
	    }
	  else{
	    kindx=-1;
	    if(j==l){
	      lindx=-1;
	      PerSig[j][n_dim-1][n_dim-1]=Sigma[j][j];
	    }
	    else
	      PerSig[j][n_dim-1][l+lindx]=Sigma[k][l];
	  }
	}
      }
      dinv(PerSig[j],n_dim,PerSig[j]);
    }

    /** Truncated Multinomial Normal Sampling for W **/
    if(*piMoP) { /* Multinomial Probit with Ordered Preferences */
      for(i=0;i<n_samp;i++){
	for(j=0;j<n_dim;j++){
	  Xbeta[i][j]=0;
	  for (k=0;k<n_cov;k++) Xbeta[i][j]+=X[i*n_dim+j][k]*beta[k];
	} 
	for(j=0;j<n_dim;j++){
	  if(Y[i][j]!=-1) {
	    itempMax=0; itempMin=0;
	    for (k=0;k<=n_dim;k++) 	  		  		  
	      if((j!=k) && (Y[i][k]!=-1)) {
		if(Y[i][k]==(Y[i][j]+1)) {
		  if(itempMax==0) {
		    maxw=W[i][k];
		    itempMax++;
		  }
		  if((itempMax>0) && (maxw > W[i][k]))
		    maxw=W[i][k];
		}
		if(Y[i][k]==(Y[i][j]-1)) {
		  if(itempMin==0) {
		    minw=W[i][k];
		    itempMin++;
		  }
		  if((itempMin>0) && (minw < W[i][k]))
		    minw=W[i][k];
		}
	      }
	  }
	  itemp=0; 
	  for (k=0;k<n_dim;k++) 
	    if(j!=k) vtemp[itemp++]=W[i][k]-Xbeta[i][k];
	  /* conditional mean and variance */
	  cmean=Xbeta[i][j];
	  cvar=1/PerSig[j][n_dim-1][n_dim-1];
	  for(k=0;k<(n_dim-1);k++) 
	    cmean-=PerSig[j][n_dim-1][k]*vtemp[k]*cvar;
	  /* sampling each W[i][j] conditionally on the other elements */
	  /* printf("%5d%5d%5d%5d%5d%14g%14g%14g", i, Y[i][0], Y[i][1],
	     Y[i][2], Y[i][3], W[i][0], W[i][1], W[i][2]); */
	  if(Y[i][j]==-1) 
	    W[i][j] = cmean + norm_rand()*sqrt(cvar);
	  else {
	    if(Y[i][j]==0) minw = cmean - 1000*sqrt(cvar);
	    if(Y[i][j]==maxy[i]) maxw = cmean + 1000*sqrt(cvar);
	    W[i][j]=TruncNorm(minw,maxw,cmean,cvar,*invcdf); 
	  }
	  /* printf("%14g\n", W[i][j]); */
	  X[i*n_dim+j][n_cov]=W[i][j];
	  X[i*n_dim+j][n_cov]*=sqrt(alpha2);
	}
      }
    }
    else { /* standard multinomial probit */
      for(i=0;i<n_samp;i++){
	for(j=0;j<n_dim;j++) {
	  Xbeta[i][j]=0;
	  for(k=0;k<n_cov;k++) Xbeta[i][j]+=X[i*n_dim+j][k]*beta[k];
	}
	for(j=0;j<n_dim;j++){
	  maxw=0.0; itemp=0; 
	  for(k=0;k<n_dim;k++) {	  		  		  
	    if(j!=k) {
	      maxw = fmax2(maxw, W[i][k]);
	      /* if(maxw<W[i][k]) maxw=W[i][k]; */
	      vtemp[itemp++]=W[i][k]-Xbeta[i][k];
	    }
	  }
	  /* conditional mean and variance */
	  cmean=Xbeta[i][j];
	  cvar=1/PerSig[j][n_dim-1][n_dim-1];
	  for(k=0;k<(n_dim-1);k++) 
	    cmean-=PerSig[j][n_dim-1][k]*vtemp[k]*cvar;
	  /* sampling each W[i][j] conditionally on the other elements */
	  if(y[i]==-1)
	    W[i][j]=cmean+norm_rand()*sqrt(cvar);
	  else if(y[i]==(j+1)) 
	    W[i][j]=TruncNorm(maxw,cmean+1000*sqrt(cvar),cmean,cvar,*invcdf); 
	  else
	    W[i][j]=TruncNorm(cmean-1000*sqrt(cvar),maxw,cmean,cvar,*invcdf);
	  X[i*n_dim+j][n_cov]=W[i][j];
	  X[i*n_dim+j][n_cov]*=sqrt(alpha2);
	}
      }
    }
    
    /* trace(S*SigInv) */
    ss=0;
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) mtemp1[j][k]=0; 
    for(i=0;i<n_dim;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++) mtemp1[j][k]+=S[j][i]*SigInv[i][k];
    for(j=0;j<n_dim;j++) ss+=mtemp1[j][j];
    
    /* multiply X and W by the Inverse of the Cholesky factor */
    dcholdc(SigInv,n_dim,mtemp1);
    for(i=0;i<n_samp*n_dim;i++)
      for(j=0;j<=n_cov;j++) Xstar[i][j]=0;
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++) 
	  for(l=0;l<=n_cov;l++)
	    Xstar[i*n_dim+k][l]+=mtemp1[j][k]*X[i*n_dim+j][l];
    
    /* construct SS matrix for SWEEP */
    for(j=0;j<=n_cov;j++)
      for(k=0;k<=n_cov;k++) SS[j][k]=0;
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<=n_cov;k++)
	  for(l=0;l<=n_cov;l++) 
	    SS[k][l]+=Xstar[i*n_dim+j][k]*Xstar[i*n_dim+j][l];
    for(j=0;j<n_cov;j++)
      for(k=0;k<=n_cov;k++)
	for(l=0;l<=n_cov;l++) 
	  SS[k][l]+=Xstar[n_samp*n_dim+j][k]*Xstar[n_samp*n_dim+j][l];
    
    /* SWEEP to get posterior mean and variance for beta */
    for(j=0;j<n_cov;j++) SWP(SS,j,n_cov+1);
    /* draw alpha2 given Sigma and W */
    ss+=SS[n_cov][n_cov];   
    alpha2=ss/(double)rchisq((double)(n_samp+nu0)*n_dim);
    
    /* draw beta given Sigma and W */
    for(j=0;j<n_cov;j++){ 
      mbeta[j]=SS[j][n_cov];
      for(k=0;k<n_cov;k++) Vbeta[j][k]=-SS[j][k]*alpha2;
    }
    rMVN(beta,mbeta,Vbeta,n_cov);  
    /* draw Sigma given beta and W */
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++){
	epsilon[i*n_dim+j]=X[i*n_dim+j][n_cov]; 
	for(k=0;k<n_cov;k++)
	  epsilon[i*n_dim+j]-=X[i*n_dim+j][k]*beta[k];
      } 
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) R[j][k]=0;
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++) R[j][k]+=epsilon[i*n_dim+j]*epsilon[i*n_dim+k];
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) mtemp1[j][k]=S[j][k]+R[j][k];
    dinv(mtemp1,n_dim,mtemp2);
    rWish(SigInv,mtemp2,nu0+n_samp,n_dim);
    dinv(SigInv,n_dim,Sigma);
    
    /* recompute some quantities using the updated alpha2 */
    for(j=0;j<n_cov;j++) beta[j]/=sqrt(alpha2);
    if (*itrace) {
      alpha2=0;
      for(k=0;k<n_dim; k++) 
	alpha2+=Sigma[k][k];
      alpha2 = alpha2/n_dim;
    } else {
      alpha2=Sigma[0][0];
    }
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) {
	Sigma[j][k]/=alpha2;
	SigInv[j][k]*=alpha2;
      }
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	W[i][j]=X[i*n_dim+j][n_cov]/sqrt(alpha2);
  
    /* print Gibbs draws for all the schmes! */
    R_CheckUserInterrupt();
    if(main_loop > *piBurnin) {
      if(keep==*piKeep) {
	for(j=0;j<n_cov;j++) 
	  pdStore[itempS++]=beta[j];
	for(j=0;j<n_dim;j++) 
	  for(k=0;k<n_dim;k++) 
	    if(j<=k) 
	      pdStore[itempS++]=Sigma[j][k];
	if (*latent) 
	  for (i=0;i<n_samp;i++) 
	    for (j=0;j<n_dim;j++) 
	      pdStore[itempS++]=W[i][j];
	keep=1;
      }
      else
	keep++;
    }
    if(*verbose) {
      if(main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	itempP+=ftrunc((double) n_gen/10); progress++;
	R_FlushConsole(); 
      }
    }
  } /* end of Gibbs sampler */
  
  /** write out the random seed **/
  PutRNGstate();
  
  /** freeing memory **/
  if(*piMoP)
    FreeintMatrix(Y, n_samp);
  FreeMatrix(W, n_samp);
  FreeMatrix(X, n_samp*n_dim+n_cov);
  FreeMatrix(Xbeta, n_samp);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(Sigma, n_dim);
  FreeMatrix(SigInv, n_dim);
  free(epsilon);
  FreeMatrix(R, n_dim);
  free(beta);
  free(mbeta);
  FreeMatrix(Vbeta, n_cov);
  FreeMatrix(A0, n_cov);
  FreeMatrix(S, n_dim);
  free(vtemp);
  free(maxy);
  free(ivtemp);
  FreeMatrix(Xstar, n_samp*n_dim+n_cov);
  FreeMatrix(mtemp, n_cov);
  FreeMatrix(mtemp1, n_dim);
  FreeMatrix(mtemp2, n_dim);
  Free3DMatrix(PerSig, n_dim, n_dim);

} /* main */


/* unitility function */
void R_max_col2(double **matrix, 
	       int nr, 
	       int nc, 
	       int *maxes, 
	       int ties_meth) {
  
  int *ncol = intArray(1);
  int *nrow = intArray(1);
  int *ties = intArray(1);
  int *itmp = intArray(1);
  double *tmp = doubleArray(nr*nc);
  int i, j, k;

  ncol[0] = nc; nrow[0] = nr; ties[0] = ties_meth;
  i = 0;
  for (k = 0; k < nc; k++)
    for (j = 0; j < nr; j++)
      tmp[i++] = matrix[j][k];

  R_max_col(tmp, nrow, ncol, maxes, ties);

  free(ncol);
  free(nrow);
  free(itmp);
  free(tmp);

}
 

void predict(double *dX,     /* X matrix */
	     int *nobs,      /* number of observations */
	     double *dcoef,  /* coefficients */
	     double *dSigma, /* covariances */
	     int *ndims,     /* number of dimensions */
	     int *ncovs,     /* number of covariates */
	     int *ndraws,    /* number of MCMC draws */
	     int *moredraws, /* number of extra draws */
	     int *verbose,
	     double *prob,   /* probability output */
	     double *choice, /* choice output */
	     double *order  /* order output */
	     ) {

  int n_samp = *nobs;
  int n_draw = *ndraws;
  int n_dim = *ndims;
  int n_cov = *ncovs;
  int n_extra = *moredraws; 

  double **X = doubleMatrix(n_samp*n_dim, n_cov);
  double *Xbeta = doubleArray(n_cov);
  double *vtemp = doubleArray(n_dim);
  double **W = doubleMatrix(n_extra, n_dim+1);
  double **beta = doubleMatrix(n_draw, n_cov);
  double **mtemp = doubleMatrix(n_dim, n_dim);
  double ***Sigma = doubleMatrix3D(n_draw, n_dim, n_dim);

  int i, j, k, main_loop, itemp, itempP, itempO, itempC;
  int total = n_extra*n_samp*n_draw;
  int count, progress = 1; 
  int itempQ = ftrunc((double) total/10);
  int *maxdim = intArray(n_extra);
  int *ind = intArray(n_dim+1);
  int *sumorder = intArray(n_dim+1);
  int *probTemp = intArray(n_dim+1);

  /** reading the data */
  itemp = 0;
  for (k = 0; k < n_cov; k++) 
    for (j = 0; j < n_dim; j++) 
      for (i = 0; i < n_samp; i++) 
	X[i*n_dim+j][k] = dX[itemp++];  

  /** reading the MCMC draws **/
  itemp = 0;
  for (k = 0; k < n_cov; k++)
    for (j = 0; j < n_draw; j++)
      beta[j][k] = dcoef[itemp++];
  
  itemp = 0;
  for (k = 0; k < n_dim; k++)
    for (j = k; j < n_dim; j++)
      for (i = 0; i < n_draw; i++) {
	Sigma[i][j][k] = dSigma[itemp++];
	Sigma[i][k][j] = Sigma[i][j][k];
      }

  /** get random seed **/
  GetRNGstate();

  /** Posterior predictive simulations **/
  itempC = 0; itempO = 0; itempP = 0; count = 0; 
  itempQ = 0;

  /* loop for observations */
  for (i = 0; i < n_samp; i++) { 
    if (n_extra == 1) {
      for (j = 0; j <= n_dim; j++) 
	probTemp[j] = 0;
    }
    /* loop for MCMC draws */
    for (main_loop = 0; main_loop < n_draw; main_loop++) {
      if (n_extra > 1) {
	for (j = 0; j <= n_dim; j++) 
	  probTemp[j] = 0;
      }
      /* compute the mean for each dimension */
      for (j = 0; j < n_dim; j++) {
	Xbeta[j] = 0;
	for (k = 0; k < n_cov; k++)
	  Xbeta[j] += X[i*n_dim+j][k] * beta[main_loop][k];
      }
      /* PdoubleArray(Xbeta, n_dim); */
      /* sample W */
      for (j = 0; j < n_extra; j++) {
	/*dinv(Sigma[main_loop], n_dim, mtemp);*/
	rMVN(vtemp, Xbeta, Sigma[main_loop], n_dim);
	for (k = 0; k < n_dim; k++)
	  W[j][k+1] = vtemp[k];
	W[j][0] = 0;
      }
      /* which dimension is max for each of n_extra W draws? 
	 PdoubleMatrix(W, n_extra, n_dim+1); */
      R_max_col2(W, n_extra, n_dim+1, maxdim, 1);
      /* PintArray(maxdim, n_extra); */
      /* order */
      for (j = 0; j < n_extra; j++) {
	for (k = 0; k <= n_dim; k++) {
	  ind[k] = k; sumorder[k] = 0;
	}
	revsort(W[j], ind, n_dim+1);
	/* PintArray(ind, n_dim+1); */
	for (k = 0; k <= n_dim; k++)
	  sumorder[ind[k]] += (k+1);
	if(*verbose) {
	  if(count == itempQ) {
	    Rprintf("%3d percent done.\n", progress*10);
	    itempQ += ftrunc((double) total/10); progress++;
	    R_FlushConsole(); 
	  }
	  count++;
	}
      }
      if (n_extra > 1) {
	/* store probability and mean order */
	for (j = 0; j <= n_dim; j++) {
	  itemp = 0;
	  for (k = 0; k < n_extra; k++)
	    if (maxdim[k] == (j+1))
	      itemp++;
	  prob[itempP++] = ((double) itemp / (double) n_extra);
	  order[itempO++] = ((double) sumorder[j] / (double) n_extra);
	}
      } else {
	/* store choice */
	for (j = 0; j <= n_dim; j++) {
	  if (maxdim[0] == (j+1)) {
	    choice[itempC++] = j;
	    probTemp[j]++;
	  }
	  order[itempO++] = sumorder[j];
	}
      }
    }
    if (n_extra == 1)
      for (j = 0; j <= n_dim; j++)
	prob[itempP++] = ((double) probTemp[j] / (double) n_draw);
  }

  /** write out the random seed **/
  PutRNGstate();

  /* freeing memory */
  FreeMatrix(X, n_samp*n_dim);
  free(vtemp);
  free(Xbeta);
  FreeMatrix(W, n_extra);
  FreeMatrix(beta, n_draw);
  FreeMatrix(mtemp, n_dim);
  Free3DMatrix(Sigma, n_draw, n_dim);
  free(maxdim);
  free(ind);
  free(sumorder);
  free(probTemp);
}
