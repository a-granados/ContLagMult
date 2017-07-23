#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>  /*per sin i cos*/
#include "rk78.h"
#include "memoria.h"
#include "lu.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


void vfield(double t, double *x, int ndim, double *dx);
void vfield_varAb(double t, double *x, int ndim, double *dx);
void F(double *xa, double *Fval, double **DF);
int Newton_cond_bound(double *x,double **DF, int m, void (*F)(double *xa,double *Fval, double **DF));
void compute_Deltaw(int n, int m, double **A, double *b,double *Deltaw);
void kernel (double *v,double **A,int m);

using std::cout;
//using std::cerr;
using std::endl;

int main(int argc, char * argv[])
{
  cout.precision(20);
  double DeltaS;
  int res,i;
  int n=2;//Number of equations
  int m=n+1;//Number of unknowns=n+1
  double *x=new double[m];
  double *xprev=new double[m];
  double *v=new double[m];
  double **DF=new double*[n];
  for (i=0;i<n;i++){
    DF[i]=new double[m];
  }

  //Initial guess:
  x[0]=5;
  x[1]=0;
  x[2]=0;

  res=Newton_cond_bound(x,DF,m,F);
  if (res==0){
    cout <<"-1 ";
    for (i=0;i<m;i++){
      cout <<x[i]<<" ";
    }
    cout <<endl;
  }
  DeltaS=0.1;
  kernel(v,DF,m);//Returns a vector belonging to the kernel of DG, with modulus 1.
  while (fabs(DeltaS)>1e-8){
    for (i=0;i<m;i++){
      xprev[i]=x[i];
      x[i]+=DeltaS*v[i];
    }
    /*
    //Plot the predicted point:
    cout <<"-2 ";
    for (i=0;i<m;i++){
      cout <<x[i]<<" ";
    }
    cout <<endl;
    */
    res=Newton_cond_bound(x,DF,m,F);
    cout <<"#res= "<<res<<endl;
    if (res==0){
    //if (res==0){
      cout <<"-1 ";
      for (i=0;i<m;i++){
	cout <<x[i]<<" ";
      }
      cout <<endl;
      DeltaS=1.01*DeltaS;
      /*
      cout <<"#DF: "<<endl;
      int j;
      for (i=0;i<n;i++){
	for (j=0;j<m;j++){
	  cout <<DF[i][j]<<" ";
	}
	cout <<endl;
      }
      */
      kernel(v,DF,m);
      //Choose the direction of the vector:
      if (v[0]>0){
	for (i=0;i<m;i++){
	  v[i]=-v[i];
	}
      }
      /*
      cout <<"#Eigvector:"<<endl;
      for (i=0;i<m;i++){
	cout <<v[i]<<" ";
      }
      cout <<endl;
      */
    }
    else{
      for(i=0;i<m;i++){
	x[i]=xprev[i];
      }
      //DeltaS=DeltaS/2.;
      DeltaS=0.9*DeltaS;
    }
  }
  
  for (i=0;i<n;i++){
    delete[] DF[i];
  }
  delete[] DF;
  delete[] x;
  delete[] xprev;
  delete[] v;
}

void F(double *xa, double *Fval, double **DF){
  //This evaluates the function we whish to solve using conditional
  //boundaries to continue the curve. It has to return the value of
  //the function (an m-vector) in Fval and the differential in DF.
  double tw = 5;
  double ts = 1;
  double gsw = -2;
  double gws = -2;
  double Wmax = 6.5;
  double bw = -0.3;
  double aw = 0.5;
  double Smax = 5;
  double as = 0.25;
  double ys = 4;
  double yw = 5;

  double fw=xa[0];
  double fs=xa[1];
  double H=xa[2];

  Fval[0]=(1/tw)*(-fw+((Wmax/2)*(1+tanh((gsw*tanh(fs/ys)-bw)/aw))));
  Fval[1]=(1/ts)*(-fs+((Smax/2)*(1+tanh((gws*tanh(fw/yw) + 1.5*H)/as))));
  /*
  DF[0][0]=-1/tw;
  DF[0][1]=(1/2)*Wmax*gsw*(1-pow(tanh(fs/ys),2))*(1-pow(tanh((gsw*tanh(fs/ys)-bw)/aw),2))/(tw*ys*aw);
  DF[0][2]=0;
  DF[1][0]=(1/2)*Smax*gws*(1-pow(tanh(fw/yw),2))*(1-pow(tanh((gws*tanh(fw/yw)+1.5*H)/as),2))/(ts*yw*as);
  DF[1][1]=-1/ts;
  DF[1][2]=.75*Smax*(1-pow(tanh((gws*tanh(fw/yw)+1.5*H)/as),2))/(ts*as);
  */
DF[0][0] = -1/tw;
DF[0][1]= (0.5)*Wmax*gsw*(1-pow(tanh(fs/ys),2))*(1-pow(tanh((gsw*tanh(fs/ys)-bw)/aw),2))/(tw*ys*aw);
DF[1][0]=(0.5)*Smax*gws*(1-pow(tanh(fw/yw),2))*(1-pow(tanh((gws*tanh(fw/yw)+1.5*H)/as),2))/(ts*yw*as);
DF[1][1]= -1/ts;
DF[0][2] = 0;
//DF[1][2] =(1/2)*Smax*1.5*(1-pow(tanh(gws*tanh((fw/yw)+1.5*H)/as),2))/(as*ts);
DF[1][2] =(0.5)*Smax*1.5*(1-pow(tanh((gws*tanh(fw/yw)+1.5*H)/as),2))/(as*ts);
  //DF[1][2]=.75*Smax*(1-pow(tanh((gws*tanh(fw/yw)+1.5*H)/as),2))/(ts*as);
}

int Newton_cond_bound(double *x,double **DF, int m, void (*F)(double *xa,double *Fval, double **DF)){
  //This performs a Newton method with conditioned boundaries: finds
  //to the closest point to x (by minimizing least squares) belonging
  //to the bifurcation curve.
  //The function returns 
  //1 if Newton method diverges
  //0 if a feasible fixed point was found
  int *perm;
  int n=m-1;//Number of equations
  int nite=0;
  double *xa=new double[m];
  double *Deltaw=new double[m];
  double *Fval=new double[n];
  double dist;
  double Newtol=1e-10;
  int maxiter=50;
  int i,j;
  
  for (i=0;i<m;i++){
    xa[i]=x[i];
  }
  F(xa,Fval,DF);
  dist=0;
  for (i=0;i<n;i++){
    dist+=Fval[i]*Fval[i];
  }
  dist=sqrt(dist);

  while (dist>Newtol &&nite<maxiter){
    compute_Deltaw(n,m,DF,Fval,Deltaw);
    for (i=0;i<m;i++){
      xa[i]=xa[i]+Deltaw[i];
    }
    F(xa,Fval,DF);
    /*
    cout <<"#xnext "<<endl;
    for (i=0;i<m;i++){
      cout <<xa[i]<<" ";
    }
    cout <<endl;
    cout <<"#Fval "<<endl;
    for (i=0;i<n;i++){
      cout <<Fval[i]<<" ";
    }
    cout <<endl;
    cout <<"#DF: "<<endl;
    for (i=0;i<n;i++){
      for (j=0;j<m;j++){
	cout<<DF[i][j]<<" ";
      }
      cout <<endl;
    }
    */
    dist=0;
    for (i=0;i<n;i++){
      dist+=Fval[i]*Fval[i];
    }
    dist=sqrt(dist);
    //cout<<"Dist= "<<dist<<endl;
    nite++;
  }
  /*
  if (nite==maxiter && dist > Newtol){
   cout<<"#massa iteracions newton"<<endl;
   exit(1);
  }
  */

  delete [] Fval;
  delete [] Deltaw;
  if (nite==maxiter || dist> Newtol){
    return 1;
    delete[] xa;
  }
  else {
    for (i=0;i<m;i++){
      x[i]=xa[i];
    }
    delete[] xa;
    return 0;
  }
}

void compute_Deltaw(int n, int m, double **A, double *b,double *Deltaw){
  //This computes the correction:
  //Deltaw=-A^T (AA^T)^-1b
  //n number of equations
  //m=m+1 number of unkowns
  double **tmp=new double*[n];
  double **AAT=new double*[n];
  double **LU=new double*[n];
  double **invAAT=new double*[n];
  double *baux=new double[n];
  double *aux=new double[n];
  int *perm=new int[n];
  int i,j,k,result;
  double tol=1e-10;
  for (i=0;i<n;i++){
    tmp[i]=new double[m];
    AAT[i]=new double[n];
    LU[i]=new double[n];
    invAAT[i]=new double[n];
  }
  //We compute AAT=A\cdot A^T
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      AAT[i][j]=0.;
      for (k=0;k<m;k++){
	AAT[i][j]+=A[i][k]*A[j][k];
	LU[i][j]=AAT[i][j];
      }
    }
  }
  //We now invert AAT
  for (i=0;i<n;i++){
    perm[i]=i;
  }
  result=lu(LU,n,perm,tol);
  if (result==0){
    cout <<"LU failed"<<endl;
    exit(1);
  }
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      if (j==i){
	baux[j]=1.;
      }
      else{
	baux[j]=0.;
      }
    }
    resol(LU,aux,baux,n,perm);
    for (j=0;j<n;j++){
      if (fabs(aux[j])<1e-13){
	aux[j]=0;
      }
      invAAT[j][i]=aux[j];
    }
  }
  //We now multiply inv(AAT) by b:
  for (i=0;i<n;i++){
    baux[i]=0.;
    for (j=0;j<n;j++){
      baux[i]+=invAAT[i][j]*b[j];
    }
  }
  //We finally multiply AT by baux:
  for (i=0;i<m;i++){
    Deltaw[i]=0.;
    for (j=0;j<n;j++){
      Deltaw[i]+=A[j][i]*baux[j];
    }
    Deltaw[i]=-Deltaw[i];
  }

  for (i=0;i<n;i++){
    delete[] tmp[i];
    delete[] AAT[i];
    delete[] LU[i];
    delete[] invAAT[i];
  }
  delete[] tmp;
  delete[] AAT;
  delete[] LU;
  delete[] invAAT;
  delete [] baux;
  delete [] perm;
  delete [] aux;
}

void kernel (double *v,double **A,int m){
  gsl_matrix *M=gsl_matrix_alloc(m,m);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (m);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (m, m);
  int i,j;
  double modv;

  for (i=0;i<m-1;i++){
    for (j=0;j<m;j++){
      gsl_matrix_set(M,i,j,A[i][j]);
    }
  }
  //The matrix has indeed dimension m-1,m, so we complete the last row
  //with zeros:
  for (j=0;j<m;j++){
    gsl_matrix_set(M,m-1,j,0);
  }


  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (m);
  
  gsl_eigen_nonsymmv (M, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);
  gsl_matrix_free(M);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_ASC);
  //The first eigenvector is the one of the kernel:
  //cout <<"# Vap (should be 0): "<<GSL_REAL(gsl_vector_complex_get(eval,0))<<" +"<<GSL_IMAG(gsl_vector_complex_get(eval,0))<<" i"<<endl;
  modv=0;
  for (i=0;i<m;i++){
    //v[i]=GSL_REAL(gsl_vector_complex_get(evec,i));
    v[i]=GSL_REAL(gsl_matrix_complex_get(evec,i,0));
    modv+=v[i]*v[i];
  }
  modv=sqrt(modv);
  for (i=0;i<m;i++){
    v[i]=v[i]/modv;
  }
  
  /*
  {

    for (i = 0; i < m; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %g + %gi\n",
         GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
	//cout <<GSL_REAL(eval_i)<<" + "<<GSL_IMAG(eval_i)<<" i"<<endl;
        for (j = 0; j < m; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
         printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
          }
      }
  }
  */
}

