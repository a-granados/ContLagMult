
#include <math.h>


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

