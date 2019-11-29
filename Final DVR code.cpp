#include <iostream>
#include <complex>
#include <cmath>
#include<fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

const  double stretchX = 1.8e-10;
const  double stretchY = 1.95e-10;

const double  dx = 0.075e-10;
const double  dy = 0.075e-10;

const int gridsizeX = 2.0*(stretchX/dx) + 1;
const int gridsizeY = 2.0*(stretchY/dx) + 1;
const int K = gridsizeX*gridsizeY;

//For the DVR matrix:
double pi2 = 3.1415926;
double hcross2  = 1.0545718e-34;
double massHydrogen2 = 1.6735575e-27;
double mass = 1.6735575e-27;
double KineticHamiltonianConst = (hcross2*hcross2)/(2.0*massHydrogen2*dx*dx);

//XY values and Psy;
MatrixXd XYvalues(K,2) ;
MatrixXcd PsyK (K,1);

void XYkFill (MatrixXd(&K));
void PsyFillk (MatrixXd MxyVal);
double Vofk ( double k);
double Delta( int i, int j);
double PotentialSurface (double x, double y);

complex<double> iota (0,1);
const double omega = 1.0e14;

//************************************************************************************

double Vc = 9.932237219250332e-022 ;
double eps1 =1.986447443850067e-022;
double eps2 =-1.986447443850066e-021;


double omg1 =    37673031346177.1;
double omg2 =   39556682913485.9;
double omg3 = 37673031346177.1;





int main()
{

 XYkFill (XYvalues); //XY values are filled in a linear K form
 PsyFillk (XYvalues);// Psy is filled in a linear K form


MatrixXd Hamiltonian(K,K);
MatrixXd DVRx(gridsizeX,gridsizeX);

for (int i=0;i<gridsizeX;i++)
     {
       for (int j=0;j<gridsizeX;j++)
                {
                  if (j==i) {DVRx(i,j) = (KineticHamiltonianConst* pow(-1,i-j)*pi2*pi2)/(3.0);}
                   else     {DVRx(i,j) = (KineticHamiltonianConst* pow(-1,i-j)*(2/pow(i-j,2))); }
                  }
         }

MatrixXd DVRy(gridsizeY,gridsizeY);

for (int i=0;i<gridsizeY;i++)
     {
       for (int j=0;j<gridsizeY;j++)
                {
                  if (j==i) {DVRy(i,j) = (KineticHamiltonianConst* pow(-1,i-j)*pi2*pi2)/(3.0);}
                   else     {DVRy(i,j) = (KineticHamiltonianConst* pow(-1,i-j)*(2/pow(i-j,2))); }
                  }
         }
/*
for (int k=0;k<K;k++) //DEBUGGER LOOP
     {

         //here i is the k
         // j+1 and ii+1 are the answers for indexing starting from 1.
         //j and ii are answers for indexing starting from 0;
         //k starts from 0;
        int j = (k)% gridsizeY;
        int i =(k-j)/ gridsizeY;
       // cout << k+1 <<". "<< XYvalues(k,0) << " " << XYvalues(k,1) <<" i= " << i  << " j="<<j <<endl ;
         cout << k+1 <<". "<< XYvalues(k,0) << " " << XYvalues(k,1) <<" "<<Vofk(k)<<endl;
     } */

//***************************************************************************************************//
   for (int k=0;k<K;k++)  // Filling the Hamiltonian
   {
       for (int kk=0;kk<K;kk++)

       {

           int j = (k)% gridsizeY;
           int i =(k-j)/ gridsizeY;

           int jj = (kk)% gridsizeY;
           int ii =(kk-jj)/ gridsizeY;

           //if (k==kk){ Hamiltonian(k,kk) = Vofk(k) + DVRx(i,ii) + DVRy(j,jj); }
           //else {Hamiltonian(k,kk)= Delta(j,jj)*DVRx(i,ii) + Delta(i,ii)*DVRy(j,jj);}

           if (k==kk){ Hamiltonian(k,kk) = PotentialSurface(XYvalues(k),XYvalues(k)) + DVRx(i,ii) + DVRy(j,jj); }
           else {Hamiltonian(k,kk)= Delta(j,jj)*DVRx(i,ii) + Delta(i,ii)*DVRy(j,jj);}

           }
   }

//****************************************************************************************
//DIAGONALIZATOIN:

SelfAdjointEigenSolver<MatrixXd> eigensolver(Hamiltonian);
   if (eigensolver.info() != Success) abort();
   else { cout<< "Success"<<endl;}


//****************************************************************************************
//Decomposition of Psy as a linear combination of its basis

MatrixXcd PsyColumn = PsyK.transpose();
MatrixXcd Constants = PsyColumn*eigensolver.eigenvectors();
//cout << Constants;   // Constants are a column matrix.
MatrixXcd DecomposedPsy;
MatrixXcd PsyAtAnyT;
MatrixXcd ConstantsWithTime(1,K);

//*******************************************************
//Storing values in a matrix
double t = 0.0;   //Deity
int snapshots =10;  //Photu
double dt = 0.02e-15; //samay

MatrixXd Storage (K,snapshots+3);

for (int i=0;i<K;i++)
{
     Storage(i,0) = XYvalues(i,0);
     Storage(i,1) = XYvalues(i,1);
}

for (int i=0;i<K;i++)
{
     Storage(i,2) = PotentialSurface(XYvalues(i,0),XYvalues(i,1));

}




for (int s =0;s < snapshots; s++)
{
     for (int i=0;i<K;i++)
     {ConstantsWithTime(0,i) = Constants(0,i)*exp((-1.0*iota*eigensolver.eigenvalues()(i,0)*t)/hcross2);
      t=t+dt;
        }
     PsyAtAnyT = eigensolver.eigenvectors()*ConstantsWithTime.transpose();

     for (int ii=0;ii<K;ii++)

    {

       // FileM << abs(PsyAtAnyT.transpose()(0,ii))<< endl;
       Storage(ii,s+3) = abs(PsyAtAnyT.transpose()(0,ii));
    }




   }






//**************************************************************************************
//From Storgae Matrix to a file grpahical.


ofstream PsyAtTimes ("D.dat"); //NAAM
for (int i=0;i<K;i++)
{
    for (int j=0;j<snapshots+3;j++)
  {
     PsyAtTimes << Storage (i,j) << " ";
  }
  PsyAtTimes <<endl;
  if((i+1)%gridsizeY==0&&i!=0) {PsyAtTimes << endl;}
}


//****************************************************************************************
/*
for (int i=0;i<K;i++)
{
    cout << Constants(0,i)  <<endl;
} */
//*****************************************************************************************
//Finding reactant and product populations

ofstream Ratios ("RatiosD.csv") ;  // NAAM2

MatrixXcd ConstantsWithTime2(1,K);
MatrixXcd PsyAtAnyT2;

for (double t2=0.0 ;t2<= 300e-15;t2= t2 + 1.0e-15)
{


for (int i=0;i<K;i++)
{ConstantsWithTime2(0,i) = Constants(0,i)*exp((-1.0*iota*eigensolver.eigenvalues()(i,0)*t2)/hcross2);}

PsyAtAnyT2 = eigensolver.eigenvectors()*ConstantsWithTime2.transpose();
double reactants=0;
double products=0;
for (int i=0;i<K;i++)
{
 if (XYvalues(i,0) < 0) {

                            reactants = reactants + abs(PsyAtAnyT2(i,0))*abs(PsyAtAnyT2(i,0))*dx*dy;
                                  }

 else {
         products = products + abs(PsyAtAnyT2(i,0))*abs(PsyAtAnyT2(i,0))*dx*dy;
             }

}
double ratio = products/reactants;

//cout<< products + reactants <<endl;

Ratios << t2 <<","<< ratio<<endl;

}

} //int main

double Vofk ( double k)
{

        //  return XYvalues(k,0) + XYvalues(k,1); //x + y
          return 0.5*massHydrogen2*(XYvalues(k,1)*XYvalues(k,1) +XYvalues(k,0)*XYvalues(k,0))*omega*omega;
        // return 0;
              }

void XYkFill (MatrixXd(&K))
{

  double tempX= -1*stretchX;
  double tempY= -1*stretchY;


 for (int k=0;k< gridsizeX*gridsizeY ; k++)

{        K(k,0) = tempX;
         K(k,1) = tempY;
         tempY= tempY + dy;

         if(((k+1)%gridsizeY)==0) { tempX = tempX + dx; tempY= -1*stretchY;}


          }

}

void PsyFillk (MatrixXd MxyVal)
{

for (int k=0;k< K;k++)
 {

 //PsyK[k] = K[i].x + K[i].y;
//  PsyK(k,0) = MxyVal(k,0) + MxyVal(k,1) ;   //x + y;
//psi(x) = (mw/(pi\hbar))^1/4 * exp(-mwx^2/(2\hbar))
double Nml = (massHydrogen2*omg1)/(2.0*pi2*hcross2);

double Psy1 = exp((-1.0*massHydrogen2*omg1*(((MxyVal(k,0)+5.50844e-011)*(MxyVal(k,0)+5.50844e-011))+((MxyVal(k,1)-3.16795e-011)*(MxyVal(k,1)-3.16795e-011))))/(2.0*hcross2));
double Psy2 = exp((-1.0*massHydrogen2*omg2*(((MxyVal(k,0)+5.50844e-011)*(MxyVal(k,0)+5.50844e-011))+((MxyVal(k,1)+3.16795e-011)*(MxyVal(k,1)+3.16795e-011))))/(2.0*hcross2));
PsyK(k,0) = Nml*(Psy1 + Psy2);//SAI
//;
 }// findPsy



}

double Delta( int i, int j)

{
     if (i==j){return 1.0;}
     else {return 0.0;}
}



double PotentialSurface (double x, double y)
{

double Vb12 =   1.191868466310040e-021;
double Vb13= 5.363408098395179e-021;
double Vb3 = 6.356631820320213e-021;



double y0=sqrt(Vb12*2/(mass*omg1*omg1));
double x01=-sqrt((Vb13-Vb12-eps1)*2/(mass*omg2*omg2));
double x02=sqrt((Vb3-eps2)*2/(mass*omg3*omg3));

MatrixXd PotentialSurface(3,3);

for (int i =0;i<3;i++)
for (int j=0;j<3;j++)
 {

   if(i!=j)  { PotentialSurface(i,j)  = Vc; }
 }

PotentialSurface(0,0)=0.5*mass*omg1*omg1*(x-x01)*(x-x01)+0.5*mass*omg1*omg1*(y-y0)*(y-y0);
PotentialSurface(1,1)=0.5*mass*omg2*omg1*(x-x01)*(x-x01)+0.5*mass*omg2*omg2*(y+y0)*(y+y0)+eps1;
PotentialSurface(2,2)=0.5*mass*omg3*omg3*(x-x02)*(x-x02)+0.5*mass*omg3*omg3*y*y+eps2;



SelfAdjointEigenSolver<MatrixXd> eigenPotential(PotentialSurface);

//cout << eigenPotential.eigenvalues()<<endl;

double minima1 = min (eigenPotential.eigenvalues()(0,0),eigenPotential.eigenvalues()(1,0));
double minima2 = min (eigenPotential.eigenvalues()(1,0),eigenPotential.eigenvalues()(2,0));
double minima3 = min (minima1,minima2);
return minima3;
}

