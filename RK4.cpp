#include <iostream>
#include <complex>
#include <cmath>
#include<fstream>

using namespace std;

const double  dx = 0.2e-10;
const  double stretch = 9.0e-10;
const double omega = 1.0e14;
const int gridsize = 2.0*(stretch/dx) + 1;

const int numberofsnapshots=10; //including the one at zero

complex<double> PsyX [gridsize] ;
complex<double> Angstrom (1e-10,0);
complex<double> pi (3.1415926,0);
complex<double> iota (0,1);
complex<double> PsyDoti[gridsize];
complex<double> massHydrogen (1.6735575e-27,0 );
complex<double> hcross (1.0545718e-34,0);
double hcross2  = 1.0545718e-34;
double massHydrogen2 = 1.6735575e-27;
double Angstrom2 = 1.0e-10;
double pi2 = 3.1415926;

complex<double> ihcross (0,6.582119e-16);
double KineticHamiltonianConst = (hcross2*hcross2)/(2.0*massHydrogen2*dx*dx);
double KineticHamiltonianOperator [gridsize][gridsize];
double PotentialHamiltonianOperator [gridsize] ;

complex<double> one (1.0,0.0);
complex<double> zero (0.0,0.0);

void PsyFillx (double Xvalues[],int tsize);
void PotentialHamiltonianX( double(&Vxy)[gridsize],double Xvalues[gridsize],int gsize);
void Kfunction(complex<double>(&Ki)[gridsize],complex<double> PsyX[],complex<double> kprev[],double Kdenomfactor,double dt);
void KineticHamiltonianFiller (double(&K)[gridsize][gridsize]);
void XvaluesFillfunctionOfGridsize (double(&K)[gridsize]);


int main(){

    double Xvalues[gridsize];
    XvaluesFillfunctionOfGridsize(Xvalues);
    PsyFillx(Xvalues,gridsize);

    ofstream Mydata2;
    Mydata2.open("CodeBlocks1.csv");


    //Filling Xvalues(1D grid):   (CAN BE USED IN PSI FILL)
    /*for (int i=0;i<gridsize;i++)
    {
    double ii = i;  //turning i which is int to long double so that it could be used in the formula below
    Xvalues[i]= -1.0*(10-ii);
    }*/

     //Creates Psy using X values
    PotentialHamiltonianX(PotentialHamiltonianOperator,Xvalues,gridsize); //POTENTIAL ENERGY HAMILTONIAN
    KineticHamiltonianFiller(KineticHamiltonianOperator);

complex<double>PsyT[gridsize];
complex<double>PsyTimeEvolving[gridsize][numberofsnapshots+1]; //number of snapshots include the snap at t=zero ,+1 for Xvalues.

for (int i=0;i<gridsize;i++){
                             PsyTimeEvolving[i][0] = Xvalues[i];
                             PsyTimeEvolving[i][1] = PsyX[i];
                             PsyT[i] = PsyX[i];
                        //     cout << Xvalues[i] << "  " << PsyTimeEvolving[i][0] << " " <<PsyTimeEvolving[i][1] <<endl;
                            // cout << real(PsyTimeEvolving[i]) <<endl;

                               }






complex<double> kzero[gridsize];
for (int i=0;i<gridsize;i++){   // FILLING KZERO

                             kzero[i] = zero;

                               }

complex<double> kone[gridsize];
complex<double> ktwo[gridsize];
complex<double> kthree[gridsize];
complex<double> kfour[gridsize];
complex<double> kweighted[gridsize];

double dt =0.01e-15;
int timesteps = 200;

///BIG FOR LOOP:
int snapcount =2;
int temp =0;


for (int ii = 0 ;ii <numberofsnapshots-1 ; ii++  ) //Total number of snaps it generates = number of snapshots -snapshot at t=0;
{

for(int rkfour = 0;rkfour< timesteps ; rkfour++) //SNAPSHOT AFTER TIME t GENERATOR
{

Kfunction (kone,PsyT,kzero,1.0,dt);
Kfunction (ktwo,PsyT,kone,2.0,dt);
Kfunction (kthree,PsyT,ktwo,2.0,dt);
Kfunction (kfour,PsyT,kthree,1.0,dt);




    for (int i=0;i<gridsize;i++){

                          kweighted[i] = ((kone[i]+2.0*ktwo[i]+2.0*kthree[i]+kfour[i])/6.0);


                         PsyT[i] = PsyT[i] + kweighted[i];


                                 }

}
for (int i=0; i<gridsize; i++)
{PsyTimeEvolving[i][snapcount] = PsyT[i];}

snapcount++;
temp++;
 /*for (int i=0;i<gridsize;i++)
    {


cout<< abs(PsyT[i]) << endl;
//Mydata2 << abs(PsyTimeEvolving[i]) << "," ;

       }*/
//cout << endl << "NEXT SET OF DATA" << endl;
//Mydata2 <<endl;


}
//BIG FOR LOOP


//cout<< " THIS IS ABS(PSYX at t=0)" << endl;



 //double Areas[gridsize] ;
 double temp1 = 0.0;

//for (int j=0;j<numberofsnapshots;j++)
//{

//cout << "NEXT DATA" << endl;
for (int i =0; i<gridsize;i++)
{
for (int j=0 ; j< numberofsnapshots +1 ; j++)
 {

if(j==0){
Mydata2  << real(PsyTimeEvolving[i][j]) << ","; }



else{
Mydata2  << abs(PsyTimeEvolving[i][j]) << ","; }
//Mydata2 << abs(PsyTimeEvolving[i]) << "," ;

        }
        Mydata2 <<endl;
//Mydata2  <<Areas[i] <<"," << endl;
 }
cout << "Just Before INTMAIN" <<endl;

}// int main


void PsyFillx (double Xvalues[],int tsize){ //this function is designed to fill psy

                                                   for (int i=0;i<tsize;i++)
                                                   {
                                                           complex<double> tempo (0.25,0.25);
                                                      // PsyX [i] = 10000.0*exp(-1.0*(Xvalues[i]/Angstrom2)*(Xvalues[i]/Angstrom2));
                                                     //  double normalized = pow( (massHydrogen2*2.0e14)/(pi2*hcross2) ,0.25);
                                                    //  PsyX [i] = 10000.0*exp(-1.0*(Xvalues[i]-0.5)*(Xvalues[i]-0.5)); //WORKS
                                                      //  PsyX [i] = normalized*exp((-1.0*(Xvalues[i])*(Xvalues[i])*massHydrogen2*2.0e14)/(2*hcross2));/WIGGLES
                                                        //PsyX [i] = pi;
                                                        //PsyX [i] = 10000.0*exp(-1.0*(Xvalues[i])*(Xvalues[i]));
                                                     //  PsyX [i] = sqrt(2.0*massHydrogen2*omega/hcross2)*(Xvalues[i])*exp((-1.0*(Xvalues[i])*(Xvalues[i])*massHydrogen*omega)/(2.0*hcross2));
                                                       PsyX [i] = sqrt(2.0*massHydrogen2*omega/hcross2)*(Xvalues[i])*exp((-1.0*(Xvalues[i])*(Xvalues[i])*massHydrogen*omega)/(2.0*hcross2));
                                                   }

                                                   }





void PotentialHamiltonianX(double (&Vxy)[gridsize],double Xvalues[gridsize],int gsize) //FILLS POTENTIAL HAM OPRTR
{

for (int i=0 ; i< gsize ; i++)
   {
    complex<double> temp (i,0);
       //  Vxy[i] =0.0;
        Vxy[i] = 0.5*massHydrogen2*Xvalues[i]*Xvalues[i]*omega*omega;
    }

}

void KineticHamiltonianFiller (double (&K)[gridsize][gridsize])
{
 for (int i=0;i<gridsize;i++)     //KINETIC ENERGY HAMILTONIAN
     {
       for (int j=0;j<gridsize;j++)
                {
                  if (j==i) {(KineticHamiltonianOperator [i][j]) = (KineticHamiltonianConst* pow(-1,i-j)*pi2*pi2)/(3.0)+ PotentialHamiltonianOperator[i];} //V(x) = +1/2m*w^2*x^2 diagnol elements
                   else     { (KineticHamiltonianOperator [i][j]) = (KineticHamiltonianConst* pow(-1,i-j)*(2/pow(i-j,2))); }


                  // if (j==i) {K [i][j] = one; } //V(x) = +1/2m*w^2*x^2 diagnol elements
                  // else     {K [i][j] = zero; }
                  }
         }


}

void Kfunction(complex<double>(&Ki)[gridsize],complex<double> PsyX[],complex<double> kprev[],double Kdenomfactor,double dt)
    {

complex<double> KofPsyPlusKi[gridsize];
complex<double> VofPsyPlusKi[gridsize];
//static complex<double> FPsyPlusKiTimesDt[gridsize];
for (int i=0;i<gridsize;i++)
       {
         complex<double> temp(0,0);
         for (int j=0;j<gridsize;j++)

              {temp = temp + KineticHamiltonianOperator[i][j]*(PsyX[j] + (kprev[j]/Kdenomfactor));}

                 KofPsyPlusKi[i] = temp;
                 VofPsyPlusKi[i] = PotentialHamiltonianOperator[i]*(PsyX[i] + (kprev[i]/Kdenomfactor)) ;
                // VofPsyPlusKi[i] = PotentialHamiltonianOperator[i];
               //  Ki[i] = ((KofPsyPlusKi[i] +  VofPsyPlusKi[i])/(iota*hcross))*dt;
              Ki[i] = ((KofPsyPlusKi[i] )/(iota*hcross2))*dt;
                  //cout<< Ki[i] <<endl;
                 }

          }
void XvaluesFillfunctionOfGridsize (double(&K)[gridsize])
{
   //if(gridsize%2!=0){
                          double temp= -1*stretch;
                          for (int i =0;i< gridsize;i++)
                                {


                                   K[i] = temp;
                                   temp = temp + dx;

                                   }

                    }

//WIGGLY NATURE is caused by the choice of dx!!
