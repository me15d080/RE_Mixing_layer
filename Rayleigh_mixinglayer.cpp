/****************************************************************/
/* AS5430 COURSE-PROJECT                                        */
/* LINEAR INVISCID INSTABILITY STUDY OF TANH VELOCITY PROFILE   */    
/* SUBMITTED BY: VIKAS DWIVEDI (ME15DO80)                       */
/****************************************************************/

# include <iostream>
# include <cmath>
# include <complex>
# include <cassert>
# include <fstream>
# include <iomanip>
# include <cstdlib>
using namespace std;

double** rk2_riccati (double**, double*, double*, double, double,double,int,int);

int main()
{

 char const* ch1="eig_fcn"; char const* ch2="disp_data"; 

//-->define iota

 complex<double> ii,iota; ii=-1;
 iota = sqrt(ii);

//------declaration of variables--------//

 int row_PHY=100  ; int col_PHY=2;
 int row_phy=701  ; int col_phy=2;

 double*  k_range           ; k_range = new double [9];
 double*  ci_range          ; ci_range= new double [9];
 double*  Z_span            ; Z_span  = new double [row_PHY];  
 double*  y_span            ; y_span  = new double [row_phy];
 double*  PHY_init          ; PHY_init= new double [col_PHY];
 double** PHY = 0           ; PHY     = new double*[row_PHY];
 complex <double>*  phy_init; phy_init= new complex<double> [col_phy];
 complex <double>** phy = 0 ; phy     = new complex<double>*[row_phy];

 
  for(int i=0;i<9;i++)
    {
     k_range  [i] = 0.11*(1.0+double(i)); // k = [0.11,0.22,...0.99], eigen functions are known for k=0,1 
     ci_range [i] = 0.0;
    }
   
  for (int i=0;i<row_PHY;i++)
    {
      Z_span[i]=-0.99+(double(i)/(row_PHY-1))*0.99; // Z = [-1,0] because of symmetry. 1st paragraph Page 547, Michake
    }
  // y_span depends upon the plots we want to replicate
  for (int i=0;i<row_phy;i++)
    {
        y_span[i]= 0.0+(double(i)/(row_phy-1))*7.0;// for eigenfunction plots (Round-1)
      //y_span[i]= 0.0+(double(i)/(row_phy-1))*(-5.0);// for streamfcn plots (Round-1)
      // y_span[i]= -5.0+(double(i)/(row_phy-1))*(7.2); // for streamfcn plots (Round-2)
    }

  // PHY(Z)'s 1st & 2nd column contain real & imaginary part of PHY respectively
  for (int i = 0; i < row_PHY; i++)
    {
      PHY[i] = new double[col_PHY];

      for (int j = 0; j < col_PHY; j++)
      {                  
        PHY[i][j] =0.0;
      }
    }
  // phy(y)'s 1st & 2nd column contain phy (eig_fcn) & dphydy (first deriv) respectively
  for (int i = 0; i < row_phy; i++)
    {
      phy[i] = new complex<double>[col_phy];

      for (int j = 0; j < col_phy; j++)
      {                  
        phy[i][j] =0.0;
      }
    }

  
 double dZ = Z_span[1]-Z_span[0];
 double dy = y_span[1]-y_span[0];
 double k,error;
 double ci_left, ci_right, ci_middle, ci_middle_updated, ci_correct;
 double PHYr_at_0_left,PHYr_at_0_right,PHYr_at_0_middle,PHYi_at_0;
 complex<double> num,deno;

/*-------------------declaration done----------------------------*/


/*---------------wavenumber loop starts here---------------------*/

for (int count=0;count<9;count++)
//for (int count=3;count<4;count++)
{
   
   k = k_range[count];
 
   PHY_init[0] = k;PHY_init[1] = 0.0; // See eq (19) Page 546, Michalke

   ci_left = 1e-10; ci_right= 0.5;    // make guess 

/*--part-1: Solve EVP (Riccati eqn) using bisection method--*/

// In bisection, we aim to find f(x_root)=0 when x_left < x_root < x_right
// Here, f is PHYr(0) and x is ci. See eq (21) Page 547, Michalke

  PHY = rk2_riccati (PHY, PHY_init, Z_span, k, ci_left, dZ, row_PHY, col_PHY); //find f(ci_left) by integrating 
  PHYr_at_0_left= PHY[row_PHY-1][0];// f (ci_left)
  
  PHY = rk2_riccati (PHY, PHY_init, Z_span, k, ci_right, dZ, row_PHY, col_PHY);//find f(ci_right) by integrating
  PHYr_at_0_right= PHY[row_PHY-1][0];// f (ci_right)

  assert(PHYr_at_0_left*PHYr_at_0_right < 0); // make sure ci_correct is bracketed

  ci_middle = 0.5*(ci_left+ci_right);// bisect the bracket
  error = 100;

  while(error > 0.5)  // error within 0.5%
  {
    
    PHY = rk2_riccati (PHY, PHY_init, Z_span, k, ci_middle, dZ, row_PHY, col_PHY);
    PHYr_at_0_middle= PHY[row_PHY-1][0];
    PHYi_at_0 = PHY[row_PHY-1][1];// IC for Rayleigh equation See eq 26 Page 548, Michalke

    if (PHYr_at_0_middle*PHYr_at_0_right<0)
    {
      ci_left = ci_middle;
    }
    else
    {
      ci_right = ci_middle;
    }
    ci_middle_updated = 0.5*(ci_left+ci_right);

    error=100.0*fabs((ci_middle_updated-ci_middle)/ci_middle_updated);

    ci_middle = ci_middle_updated;
  }

  ci_correct = ci_middle; 
  ci_range [count] = ci_correct; // c(k) is known
  
/*--------part-2: integrate Rayleigh's equation (complex)---------------*/

// most amplified wave:
// y Re(phi[0][0])  Im(phi[0][0])  Re(phi[0][1])  Im(phi[0][1])
//-5 0.16464931    -0.070682136    0.084842613   -0.060205098

//-> (Round 1)
phy_init[0] =1.0 ;
phy_init[1] = iota*(PHYi_at_0); 

//-> (Round 2)
// phy_init[0] = 0.16464931  -0.070682136*iota;
// phy_init[1] = 0.084842613 -0.060205098*iota; 

for (int j=0;j<col_phy;j++ )
{
  phy[0][j]=phy_init[j];
}


/*--intermediate variables--*/
complex<double>* dphydy;
complex<double>* phy_guess;
complex<double>* dphydy_guess;

dphydy       = new complex<double>[col_phy];
phy_guess    = new complex<double>[col_phy];
dphydy_guess = new complex<double>[col_phy];

/*--integrate using Heun's method (RK2)--*/
 
for (int i=0;i<row_phy-1;i++)
{

//use Eq.3 [Rayleigh equation]

  dphydy[0]= phy[i][1];
         
  num = tanh(y_span[i])/(cosh(y_span[i])*cosh(y_span[i]));
  deno = 0.5*tanh(y_span[i])-ci_correct*iota;

  dphydy[1]= (k*k-(num/deno))*phy[i][0] ;// See eq 3, pg 544, Michalke

//predictor step1: guess the value of y(i+1) 

  for(int j=0;j<col_phy;j++)
  {
   phy_guess[j]=phy[i][j]+dphydy[j]*dy; 
  }

//predictor step2: guess the value of dphydx(i+1) using Eq.3 [Rayleigh equation]

  dphydy_guess[0]= phy_guess[1];

  num = tanh(y_span[i+1])/(cosh(y_span[i+1])*cosh(y_span[i+1]));
  deno = 0.5*tanh(y_span[i+1])-ci_correct*iota;

  dphydy_guess[1]= (k*k-(num/deno))*phy_guess[0];

//corrector step: average the slope at two ends

 for( int j=0;j<col_phy;j++)
 {
  dphydy[j]=(dphydy_guess[j]+dphydy[j])*0.5;
 }

//calculate y(i+1) with average slope
 for( int j=0;j<col_phy;j++)
 {
  phy[i+1][j]=phy[i][j]+dphydy[j]*dy;
 }

}

/*--clean up intermediate variables--*/
delete [] dphydy;
delete [] phy_guess;
delete [] dphydy_guess;

//--> create wavenumber folder
 char command[50];
 sprintf(command,"mkdir /home/vikas/Desktop/AS_5430_project_ME15D080/output/%d",count+1);
 system (command);

//--> write eigen functions

 char fname1[50];
 sprintf(fname1, "/home/vikas/Desktop/AS_5430_project_ME15D080/output/%d/%s",count+1,ch1);
 ofstream fout1;
 fout1.open(fname1);

 for(int i = 0; i < row_phy; i++)
    {
     fout1.precision(8); 
     fout1 << y_span[i] << " "<<real(phy[i][0]) << " " << imag(phy[i][0]) 
                        << " "<<real(phy[i][1]) << " " << imag(phy[i][1]) << "\n";
    }

 fout1.close();

}
/*---------- wavenumber-loop ends here-------------*/

//--> write dispersion data
 char fname2[50];
 sprintf(fname2, "/home/vikas/Desktop/AS_5430_project_ME15D080/output/%s",ch2);
 ofstream fout2;
 fout2.open(fname2);

 for(int i = 0; i < 9; i++)
    {
     fout2.precision(8); 
     fout2 << k_range[i] << " " << k_range[i]*ci_range[i] <<" " << ci_range[i] << "\n";            
    }

fout2.close();

/*----clean up-------*/

 delete [] k_range ;
 delete [] ci_range; 
 delete [] Z_span  ;  
 delete [] y_span  ; 
 delete [] PHY_init; 
 delete [] phy_init;

 for (int  i = 0; i < row_PHY; i++)
    {
     delete [] PHY[i];
    }
 delete [] PHY;
 for (int  i = 0; i < row_phy; i++)
    {
     delete [] phy[i];
     
    }
 delete [] phy;
 


/*--end of program--*/
return 0;
}

/*-----------------------------------*/
/*---- function description----------*/
/*-----------------------------------*/


double** rk2_riccati (double** y, double* y_init, double* x_span, double k, double ci,double dx,int row, int col)
{

/*--specify initial condition--*/

for (int j=0;j<col;j++ )
{
y[0][j]=y_init[j];
}

/*--intermediate variables--*/
double* dydx;
double* y_guess;
double* dydx_guess;

dydx       = new double[col];
y_guess    = new double[col];
dydx_guess = new double[col];

/*--integrate using Heun's method (RK2)--*/
 
for (int i=0;i<row-1;i++)
{
//See eq 18 [Riccati equation]  Pg 546, Michalke

  dydx[0]= ((k*k + y[i][1]*y[i][1]-y[i][0]*y[i][0])/(1-x_span[i]*x_span[i]))
           -((2*x_span[i]*x_span[i])/(x_span[i]*x_span[i]+4*ci*ci));

  dydx[1]= -((2*y[i][0]*y[i][1])/(1-x_span[i]*x_span[i]))
             -4*ci*x_span[i]/(x_span[i]*x_span[i]+4*ci*ci);

//predictor step1: guess the value of y(i+1) 
  for(int j=0;j<col;j++)
  {
   y_guess[j]=y[i][j]+dydx[j]*dx; 
  }

//predictor step2: guess the value of dydx(i+1) using Eq.18 [Riccati equation]

  dydx_guess[0]= ((k*k + y_guess[1]*y_guess[1]-y_guess[0]*y_guess[0])/(1-x_span[i+1]*x_span[i+1]))
                -((2*x_span[i+1]*x_span[i+1])/(x_span[i+1]*x_span[i+1]+4*ci*ci));
  dydx_guess[1]=-((2*y_guess[0]*y_guess[1])/(1-x_span[i+1]*x_span[i]))
                -4*ci*x_span[i+1]/(x_span[i+1]*x_span[i+1]+4*ci*ci);

//corrector step: average the slopes at two ends

 for( int j=0;j<col;j++)
 {
  dydx[j]=(dydx_guess[j]+dydx[j])*0.5;
 }

//calculate y(i+1) with average slope
 for( int j=0;j<col;j++)
 {
  y[i+1][j]=y[i][j]+dydx[j]*dx;
 }

}


/*--clean up intermediate variables--*/
delete [] dydx;
delete [] y_guess;
delete [] dydx_guess;


return y;
}
