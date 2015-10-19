//   This is a simple program which tests the trapezoidal rule, Simpsons' rule,
//   and Gaussian quadrature using Legendre and Laguerre polynomials
//   It integrates the simple function exp(-2*2(r1+r2))* 1 / abs(r1-r2) for the interval
//   x \in [0,infty). The exact result is 1. For Legendre based quadrature a
//   tangent mapping is also used.

#include <QCoreApplication>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;


//     Here we define various functions called by the main program

double int_function2(double i,double j,double k,double l,double m,double n);
double int_function(double i,double j,double k,double l,double m,double n);
void gauss_laguerre(double *, double *, int, double);
double trapezoidal_rule ( double, double, int, double (*func)(double) );
double simpson ( double, double, int, double (*func)(double) );
void gauleg(double, double, double *, double *, int);
double gammln(double);
double pi = 3.1415926536;

//   Main function begins here
int main()
{
     int N;
     double a, b;
     cout << "Read in the number of integration points" << endl;
     cin >> N;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *x = new double [N];
     double *w = new double [N];
     // Gauss-Laguerre is old-fashioned translation of F77 --> C++
     // arrays start at 1 and end at n
     double *xgl = new double [N+1];
     double *wgl = new double [N+1];
     double *xphi = new double [N];
     double *xtheta = new double [N];
     double *wphi = new double [N];
     double *wtheta = new double [N];

     // These arrays are used for improved Gauss-Legendre, mapping of
     // x \in [-1,1] to x \in [0, infinity)
     double *r = new double [N];
     double *s = new double [N];

//   set up the mesh points and weights
     gauleg(a,b,x,w, N);

     //   set up the mesh points and weights and the power of x^alf
          double alf = 2.0;
          gauss_laguerre(xgl,wgl, N, alf);
          gauleg(0, 2*pi, xphi, wphi, N);
          gauleg(0, pi, xtheta, wtheta, N);

//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum. Here brute force gauleg
     double int_gauss = 0.;
     for ( int i = 0;  i < N; i++){
         for ( int j = 0;  j < N; j++){
         for ( int k = 0;  k < N; k++){
         for ( int l = 0;  l < N; l++){
         for ( int m = 0;  m < N; m++){
         for ( int n = 0;  n < N; n++){
        int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
         }}}}}
     }
     //   evaluate the integral with the Gauss-Laguerre method
     //   Note that we initialize the sum
          double int_gausslag = 0.;
          for ( int i = 1;  i <= N; i++){
              for ( int j = 0;  j <= N; j++){
              for ( int k = 0;  k < N; k++){
              for ( int l = 0;  l < N; l++){
              for ( int m = 0;  m < N; m++){
              for ( int n = 0;  n < N; n++){

             int_gausslag += wgl[i]*wgl[j]*wphi[k]*wphi[l]*wtheta[m]*wtheta[n]*int_function2(xgl[i],xgl[j],xphi[k],xphi[l],xtheta[m],xtheta[n])*(sin(xtheta[m])*(sin(xtheta[n])))/16/16/16;
              }}}}}

          }


      cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss << endl;
      cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15) << int_gausslag << endl;
      //cout << "Gaussian-Legendre improved quad = " << setw(20) << setprecision(15) << int_gausslegimproved << endl;
      delete [] x;
      delete [] w;
      delete [] wgl;
      delete [] xgl;
      delete [] xphi;
      delete [] xtheta;
      delete [] wphi;
      delete [] wtheta;

      return 0;
}  // end of main program
//  this function defines the function to integrate

double int_function(double i,double j,double k,double l,double m,double n){
    double alpha = 2.;
 // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(i*i+j*j+k*k);
    double exp2=-2*alpha*sqrt(l*l+m*m+n*n);
    double deno=sqrt(pow((i-l),2)+pow((j-m),2)+pow((k-n),2));
 // Cheating!!
   if(deno <pow(10.,-6.)) { return 0;}
   else return exp(exp1+exp2)/deno;
 } // end of function to evaluate

double int_function2(double r1,double r2, double theta1, double theta2, double phi1, double phi2){
    double beta =(cos(theta1)*cos(theta2))+(sin(theta1)*sin(theta2)*cos(phi1 - phi2));
    double root = sqrt(((r1*r1)/16)+((r2*r2)/16)-((2*r1*r2*beta)/16));
    if(root > pow(10.,-6.)){
        return 1/root;
    }else{
        return 0;
    }

}



       /*
       ** The function
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


/*     function to integrate a function func over the                   */
/*     interval [a,b] with input a, b, and the number of steps          */
/*     n.  it returns the sum as the variable trapez_sum                */
/*     the trapezoidal rule is used                                     */


double trapezoidal_rule(double a, double b, int n, double (*func)(double))
{

      double trapez_sum;
      double fa, fb, x, step;
      int    j;

      step=(b-a)/((double) n);
      fa=(*func)(a)/2. ;
      fb=(*func)(b)/2. ;
      trapez_sum=0.;

      for (j=1; j <= n-1; j++){
         x=j*step+a;
         trapez_sum+=(*func)(x);
      }

      trapez_sum=(trapez_sum+fb+fa)*step;
      return trapez_sum;

}  /* end trapezoidal_rule  */

/*     function to integrate a function func over the          */
/*     interval [a,b] with input a, b, and the number of steps */
/*     n.  it returns the sum as the variable simpson_sum      */
/*     simpson's method is used                                */


double simpson(double a, double b, int n, double (*func)(double))
{
      double simpson_sum;
      double fa, fb, x, step, fac;
      int    j;

      step = (b-a)/((double) n);
      fa=(*func)(a) ;
      fb=(*func)(b) ;
      simpson_sum=fa ;
      fac=2.;

      for (j=1; j <= n-1 ; j++){
           if ( fac == 2.){
               fac = 4.;
           }
           else{
               fac = 2.;
           }  /* end of if test */
          x=j*step+a;
          simpson_sum+=(*func)(x)*fac;
      }  /* end of for loop */

      simpson_sum=(simpson_sum+fb)*step/3.;
      return simpson_sum;

}  /*    end function simpson   */


void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
//#undef EPS
//#undef MAXIT
