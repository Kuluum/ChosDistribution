#include "chosdistribution.h"
#include <math.h>
#include <complex>

#include <complex>

complex<double> ChosDistribution::cgamma(complex<double> z,int OPT) {
    complex<double> g,z0,z1;
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na,t,x1,y1,sr,si;
    int i,j,k;

    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};

    x = real(z);
    y = imag(z);
    if (x > 171) return complex<double>(1e308,0);
    if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return complex<double>(1e308,0);
    else if (x < 0.0) {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
    }
    x0 = x;
    if (x <= 7.0) {
        na = (int)(7.0-x);
        x0 = x+na;
    }
    q1 = sqrt(x0*x0+y*y);
    th = atan(y/x0);
    gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
    gi = th*(x0-0.5)+y*log(q1)-y;
    for (k=0;k<10;k++){
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
    }
    if (x <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) {
            gr1 += (0.5*log((x+j)*(x+j)+y*y));
            gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (x1 <= 0.0) {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
        th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
    }
    if (OPT == 0) {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
    }
    g = complex<double>(gr,gi);
    return g;
}

complex<double> ChosDistribution::cbeta(complex<double> z1, complex<double> z2) {
    return cgamma(z1, 0)*cgamma(z2, 0)/cgamma(z1+z2, 0);
}

double ChosDistribution::A(double beta, double ny) {
    double squareDiff = beta*beta - ny*ny;
    int sig = (squareDiff >= 0) ? 1 : -1;

    if(squareDiff > 0) {
        return exp(atan(2*beta*ny/squareDiff));
    }
    else if (squareDiff < 0) {
        return exp(atan(2*beta*ny/squareDiff+sig*M_PI));
    }
    else {
        return exp(M_PI_2*sig);
    }
}

double ChosDistribution::value(double x, double m, double a, double beta, double ny) {
   complex<double> z1 = complex<double>(m/2, -1*m*(x+ny-a)/2/beta);
   complex<double> z2 = complex<double>(m/2, m*(x+ny-a)/2/beta);
   complex<double> betaFunc = cbeta(z1, z2);

   double betaRe = betaFunc.real();
   double fA = pow(A(beta, ny), m*(x+ny-a)/2/beta);

   double value = pow(2, m-2) * m * pow(beta, m-1) * fA * betaRe / M_PI / pow(beta*beta+ny*ny, m/2);

   return value;
}
