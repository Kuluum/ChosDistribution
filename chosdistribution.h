#ifndef CHOSDISTRIBUTION_H
#define CHOSDISTRIBUTION_H

#include <complex>
using namespace std;

class ChosDistribution
{
public:
    static double value(double x, double m, double a, double beta, double ny);

private:
    static double A(double beta, double ny);
    static complex<double> cgamma(complex<double> z,int OPT);
    static complex<double> cbeta(complex<double> z1, complex<double> z2);
};

#endif // CHOSDISTRIBUTION_H
