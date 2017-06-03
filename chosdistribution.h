#ifndef CHOSDISTRIBUTION_H
#define CHOSDISTRIBUTION_H

#include <complex>
#include <vector>
#include "DataModel/distributiondata.h"
#include "utility"
using namespace std;

class ChosDistribution
{
public:
    ChosDistribution();
    //void setDistribution(DistributionData *distribution);
    void setInitialParams(vector<double> initialParams);

    void setPoints(PointsVector points);

    static double RSS(PointsVector dataVector, double mean, double sig, double as, double ex);

    double currentRss();

    static double value(double x, double m, double a, double beta, double ny);

    vector<double> gradDescent(size_t shakeCount);
    vector<double> gradLinDens(vector<double> params, double dens, double a, double b);

   // pair<vector<double>, vector<vector<double>> > derevetives();

    static double valueWithDistrParams(double x, double mean, double sig, double as, double ex);

    vector<double>descentProgress;


private:
    static double A(double beta, double ny);
    static complex<double> cbeta(complex<double> z1, complex<double> z2);
    static complex<double> gamma(complex<double> z);

    static vector<double> functionGradient(double x, double mean, double sig, double as, double ex);
    vector<double> gradLin(vector<pair<double, double>> points, vector<double> params);
    vector<double> gradQuadr(vector<pair<double, double>> points, vector<double> params);

    pair<double, vector<double>>  gradDescent(vector<double> initParams);
    pair<double, vector<double>>  gradDescentDensity(vector<double> initParams, double densityShift, double a, double b);
    //DistributionData *distribution;

    PointsVector points;

    vector<double> initialParams;
    vector<double> currParams;
    double currRss;

    vector<double> shakeParams(double m, double s);
};

#endif // CHOSDISTRIBUTION_H
