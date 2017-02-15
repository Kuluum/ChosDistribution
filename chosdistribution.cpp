#include "chosdistribution.h"
#include <math.h>
#include <complex>
#include <Algorithm/algorithms.h>

#include "QDebug"

ChosDistribution::ChosDistribution() {
}

void ChosDistribution::setDistribution(DistributionData *distribution) {
    this->distribution = distribution;
}

void ChosDistribution::setInitialParams(vector<double> initialParams) {
    this->initialParams = initialParams;
    currParams = initialParams;
}

//https://en.wikipedia.org/wiki/Lanczos_approximation
complex<double> ChosDistribution::gamma(complex<double> z)
{
    double p[8] = {676.5203681218851,
    -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7};

    double pi = M_PI;
    if ( real(z)<0.5 ) {
        return pi / (sin(pi*z)*gamma(1.0-z));
    }

    z -= 1.0;
    complex<double>x = 0.99999999999980993;

    for (int i=0; i<8; i++) {
        x += p[i]/(z+complex<double>(i+1,0));
    }

    complex<double>t = z + (8 - 0.5);
    return sqrt(2*pi) * pow(t,z+0.5) * exp(-t) * x;
}

complex<double> ChosDistribution::cbeta(complex<double> z1, complex<double> z2) {
    return gamma(z1) * gamma(z2) / gamma(z1+z2);

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

double ChosDistribution::valueWithDistrParams(double x, double mean, double sig, double as, double ex) {
    double a = mean;
    double m = 2 / (ex - as*as);
    double ny = m*sig*as/2;
    double beta = sqrt(m*sig*sig - ny*ny);

    return ChosDistribution::value(x, m, a, beta, ny);
}

vector<double> ChosDistribution::iterate()
{
    for (int i = 0; i < 10; i++) {
        vector<double> grad = dGrad();
        qDebug() << "grad = " << grad;
        for (int i = 0; i < 4; i ++) {
            if (fabs(grad[i]) > 1.0e-6) {
                currParams[i] += grad[i];
            }
        }
        qDebug() << "params = " << currParams;
        qDebug() << "RSS = " << RSS(currParams[0], currParams[1], currParams[2], currParams[3]);
        qDebug() << "\n";
    }

    return currParams;
}

double ChosDistribution::RSS(double mean, double sig, double as, double ex) {
    if (distribution == nullptr) {
        return 1e+13;
    }

//    double rss = 0;
//    for (QPair<double, double> p : distribution->getStepRelativePoints()) {
//        double diff = p.second - valueWithDistrParams(p.first, mean, sig, as, ex);
//        rss += diff * diff;
//    }
//    return rss;
    double rss = Algorithms::RSS(distribution->getStepRelativePoints(), [mean, sig, as, ex](double x) {
        return valueWithDistrParams(x, mean, sig, as, ex);
    });

    return rss;
}

vector<double> ChosDistribution::dGrad()
{
    vector<double>resGrad({0.0, 0.0, 0.0, 0.0});

    for (QPair<double, double> p : distribution->getStepRelativePoints()) {
        vector<double> grad = gradDistr(p.first, currParams[0], currParams[1], currParams[2], currParams[3]);
        double diff = valueWithDistrParams(p.first, currParams[0], currParams[1], currParams[2], currParams[3]) - p.second;
        for (int i = 0; i < 4; i++) {
            double currG = grad[i];
            resGrad[i] += diff * currG;
        }
    }

    for (int i = 0; i < 4; i++) {
        resGrad[i] *= 2;
    }
    return resGrad;
}

std::vector<double> ChosDistribution::gradDistr(double x, double mean, double sig, double as, double ex)
{
    double h = 0.0001;

    //double derValX = (valueWithDistrParams(x+h, mean, sig, as, ex) - valueWithDistrParams(x-h, mean, sig, as, ex)) / (2*h);
    double derValMean = (RSS(mean+h, sig, as, ex) - RSS(mean-h, sig, as, ex)) / (2*h);
    double derValSig = (RSS(mean, sig+h, as, ex) - RSS(mean, sig-h, as, ex)) / (2*h);
    double derValAs = (RSS(mean, sig, as+h, ex) - RSS(mean, sig, as-h, ex)) / (2*h);
    double derValEx = (RSS(mean, sig, as, ex+h) - RSS(mean, sig, as, ex-h)) / (2*h);

    return std::vector<double>({derValMean, derValSig, derValAs, derValEx});
}

vector<vector<double>> ChosDistribution::dHess()
{

    vector<vector<double>> resHess = {{0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0, 0.0}};


    for (QPair<double, double> p : distribution->getStepRelativePoints()) {
         vector<double> grad = gradDistr(p.first, currParams[0], currParams[1], currParams[2], currParams[3]);
         double diff = valueWithDistrParams(p.first, currParams[0], currParams[1], currParams[2], currParams[3]) - p.second;
         vector<vector<double>> hess = hessianDistr(p.first, currParams[0], currParams[1], currParams[2], currParams[3]);
         for (int i = 0; i < 4; i ++) {
             for (int j = 0; j < 4; j ++) {
                 //resHes[i][j] += 2 * (grad[i] * grad[j] - diff*hess[i][j]);
             }
         }
    }

}

vector<vector<double>> ChosDistribution::hessianDistr(double x, double mean, double sig, double as, double ex)
{
    vector<double> params = {mean, sig, as, ex};

    double h = 0.0001;
    double arg4 = valueWithDistrParams(x, params[0], params[1], params[2], params[3]);

    vector<vector<double>> resHessian;
    for (int i = 0; i < 4; i++) {
        vector<double> jVect;
        for (int j = 0; j < 4; j++) {
            vector<double> dParams = params;
            dParams[i] += h;
            dParams[j] += h;
            double arg1 = valueWithDistrParams(x, dParams[0], dParams[1], dParams[2], dParams[3]);
            dParams[i] -= h;
            double arg2 = valueWithDistrParams(x, dParams[0], dParams[1], dParams[2], dParams[3]);
            dParams[j] -= h;
            dParams[i] += h;
            double arg3 = valueWithDistrParams(x, dParams[0], dParams[1], dParams[2], dParams[3]);

            jVect.push_back((arg1 - arg2 - arg3 + arg4)/h/h);
        }
        resHessian.push_back(jVect);
    }
    return resHessian;//std::vector<double>({derValMean, derValSig, derValAs, derValEx});
}
