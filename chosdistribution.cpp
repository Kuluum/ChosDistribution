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


void ChosDistribution::gradDescent() {
    PROFILE_BLOCK("grad descent");
    double rssMin = 0.0001;
    bool descentSuccess = false;
    double previousRss = 1e+6;
//    double prevPreviousRss = 1e+6;

    vector<pair<double, vector<double>> >shakeVector;

    double learnRate = 1;

    vector<double> newParams(currParams);
    int step = 0;
    while(true)
    {
        step++;
        descentSuccess = false;
        vector<double> grad = gradLin();

        // Seve small grad elements.
        for (int i = 0; i < 4; i ++) {
            if (fabs(grad[i]) > 0.001) {
                newParams[i] -= learnRate * grad[i];
                descentSuccess = true;
            }
        }

        double r = RSS(distribution->getStepRelativePoints(), newParams[0], newParams[1], newParams[2], newParams[3]);

        qDebug() << step << "learn rate = " << learnRate << " rss = " << r;
        qDebug() << "params = " << newParams << "\n";

        // Bad descent step.
        if (isnan(r) || r >= previousRss) {
            // Bad descent fine.
            learnRate -= 0.5;

            // There were too much bad steps,assume that we get curent local minimum. Lets make shake.
            if (learnRate <= 0) {
                learnRate = 1;
                shakeVector.push_back(make_pair(previousRss, currParams));
                newParams = shakeParams();
                copy(newParams.begin(), newParams.end(), currParams.begin());
                previousRss = 1e+6;
            }

            copy(currParams.begin(), currParams.end(), newParams.begin());
            qDebug() << "bad descent";
        }
        // Good descent step.
        else {
            // Good descent promotion.
            learnRate += 0.1;

            previousRss = r;
            copy(newParams.begin(), newParams.end(), currParams.begin());
        }

        // Get target.
        if(r <= rssMin) {
            break;
        }

        // Last descent was too small, we get local minimum. Lets make shake to continue searching.
        if (!descentSuccess) {
            learnRate = 1;
            shakeVector.push_back(make_pair(previousRss, currParams));
            newParams = shakeParams();
            copy(newParams.begin(), newParams.end(), currParams.begin());
            previousRss = 1e+6;
        }

        if (shakeVector.size() >= 50) {
            sort(shakeVector.begin(), shakeVector.end(),
                [](const pair<double, vector<double>> &a, const pair<double, vector<double>> & b) -> bool
            {
                return a.first > b.first;
            });
            break;
        }
    }

}

double fRand(double fMin, double fMax)
{
    srand(time(0));
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

vector<double> ChosDistribution::shakeParams()
{

    double minMean = -10;
    double maxMean = 10;

    double minSig = 0.4;
    double maxSig = 10;

    double minAs = -2.5;
    double maxAs = 2.5;

    double maxEx = 20;
    srand(time(0));
    double mean = fRand(minMean, maxMean);
    double sig = fRand(minSig, maxSig);
    double as = fRand(minAs, maxAs);
    double ex = fRand (1.5 * as * as, maxEx);

    return vector<double>({mean, sig, as, ex});
}

vector<double> ChosDistribution::gradLin()
{
    vector<double>resGrad(4, 0.0);

    for (pair<double, double> p : distribution->getStepRelativePoints())
    {
        vector<double> funcGrad = functionGradient(p.first, currParams[0], currParams[1], currParams[2], currParams[3]);
        double predictionError = valueWithDistrParams(p.first, currParams[0], currParams[1], currParams[2], currParams[3]) - p.second;
        for (int i = 0; i < 4; i++) {
            resGrad[i] += 2 * predictionError * funcGrad[i];
        }
    }

    return resGrad;
}


vector<double> ChosDistribution::gradQuadr()
{
    vector<double> grad(4, 0.0);
    vector<vector<double>> hess(4, vector<double>(4, 0.0));

    for (pair<double, double> p : distribution->getStepRelativePoints())
    {
        vector<double> funcGrad = functionGradient(p.first, currParams[0], currParams[1], currParams[2], currParams[3]);
        double predictionError = valueWithDistrParams(p.first, currParams[0], currParams[1], currParams[2], currParams[3]) - p.second;
        for (int i = 0; i < 4; i++) {
            grad[i] += 2 * predictionError * funcGrad[i];
            for (int j = 0; j < 4; ++j) {
                hess[i][j] += 2 * funcGrad[i] * funcGrad[j];
            }
        }
    }

    vector<vector<double>> invHess = Algorithms::mesaInvertMatrix(hess);

    vector<double> resGrad = Algorithms::hessXgrad(invHess, grad);

    return resGrad;
}

double ChosDistribution::RSS(DisVector dataVector,double mean, double sig, double as, double ex) {
    if (distribution == nullptr) {
        return 1e+13;
    }

    double rss = Algorithms::RSS(dataVector, [mean, sig, as, ex](double x) {
        return valueWithDistrParams(x, mean, sig, as, ex);
    });

    return rss;
}


vector<double> ChosDistribution::functionGradient(double x, double mean, double sig, double as, double ex) {
    double h = 0.0001;

    double dMean = (valueWithDistrParams(x, mean+h, sig, as, ex) - valueWithDistrParams(x, mean-h, sig, as, ex)) / (2*h);
    double dSig = (valueWithDistrParams(x, mean, sig+h, as, ex) - valueWithDistrParams(x, mean, sig-h, as, ex)) / (2*h);
    double dAs = (valueWithDistrParams(x, mean, sig, as+h, ex) - valueWithDistrParams(x, mean, sig, as-h, ex)) / (2*h);
    double dEx = (valueWithDistrParams(x, mean, sig, as, ex+h) - valueWithDistrParams(x, mean, sig, as, ex-h)) / (2*h);

    return vector<double>({dMean, dSig, dAs, dEx});
}

