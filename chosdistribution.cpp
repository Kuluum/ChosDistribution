#include "chosdistribution.h"
#include <math.h>
#include <complex>
#include <Algorithm/algorithms.h>
#include "QDebug"

ChosDistribution::ChosDistribution() {
}

void ChosDistribution::setInitialParams(vector<double> initialParams) {
    this->initialParams = initialParams;
    currParams = initialParams;
}

void ChosDistribution::setPoints(PointsVector points) {
    this->points = points;
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
    int sign = (ny >= 0) ? 1 : -1;

    if(squareDiff > 0) {
        return exp(atan(2*beta*ny/squareDiff));
    }
    else if (squareDiff < 0) {
        return exp(atan(2*beta*ny/squareDiff) + sign*M_PI);
    }
    else {
        return exp(M_PI_2 * sign);
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
    if (sig < 0.1) {
        return nan("");
    }

    double a = mean;
    double m = 2 / (ex - as*as);
    double ny = m*sig*as/2;
    double beta = sqrt(m*sig*sig - ny*ny);

    return ChosDistribution::value(x, m, a, beta, ny);
}


vector<double> ChosDistribution::gradDescent(size_t shakeCount)
{
//    PROFILE_BLOCK("grad descent");
    double bestRss = 100000;
    vector<double> bestParams;
    double rssMin = 0.005;

    vector<double> initParams = this->initialParams;
    if (shakeCount > 0) {
       for (int i = 0; i < shakeCount+1; i++)
       {
           auto descentRes = gradDescent(initParams);
           if (descentRes.first < bestRss) {
               bestRss = descentRes.first;
               bestParams = descentRes.second;
           }
//           if (bestRss <= rssMin)
//           {
//               break;
//           }
           initParams = shakeParams(initialParams[0], initialParams[1]);
       }
    }
    else
    {
        auto descentRes = gradDescent(initParams);
        bestParams = descentRes.second;
    }
    return bestParams;
}

pair<double, vector<double>> ChosDistribution::gradDescent(vector<double> initParams)
{
    bool descentSuccess = false;
    double rssMin = 0.001;
    double previousRss = 1e+6;
    double learnRate = 1;

    vector<double> currentParams(initParams);
    vector<double> bestParams(4, 0.0);

    qDebug() << "start gradient with init params " << initParams;
    double r = RSS(points, currentParams[0], currentParams[1], currentParams[2], currentParams[3]);
    if (isnan(r)) {
        currentParams[2] = 0;
        currentParams[3] = 0.1;
    }
    // Get the target.
    else if(r <= rssMin) {
        qDebug() << "Good match riched with rss " << r;
        bestParams = currentParams;
        qDebug() << "BEST PARAMS: " << bestParams;
        currRss = r;
        return make_pair(r, bestParams);
    }

    vector<double> newParams(currentParams);
    int step = 0;
    while(true)
    {
        step++;
        descentSuccess = false;
        vector<double> grad = gradLin(points, currentParams);

        // Seve small grad elements.
        for (int i = 0; i < 4; i ++) {
            if (fabs(grad[i]) > 0.001) {
                newParams[i] -= learnRate * grad[i];
                descentSuccess = true;
            }
        }

        double r = RSS(points, newParams[0], newParams[1], newParams[2], newParams[3]);
        descentProgress.push_back(r);
        qDebug() << step << ") learn rate = " << learnRate << " rss = " << r;
        qDebug() << "params = " << newParams;

        // Last descent was too small, we get local minimum. Lets make shake to continue searching.
        if (!descentSuccess) {
            qDebug() << "Descent too small. Local minimum riched.";
                bestParams = currentParams;
                currRss = previousRss;
                break;
        }

        // Get target.
        if(r <= rssMin) {
            qDebug() << "Good match riched.";
            bestParams = newParams;
            currRss = r;
            break;
        }
        // Bad descent step.
        else if (isnan(r) || r >= previousRss) {
            qDebug() << "bad descent";

            // Bad descent fine.
            learnRate -= 0.5;
            if (learnRate == 0) {
                learnRate = 0.05;
            }
            // There were too much bad steps, assume that we get curent local minimum. Lets make shake.
            if (learnRate < 0) {
                    bestParams = currentParams;
                    currRss = previousRss;
                    break;
            }
            else
            {
                newParams = currentParams;
            }
        }
        // Good descent step.
        else {
            // Good descent promotion.
            learnRate += 0.1;

            previousRss = r;
            currentParams = newParams;
        }

    }

    qDebug() << "BEST PARAMS: " << bestParams << " rss: " << currRss;
    currentParams = bestParams;
    return make_pair(currRss, bestParams);
}

pair<double, vector<double>> ChosDistribution::gradDescentDensity(vector<double> initParams, double densityShift, double a, double b)
{
    bool descentSuccess = false;
    double rssMin = 0.005;
    double previousRss = 1e+6;
    double learnRate = 1;

    vector<double> currentParams(initParams);
    vector<double> bestParams(4, 0.0);

    qDebug() << "start gradient with init params " << initParams;

    vector<double> newParams(currentParams);
    int step = 0;
    while(true)
    {
        step++;
        descentSuccess = false;
        vector<double> grad = gradLinDens(currentParams, densityShift, a, b);

        // Seve small grad elements.
        for (int i = 0; i < 4; i ++) {
            if (fabs(grad[i]) > 0.001) {
                newParams[i] -= learnRate * grad[i];
                descentSuccess = true;
            }
        }

        double integralValue = Algorithms::Integral(a, b, [&](double d)->double{
            double v = ChosDistribution::value(d, newParams[0], newParams[1], newParams[2], newParams[3]);
            return v;
        });
        //double predictionError = (integralValue + densityShift) - 1;

        double r = (integralValue + densityShift) - 1;

        //qDebug() << step << ") learn rate = " << learnRate << " rss = " << r;
        //qDebug() << "params = " << newParams;

        // Last descent was too small, we get local minimum. Lets make shake to continue searching.
        if (!descentSuccess) {
           // qDebug() << "Descent too small. Local minimum riched.";
                bestParams = currentParams;
                currRss = previousRss;
                break;
        }

        // Get target.
        if(r <= rssMin) {
            qDebug() << "Good match riched. No need shake.";
            bestParams = newParams;
            currRss = r;
            break;
        }
        // Bad descent step.
        else if (isnan(r) || r >= previousRss) {
          //  qDebug() << "bad descent";

            // Bad descent fine.
            learnRate -= 0.5;
            if (learnRate == 0) {
                learnRate = 0.05;
            }
            // There were too much bad steps, assume that we get curent local minimum. Lets make shake.
            if (learnRate < 0) {
                    bestParams = currentParams;
                    currRss = previousRss;
                    break;
            }
            else
            {
                newParams = currentParams;
            }
        }
        // Good descent step.
        else {
            // Good descent promotion.
            learnRate += 0.1;

            previousRss = r;
            currentParams = newParams;
        }

    }

   // qDebug() << "BEST PARAMS: " << bestParams << " rss: " << currRss;
    currentParams = bestParams;
    return make_pair(currRss, bestParams);
}

#include <random>

double fRand(double fMin, double fMax)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(fMin, fMax);
   // double a_random_double = -1;
   // for (int i = 0; i < 10; i ++) {
    double a_random_double = dis(gen);
   // qDebug() << "frand (" << fMin << "; " << fMax << ") = " << a_random_double;
   // }
    return a_random_double;
}

vector<double> ChosDistribution::shakeParams(double m, double s)
{
    double minSig = 0.4;
    double maxSig = 10;

    double minAs = -2.5;
    double maxAs = 2.5;

    double maxEx = 20;
    double mean = m * fRand(0.33, 3);
    double sig = s * fRand(0.33, 3);
    if (sig < minSig)
    {
        sig = minSig;
    }
    else if (sig > maxSig)
    {
        sig = maxSig;
    }

    double as = fRand(minAs, maxAs);
    double minEx = 1.5 * as * as;
    double ex = fRand (minEx, maxEx);

    auto v = vector<double>({mean, sig, as, ex});
    return v;
}

vector<double> ChosDistribution::gradLin(vector<pair<double, double>> points, vector<double> params)
{
    vector<double>resGrad(4, 0.0);

    for (pair<double, double> p : points)
    {
        vector<double> funcGrad = functionGradient(p.first, params[0], params[1], params[2], params[3]);
        double predictionError = valueWithDistrParams(p.first, params[0], params[1], params[2], params[3]) - p.second;
        for (int i = 0; i < 4; i++) {
            resGrad[i] += 2 * predictionError * funcGrad[i];
        }
    }

    return resGrad;
}

vector<double> ChosDistribution::gradLinDens(vector<double> params, double dens, double a, double b)
{

    vector<double>resGrad(4, 0.0);
    double step = points[1].first - points[0].first;
    for (auto p : points) {
    vector<double> funcGrad = functionGradient(p.first, params[0], params[1], params[2], params[3]);
    double integralValue = Algorithms::Integral(p.first-step/2, p.first+step/2, [&](double d)->double{
        double v = ChosDistribution::value(d, params[0], params[1], params[2], params[3]);
        return v;
    });
    double predictionError = (integralValue + dens) - 1;
    for (int i = 0; i < 4; i++) {
        resGrad[i] += 2 * predictionError * funcGrad[i];
    }
    }
    return resGrad;
}

vector<double> ChosDistribution::gradQuadr(vector<pair<double, double>> points, vector<double> params)
{
    vector<double> grad(4, 0.0);
    vector<vector<double>> hess(4, vector<double>(4, 0.0));

    for (pair<double, double> p : points)
    {
        vector<double> funcGrad = functionGradient(p.first, params[0], params[1], params[2], params[3]);
        double predictionError = valueWithDistrParams(p.first, params[0], params[1], params[2], params[3]) - p.second;
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

double ChosDistribution::RSS(PointsVector dataVector,double mean, double sig, double as, double ex) {

    double rss = Algorithms::RSS(dataVector, [mean, sig, as, ex](double x) {
        return valueWithDistrParams(x, mean, sig, as, ex);
    });

    return rss;
}

double ChosDistribution::currentRss() {
    return currRss;
}

vector<double> ChosDistribution::functionGradient(double x, double mean, double sig, double as, double ex) {
    double h = 0.0001;

    double dMean = (valueWithDistrParams(x, mean+h, sig, as, ex) - valueWithDistrParams(x, mean-h, sig, as, ex)) / (2*h);
    double dSig = (valueWithDistrParams(x, mean, sig+h, as, ex) - valueWithDistrParams(x, mean, sig-h, as, ex)) / (2*h);
    double dAs = (valueWithDistrParams(x, mean, sig, as+h, ex) - valueWithDistrParams(x, mean, sig, as-h, ex)) / (2*h);
    double dEx = (valueWithDistrParams(x, mean, sig, as, ex+h) - valueWithDistrParams(x, mean, sig, as, ex-h)) / (2*h);

    return vector<double>({dMean, dSig, dAs, dEx});
}

