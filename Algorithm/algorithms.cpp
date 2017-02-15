#include "Algorithm/algorithms.h"

Algorithms::Algorithms()
{

}

double Algorithms::RSS(QVector<QPair<double, double>> dataVector, const std::function<double (double)> func) {
    double rss = 0;
    for (QPair<double, double> p : dataVector) {
        double diff = p.second - func(p.first);
        rss += diff * diff;
    }
    return rss;
}
