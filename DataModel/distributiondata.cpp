#include "DataModel/distributiondata.h"
#include <math.h>
DistributionData::DistributionData()
{
    pointsChanged = true;
}

void DistributionData::setPoints(QVector<QPair<double, double> > points) {
    this->points = points;
    pointsChanged = true;
}

void DistributionData::addPoint(QPair<double, double> point) {
    points.append(point);
    pointsChanged = true;
}

void DistributionData::addPoint(double x, double y) {
    auto pair = qMakePair(x, y);
    points.append(pair);
    pointsChanged = true;
}

QVector<QPair<double, double>> DistributionData::getPoints() {
    return points;
}

QVector<double> DistributionData::getDistributionParameters(int from, int to) {
    QVector<QPair<double, double>> points = getStepRelativePoints();

    double fSumm = 0.0;
    double xfSumm = 0.0;
    int m = points.size();
    for (int i = from; i < to; ++i) {
        auto p = points.at(i);
        fSumm += p.second;
        xfSumm += p.first * p.second;
    }
    double mean = xfSumm / fSumm;

    double variance = 0.0;
    double assymetry = 0.0;
    double excess = 0.0;
    for (int i = from; i < to; ++i) {
        auto p = points.at(i);
        double diff = p.first - mean;
        double square = diff * diff;
        variance += square * p.second;
        assymetry += square * diff * p.second;
        excess += square * square * p.second;
    }


    variance /= fSumm;
    double sigma = sqrt(variance);
    assymetry = assymetry / pow(sigma, 3) / fSumm;
    excess = excess/variance/variance/fSumm - 3;

    return QVector<double>({mean, sigma, assymetry, excess});
}

DisVector DistributionData::getRelativePoints() {
    if (pointsChanged) {
        double summ = 0;
        for (QPair<double, double> &p : points) {
            summ += p.second;
        }

        for (QPair<double, double> &p : points) {
            relativePoints.append(qMakePair<double, double>(p.first, p.second/summ));
        }
        pointsChanged = false;
    }

    return relativePoints;
}

DisVector DistributionData::getStepRelativePoints() {
    DisVector stepRelativePoints;

    for (auto &p : getRelativePoints()) {
        stepRelativePoints.append(qMakePair<double, double>(p.first, p.second / 0.5));
    }

    return stepRelativePoints;
}
