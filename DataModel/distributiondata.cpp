#include "DataModel/distributiondata.h"
#include <math.h>

DistributionData::DistributionData()
{
    pointsChanged = true;
}

void DistributionData::setPoints(DisVector points) {
    this->points = points;
    pointsChanged = true;
}

void DistributionData::addPoint(pair<double, double> point) {
    points.push_back(point);
    pointsChanged = true;
}

void DistributionData::addPoint(double x, double y) {
    auto pair = make_pair(x, y);
    points.push_back(pair);
    pointsChanged = true;
}

DisVector DistributionData::getPoints() {
    return points;
}

vector<double> DistributionData::getDistributionParameters(int from, int to) {
    DisVector points = getStepRelativePoints();

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

    return vector<double>({mean, sigma, assymetry, excess});
}

DisVector DistributionData::getRelativePoints() {
    if (pointsChanged) {
        double summ = 0;
        for (pair<double, double> &p : points) {
            summ += p.second;
        }

        for (pair<double, double> &p : points) {
            relativePoints.push_back(make_pair(p.first, p.second/summ));
        }
        pointsChanged = false;
    }

    return relativePoints;
}

DisVector DistributionData::getStepRelativePoints() {
    DisVector stepRelativePoints;

    for (auto &p : getRelativePoints()) {
        stepRelativePoints.push_back(make_pair(p.first, p.second / 0.5));
    }

    return stepRelativePoints;
}
