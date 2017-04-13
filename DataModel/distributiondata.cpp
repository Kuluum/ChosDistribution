#include "DataModel/distributiondata.h"
#include <math.h>

DistributionData::DistributionData()
{
    pointsChanged = true;
    step = 1;
}

void DistributionData::setStep(double step)
{
    this->step = step;
}

void DistributionData::setPoints(PointsVector points) {
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

PointsVector DistributionData::getPoints() {
    return points;
}

vector<double> DistributionData::getDistributionParameters(int from, int to) {
    PointsVector points = getStepRelativePoints();

    double fSumm = 0.0;
    double xfSumm = 0.0;
    int m = points.size();
    if (to > m) {
        to = m;
    }
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

PointsVector DistributionData::getRelativePoints() {
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

PointsVector DistributionData::getStepRelativePoints() {
    PointsVector stepRelativePoints;

    for (auto &p : getRelativePoints()) {
        stepRelativePoints.push_back(make_pair(p.first, p.second / step));
    }

    return stepRelativePoints;
}

PointsVector DistributionData::getStepRelativePoints(int from, int to) {
    PointsVector stepRelativePoints;

    if (from < 0)
    {
        from = 0;
    }

    if (to < 0)
    {
        to = 0;
    }

    if (to > points.size())
    {
        to = points.size();
    }

    auto relPoints = getRelativePoints();

    for (int i = from; i < to; i ++) {
        auto point = relPoints[i];
        stepRelativePoints.push_back(make_pair(point.first, point.second / step));
    }

    return stepRelativePoints;
}

