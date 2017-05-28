#include "DataModel/distributiondata.h"
#include <math.h>

DistributionData::DistributionData()
{
    pointsChanged = true;  
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
    PointsVector points = getRelativePoints();
   // if (normalizedRelativePoints.size() > 0) {
   //     points = normalizedRelativePoints;
   // }
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
    double first = points.front().first;
    double last = points.back().first;

    double length = last - first;
    if (length <= 10) {
        if (pointsChanged) {
            relativePoints.clear();
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
    else {
        if (pointsChanged) {
            double min = points.front().first;
            double max = points.back().first;

            normalizedRelativePoints.clear();

            double summ = 0;
            for (pair<double, double> &p : points) {
                summ += p.second;
            }

            for (pair<double, double> &p : points) {
                double normalisedX = 10 * (p.first - min)/(max - min);
                normalizedRelativePoints.push_back(make_pair(normalisedX, p.second/summ));
            }
        }

        return normalizedRelativePoints;
    }

}

PointsVector DistributionData::getStepRelativePoints() {
    PointsVector stepRelativePoints;
    double step = getStep();
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
    double step = getStep();
    for (int i = from; i < to; i ++) {
        auto point = relPoints[i];
        stepRelativePoints.push_back(make_pair(point.first, point.second / step));
    }

    return stepRelativePoints;
}

double DistributionData::getDistributionSize() {
    double summ = 0;
    for (auto &p : points) {
        summ += p.second;
    }
    return summ;
}

double DistributionData::getStep()
{
    getRelativePoints();
    if (normalizedRelativePoints.size() > 0) {
        return normalizedRelativePoints[1].first - normalizedRelativePoints[0].first;
    }
    else {
        return relativePoints[1].first - relativePoints[0].first;
    }
    //return step;
}
