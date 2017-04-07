#ifndef DISTRIBUTIONDATA_H
#define DISTRIBUTIONDATA_H

#include <vector>
#include <utility>

using namespace std;

typedef vector<pair<double, double>> DisVector;

class DistributionData
{
public:
    DistributionData();
    void addPoint(pair<double, double> point);
    void addPoint(double x, double y);
    void setPoints(DisVector points);

    DisVector getPoints();
    DisVector getRelativePoints();
    DisVector getStepRelativePoints();
    vector<double> getDistributionParameters(int, int);

private:
    DisVector points;
    DisVector relativePoints;
    double step;
    bool pointsChanged;
};

#endif // DISTRIBUTIONDATA_H
