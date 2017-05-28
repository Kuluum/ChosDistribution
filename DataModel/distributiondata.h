#ifndef DISTRIBUTIONDATA_H
#define DISTRIBUTIONDATA_H

#include <vector>
#include <utility>

using namespace std;

typedef pair<double, double> Point;
typedef vector<Point> PointsVector;

class DistributionData
{
public:
    DistributionData();
    void addPoint(pair<double, double> point);
    void addPoint(double x, double y);
    void setPoints(PointsVector points);

    PointsVector getPoints();
    PointsVector getRelativePoints();
    PointsVector getStepRelativePoints();
    PointsVector getStepRelativePoints(int from, int to);
    vector<double> getDistributionParameters(int, int);
    double getStep();
    double getDistributionSize();

private:
    PointsVector points;
    PointsVector relativePoints;
    PointsVector normalizedRelativePoints;
    bool pointsChanged;
};

#endif // DISTRIBUTIONDATA_H
