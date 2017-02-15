#ifndef DISTRIBUTIONDATA_H
#define DISTRIBUTIONDATA_H

#include <QVector>
#include <QPair>

typedef QVector<QPair<double, double>> DisVector;
class DistributionData
{
public:
    DistributionData();
    void addPoint(QPair<double, double> point);
    void addPoint(double x, double y);
    void setPoints(QVector<QPair<double, double>> points);

    DisVector getPoints();
    DisVector getRelativePoints();
    DisVector getStepRelativePoints();
    QVector<double> getDistributionParameters(int, int);

private:
    DisVector points;
    DisVector relativePoints;
    bool pointsChanged;
};

#endif // DISTRIBUTIONDATA_H
