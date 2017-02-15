#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <functional>
#include <QVector>
#include <QPair>

class Algorithms
{
public:
    Algorithms();
    static double RSS(QVector<QPair<double, double>> dataVector, const std::function<double (double)> func);
};

#endif // ALGORITHMS_H
