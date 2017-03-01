#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <functional>
#include <vector>
#include <utility>

class Algorithms
{
public:
    Algorithms();
    static double RSS(std::vector<std::pair<double, double>> dataVector, const std::function<double (double)> func);
    //static constexpr std::vector<std::vector> linedHess(std::vector grad);
};

#endif // ALGORITHMS_H
