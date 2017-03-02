#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <functional>
#include <vector>
#include <utility>

using namespace std;

typedef vector<vector<double>> matrix;

class Algorithms
{
public:
    Algorithms();
    static double RSS(vector<pair<double, double>> dataVector, const function<double (double)> func);
    static matrix multMatrix(matrix &m1, matrix &m2);
    //static constexpr std::vector<std::vector> linedHess(std::vector grad);
};

#endif // ALGORITHMS_H
