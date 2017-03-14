#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <functional>
#include <vector>
#include <utility>


#include <chrono>
#include <iostream>
struct profiler
{
    std::string name;
    std::chrono::high_resolution_clock::time_point p;
    profiler(std::string const &n) :
        name(n), p(std::chrono::high_resolution_clock::now()) { }
    ~profiler()
    {
        using dura = std::chrono::duration<double>;
        auto d = std::chrono::high_resolution_clock::now() - p;
        std::cout << name << ": "
            << std::chrono::duration_cast<dura>(d).count()
            << std::endl;
    }
};

#define PROFILE_BLOCK(pbn) profiler _pfinstance(pbn)


using namespace std;

typedef vector<vector<double>> matrix;

class Algorithms
{
public:
    Algorithms();
    static double RSS(vector<pair<double, double>> &dataVector, const function<double (double)> func);
//    static matrix multMatrix(const matrix &m1, const matrix &m2);
    static vector<double> hessXgrad(const matrix &hess, const vector<double> &grad);
    static matrix mesaInvertMatrix(const matrix &inp);
};

#endif // ALGORITHMS_H
