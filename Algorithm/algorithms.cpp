#include "Algorithm/algorithms.h"

Algorithms::Algorithms()
{

}

double Algorithms::RSS(std::vector<std::pair<double, double>> dataVector, const std::function<double (double)> func)
{
    double rss = 0;
    for (auto p : dataVector) {
        double diff = p.second - func(p.first);
        rss += diff * diff;
    }
    return rss;
}

//std::vector<std::vector> Algorithms::linedHess(std::vector grad)
//{
//    size_t size = grad.size();
//    std::vector<std::vector> hess(size, std::vector(size));

//       for(int i = 0; i < size; ++i)
//       {
//           for(int j = 0; j < size; ++j)
//           {
//               hess[i][j] ;
//           }
//       }
//}
