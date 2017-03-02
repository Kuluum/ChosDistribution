#include "Algorithm/algorithms.h"

Algorithms::Algorithms()
{

}

double Algorithms::RSS(vector<pair<double, double>> dataVector, const function<double (double)> func)
{
    double rss = 0;
    for (auto p : dataVector) {
        double diff = p.second - func(p.first);
        rss += diff * diff;
    }
    return rss;
}

matrix Algorithms::multMatrix(matrix &m1, matrix &m2)
{
    size_t m1_xSize = m1.size();
    size_t m1_ySize = m1[0].size();

    size_t m2_xSize = m2.size();
    size_t m2_ySize = m2[0].size();

    matrix result = matrix(m1_xSize, vector<double>(m2_ySize, 0.0));

    for (int i = 0; i < m1_xSize; i++)
    {
       for (int j = 0; j < m2_ySize; j++)
       {
          result[i][j] = 0;
          for (int k = 0; k < m1_ySize; k++)
          {
             result[i][j] += m1[i][k] * m2[k][j];
          }
       }
    }

    return result;
}

