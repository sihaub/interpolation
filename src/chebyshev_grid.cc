#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{
    StandardGrid::StandardGrid(size_t p)
    {
        _p = p;

        for (size_t j = 0; j <= p; j++){
            _tj.push_back(cos(j * M_PI / static_cast<double>(p)));
        } 

        for (size_t j = 0; j <= p; j++) {
            double sign = j == 0 % 2 ? +1 : -1;
            double scaling;
            if (j==0 || j == p) scaling = 0.5;
            _betaj.push_back(sign * scaling);
        }

        _Dij.resize(_p + 1, vector_d(_p+1, 0.));
        _Dij[0][0] = (2. * _p * _p + 1) / 6.;
        _Dij[_p][_p] = -_Dij[0][0];

        for(size_t j=1; j<p; j++) {
            _Dij[j][j] = -0.5 * _tj[j] / (1. - pow(_tj[j], 2));
        }

        for (size_t i = 0; i <= p; i++) {
            for (size_t j = 0; j <= p; j++) {
                if (j == i) continue;
                _Dij[i][j] = -(_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
            }
        }
    }


double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
    if(t<-1 || t>1){
        throw std::domain_error("StandardGrid::interpolate t must be in [-1,1]");
    }
    if (end - start != _p) {
        throw std::domain_error("StandardGrid::interpolate end-start should be equal to p");
    }

    double den =0.;
    for (size_t j=0; j <= _p; j++) {
        if (t == _tj[j]) return fj[j + start] 
        den += _betaj[j] / (t - _tj[j]);
    }
 
}    double res = 0.;
    for (size_t i=0; i <= _p; i++)
{
    
}
}

} // namespace Chebyshev
} // namespace Interpolation