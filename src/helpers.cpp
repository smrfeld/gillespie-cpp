//
//  helpers.cpp
//  gilsp
//
//  Created by Oliver Ernst on 2/10/20.
//

#include "helpers.hpp"

namespace gilsp {

int conc_in_um_to_no(double conc_in_um, double vol_in_l)
{
    return round(conc_in_um * pow(10, -6) * vol_in_l * AVOGADRO);
};

double no_to_conc_in_um(int no, double vol_in_l) {
    return no / (pow(10,-6) * vol_in_l * AVOGADRO);
};
  
// Returns value of Binomial Coefficient C(n, k)
int binomial_coeff(int n, int k) {
    int res = 1;
  
    if (k > n) {
        return 0;
    }
    
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k ) {
        k = n - k;
    }
  
    // Calculate value of
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
  
    return res;
}


};
