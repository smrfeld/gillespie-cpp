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

double conc_in_um_to_no_as_double(double conc_in_um, double vol_in_l)
{
    return conc_in_um * pow(10, -6) * vol_in_l * AVOGADRO;
};

double no_to_conc_in_um(int no, double vol_in_l) {
    return no / (pow(10,-6) * vol_in_l * AVOGADRO);
};

double binomial_coeff_safe(int n, int k) {
    if (k == 1) {
        return n;
    } else if (k <= n) {
        return boost::math::binomial_coefficient<double>(n,k);
    } else {
        return 0;
    }
};


std::string pad_str(int i, int n_zeros) {
    std::stringstream fname;
    fname << std::setfill('0') << std::setw(n_zeros) << i;
    return fname.str();
};

std::string prec_str(double x, int prec) {
    std::stringstream fname;
    fname << std::fixed << std::setprecision(prec) << x;
    return fname.str();
};

std::string prec_p_str(double x, int prec) {
    auto str = prec_str(x,prec);
    for (int i = 0; i < str.length(); ++i) {
        if (str[i] == "."[0]) {
            str[i] = "p"[0];
        }
    }
    
    return str;
};

int sample_poisson(double lambda) {
    double L = exp(-lambda);
    int k = 0;
    double p = 1;
    do {
        k++;
        
        // Random in [0,1]
        double u = (double) rand() / (RAND_MAX);

        p *= u;

    } while (p > L);
    
    return k-1;
}

};
