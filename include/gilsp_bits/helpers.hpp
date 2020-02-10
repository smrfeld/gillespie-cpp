//
//  helpers.hpp
//  gilsp
//
//  Created by Oliver Ernst on 2/10/20.
//

#include <math.h>       /* round, floor, ceil, trunc */

namespace gilsp {

const double AVOGADRO = 6.02214076 * pow(10,23);

int conc_in_um_to_no(double conc_in_um, double vol_in_l);

double no_to_conc_in_um(int no, double vol_in_l);

int binomial_coeff(int n, int k);
};
