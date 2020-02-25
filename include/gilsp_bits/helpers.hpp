#include <math.h>       /* round, floor, ceil, trunc */
#include <boost/math/special_functions/binomial.hpp>

namespace gilsp {

const double AVOGADRO = 6.02214076 * pow(10,23);

int conc_in_um_to_no(double conc_in_um, double vol_in_l);
double conc_in_um_to_no_as_double(double conc_in_um, double vol_in_l);

double no_to_conc_in_um(int no, double vol_in_l);

double binomial_coeff_safe(int n, int k);
};
