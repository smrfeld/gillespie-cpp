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

};
