#include "../include/gilsp_bits/rxn.hpp"

// Other headers
/*
 #include "../include/dblz_bits/ixn_param_traj.hpp"
 #include "../include/dblz_bits/ixn_param.hpp"
 #include "../include/dblz_bits/general.hpp"
 #include "../include/dblz_bits/diff_eq_rhs.hpp"
 #include "../include/dblz_bits/moment_diff.hpp"
 */

#include <iostream>
#include <fstream>

namespace gilsp {

// ***************
// MARK: - Constructor
// ***************

Rxn::Rxn(std::string name, double kr, std::vector<std::string> reactants, std::vector<std::string> products) {
    _name = name;
    _kr = kr;
    _r_list = reactants;
    _p_list = products;
    
    _update_species_and_multiplicities();
};
Rxn::Rxn(const Rxn& other) {
    _copy(other);
};
Rxn::Rxn(Rxn&& other) {
    _move(other);
};
Rxn& Rxn::operator=(const Rxn& other) {
    if (this != &other) {
        _clean_up();
        _copy(other);
    };
    return *this;
};
Rxn& Rxn::operator=(Rxn&& other) {
    if (this != &other) {
        _clean_up();
        _move(other);
    };
    return *this;
};
Rxn::~Rxn()
{
    _clean_up();
};
void Rxn::_clean_up() {
    // Nothing....
};
void Rxn::_copy(const Rxn& other) {
    _name = other._name;
    _kr = other._kr;
    _r_list = other._r_list;
    _p_list = other._p_list;
    _species = other._species;
    _species_r_multiplicity = other._species_r_multiplicity;
    _species_p_multiplicity = other._species_p_multiplicity;
};
void Rxn::_move(Rxn& other) {
    _name = other._name;
    _kr = other._kr;
    _r_list = other._r_list;
    _p_list = other._p_list;
    _species = other._species;
    _species_r_multiplicity = other._species_r_multiplicity;
    _species_p_multiplicity = other._species_p_multiplicity;
    
    // Reset the other
    other._name = "";
    other._kr = 0;
    other._r_list.clear();
    other._p_list.clear();
    other._species.clear();
    other._species_r_multiplicity.clear();
    other._species_p_multiplicity.clear();
};

// Internal update species and r/p multiplicity
void Rxn::_update_species_and_multiplicities() {
    // All species
    _species.clear();
    for (auto const &sp: _r_list) {
        if (!std::count(_species.begin(), _species.end(), sp)) {
            _species.push_back(sp);
        };
    };
    for (auto const &sp: _p_list) {
        if (!std::count(_species.begin(), _species.end(), sp)) {
            _species.push_back(sp);
        };
    };
    
    // Multiplicities
    _species_r_multiplicity.clear();
    _species_p_multiplicity.clear();
    for (auto const &sp: _species) {
        int r_count = std::count(_r_list.begin(), _r_list.end(), sp);
        if (r_count > 0) {
            _species_r_multiplicity.push_back(std::make_pair(sp,r_count));
        };
        
        int p_count = std::count(_p_list.begin(), _p_list.end(), sp);
        if (p_count > 0) {
            _species_p_multiplicity.push_back(std::make_pair(sp,p_count));
        };
    };
};

// ***************
// MARK: - Getters/setters
// ***************

std::string Rxn::get_name() const {
    return _name;
};
void Rxn::set_name(std::string name) {
    _name = name;
};

double Rxn::get_kr() const {
    return _kr;
};
void Rxn::set_kr(double kr) {
    _kr = kr;
};

const std::vector<std::string>& Rxn::get_reactants() const {
    return _r_list;
}
void Rxn::set_reactants(std::vector<std::string> reactants) {
    _r_list = reactants;
}
int Rxn::get_no_reactants() const {
    return _r_list.size();
}

const std::vector<std::string>& Rxn::get_products() const {
    return _p_list;
}
void Rxn::set_products(std::vector<std::string> products) {
    _p_list = products;
}
int Rxn::get_no_products() const {
    return _p_list.size();
}

const std::vector<std::pair<std::string,int>>& Rxn::get_reactant_multiplicity() const {
    return _species_r_multiplicity;
}
const std::vector<std::pair<std::string,int>>& Rxn::get_product_multiplicity() const {
    return _species_p_multiplicity;
}

const std::vector<std::string>& Rxn::get_species() const {
    return _species;
}

// ***************
// MARK: - Print
// ***************

void Rxn::print() const {
    std::cout << _name << std::endl;
    
    for (auto const &sp: _r_list) {
        std::cout << sp << " ";
    }
    std::cout << "-> ";
    for (auto const &sp: _p_list) {
        std::cout << sp << " ";
    }
    std::cout << "with rate: " << _kr << std::endl;
}

};
