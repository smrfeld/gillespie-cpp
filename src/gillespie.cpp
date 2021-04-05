#include "../include/gilsp_bits/gillespie.hpp"
#include "../include/gilsp_bits/helpers.hpp"

#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <numeric>      // std::accumulate
#include "assert.h"

namespace gilsp {

// ***************
// MARK: - Constructor
// ***************

Gillespie::Gillespie() {
};
Gillespie::Gillespie(const Gillespie& other) {
    _copy(other);
};
Gillespie::Gillespie(Gillespie&& other) {
    _move(other);
};
Gillespie& Gillespie::operator=(const Gillespie& other) {
    if (this != &other) {
        _clean_up();
        _copy(other);
    };
    return *this;
};
Gillespie& Gillespie::operator=(Gillespie&& other) {
    if (this != &other) {
        _clean_up();
        _move(other);
    };
    return *this;
};
Gillespie::~Gillespie()
{
    _clean_up();
};
void Gillespie::_clean_up() {
    // Nothing....
};
void Gillespie::_copy(const Gillespie& other) {
    _props = other._props;
    _props_sigma = other._props_sigma;
};
void Gillespie::_move(Gillespie& other) {
    _props = other._props;
    _props_sigma = other._props_sigma;

    // Reset the other
    other._props.clear();
};

// ***************
// MARK: - Get propensity
// ***************

double Gillespie::_get_prop(double rate, const std::vector<std::pair<std::string,int>> &species_mult, const Counts &counts) {
    double prop = rate;
    
    for (auto const &pr: species_mult) {
        prop *= binomial_coeff_safe(counts.get_count(pr.first), pr.second);
    }
    
    return prop;
}

// ***************
// MARK: - Choose next reaction
// ***************

std::pair<Rxn const*, double> Gillespie::choose_next_rxn(const std::vector<Rxn> &rxn_list, const Counts &counts) {

    if (rxn_list.size() != _props.size()) {
        _props = std::vector<double>(rxn_list.size(), 0);
    }
    for (auto i=0; i<_props.size(); i++) {
        double prev = 0.0;
        if (i != 0) {
            prev = _props[i-1];
        }
        
        if (rxn_list.at(i).get_no_reactants() != 0) {
            _props[i] = prev + _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_reactant_multiplicity(), counts);
        } else {
            _props[i] = prev + _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_product_multiplicity(), counts);
        }
    }
    
    double props_cum = _props.back();
    if (props_cum < _props_sigma) {
        return std::make_pair(nullptr, 10000000.0);
    }
    
    // Choose the reaction
    double r = props_cum * (double) rand() / (RAND_MAX);
    
    // Find the lower bound
    auto low = std::lower_bound (_props.begin(), _props.end(), r);
    int idx = low - _props.begin();
    
    // Find time
    double t = (double) rand() / (RAND_MAX);
    double dt_next = log(1.0 / t) / props_cum;
    
    return std::make_pair(&rxn_list[idx], dt_next);
}

// ***************
// MARK: - Do the reaction
// ***************

void Gillespie::do_rxn(Rxn const* rxn, Counts &counts) {
    for (auto const &pr: rxn->get_reactant_multiplicity()) {
        counts.increment_count(pr.first, -1*pr.second, true);
    }
    for (auto const &pr: rxn->get_product_multiplicity()) {
        counts.increment_count(pr.first, pr.second, true);
    }
}

// ***************
// MARK: - Run and return the counts history
// ***************

CountsHist Gillespie::run(const std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, OptionsGillespie options) {
    
    // Check args
    if (options.write_as_we_go) {
        assert (options.write_dir != "");
    }
    
    // Setup writing as we go
    // ...
    
    // Store initial counts
    Counts counts_init(counts);
    
    // Setup count storage
    CountsHist counts_hist;
    counts_hist.store_counts(0, counts);
    double t_st_next = dt_st_every;
    
    // Run
    double t_curr = 0.0;
    while (t_curr <= t_max) {
    
        // Pick next rxn
        auto pr = choose_next_rxn(rxn_list, counts);
        
        // Check if exists
        if (pr.first == nullptr) {
            break;
        }
        
        // std::cout << "Next rxn time: " << pr.second << std::endl;
        
        // Store
        while ((t_st_next < t_curr + pr.second) && (t_st_next <= t_max)) {
            
            if (options.verbose) {
                std::cout << t_st_next << " / " << t_max << std::endl;
            }
            
            // Store
            counts_hist.store_counts(t_st_next, counts);
            
            // Write?
            //...
            
            // Advance
            t_st_next += dt_st_every;
        }
                
        // Advance time
        t_curr += pr.second;
        
        // Reached the end of time?
        if (t_curr > t_max) {
            break;
        }

        // Do the reaction
        do_rxn(pr.first, counts);
        
        // Conserve species
        for (auto const &sp: options.conserved_species) {
            counts.set_count(sp, counts_init.get_count(sp));
        }
    }
    
    // Fill up the rest of the counts
    while (t_st_next <= t_max) {
        
        if (options.verbose) {
            std::cout << t_st_next << " / " << t_max << std::endl;
        }
        
        // Store
        counts_hist.store_counts(t_st_next, counts);
                
        // Advance
        t_st_next += dt_st_every;
    }
    
    return counts_hist;
}

};
