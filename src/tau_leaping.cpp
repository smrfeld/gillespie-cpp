#include "../include/gilsp_bits/tau_leaping.hpp"
#include "../include/gilsp_bits/helpers.hpp"

#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <numeric>      // std::accumulate
#include "assert.h"

namespace gilsp {

// ***************
// MARK: - Constructor
// ***************

TauLeaping::TauLeaping() {
};
TauLeaping::TauLeaping(const TauLeaping& other) {
    _copy(other);
};
TauLeaping::TauLeaping(TauLeaping&& other) {
    _move(other);
};
TauLeaping& TauLeaping::operator=(const TauLeaping& other) {
    if (this != &other) {
        _clean_up();
        _copy(other);
    };
    return *this;
};
TauLeaping& TauLeaping::operator=(TauLeaping&& other) {
    if (this != &other) {
        _clean_up();
        _move(other);
    };
    return *this;
};
TauLeaping::~TauLeaping()
{
    _clean_up();
};
void TauLeaping::_clean_up() {
    // Nothing....
};
void TauLeaping::_copy(const TauLeaping& other) {
    _rates = other._rates;
    _no_times_rxns_occur_in_tau = other._no_times_rxns_occur_in_tau;
    _rates_sigma = other._rates_sigma;
};
void TauLeaping::_move(TauLeaping& other) {
    _rates = other._rates;
    _no_times_rxns_occur_in_tau = other._no_times_rxns_occur_in_tau;
    _rates_sigma = other._rates_sigma;

    // Reset the other
    other._rates.clear();
    other._no_times_rxns_occur_in_tau.clear();
};

// ***************
// MARK: - Get propensity
// ***************

double TauLeaping::_get_prop(double rate, const std::vector<std::pair<std::string,int>> &species_mult, const Counts &counts) {
    double prop = rate;
    
    for (auto const &pr: species_mult) {
        prop *= binomial_coeff_safe(counts.get_count(pr.first), pr.second);
    }
    
    return prop;
}

// ***************
// MARK: - Choose next reaction
// ***************

std::pair<bool,double> TauLeaping::calculate_no_times_reaction_occurs_in_tau(const std::vector<Rxn> &rxn_list, const Counts &counts) {

    if (rxn_list.size() != _rates.size()) {
        _rates = std::vector<double>(rxn_list.size(), 0);
    }
    for (auto i=0; i<_rates.size(); i++) {
        if (rxn_list.at(i).get_no_reactants() != 0) {
            _rates[i] = _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_reactant_multiplicity(), counts);
        } else {
            _rates[i] = _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_product_multiplicity(), counts);
        }
    }
    
    // Normalize
    double rates_total = std::accumulate(_rates.begin(), _rates.end(),
                                   decltype(_rates)::value_type(0));
    
    if (rates_total < _rates_sigma) {
        return std::make_pair(false,0.0);
    }
    
    for (auto i=0; i<_rates.size(); i++) {
        _rates[i] /= rates_total;
    }

    // Number of times each event occurs in tau interval
    double tau = 1.0;
    if (rxn_list.size() != _no_times_rxns_occur_in_tau.size()) {
        _no_times_rxns_occur_in_tau = std::vector<int>(rxn_list.size(), 0);
    }
    for (auto i=0; i<_no_times_rxns_occur_in_tau.size(); i++) {
        _no_times_rxns_occur_in_tau[i] = sample_poisson(_rates.at(i) * tau);
    }
    
    return std::make_pair(true,tau);
}

// ***************
// MARK: - Do the reaction
// ***************

void TauLeaping::update_states_in_tau(const std::vector<Rxn> &rxn_list, Counts &counts) {
    // Update the states
    for (auto i=0; i<rxn_list.size(); i++) {
        for (auto const &pr: rxn_list.at(i).get_reactant_multiplicity()) {
            counts.increment_count(pr.first, -1 * pr.second * _no_times_rxns_occur_in_tau.at(i));
        }
        for (auto const &pr: rxn_list.at(i).get_product_multiplicity()) {
            counts.increment_count(pr.first, pr.second * _no_times_rxns_occur_in_tau.at(i));
        }
    }
}

// ***************
// MARK: - Run and return the counts history
// ***************

CountsHist TauLeaping::run(const std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, bool verbose, std::vector<std::string> conserved_species, bool write_as_we_go, std::string write_dir) {
    
    // Check args
    if (write_as_we_go) {
        assert (write_dir != "");
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
    
        auto pr0 = calculate_no_times_reaction_occurs_in_tau(rxn_list, counts);
        if (!pr0.first) {
            // Stop, no rxns left
            break;
        }
        double tau = pr0.second;
                
        // Store current counts
        while ((t_st_next < t_curr + tau) && (t_st_next <= t_max)) {
            
            if (verbose) {
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
        t_curr += tau;
        
        // Reached the end of time?
        if (t_curr > t_max) {
            break;
        }
        
        // Update states = do rxns
        update_states_in_tau(rxn_list, counts);
                
        // Conserve species
        for (auto const &sp: conserved_species) {
            counts.set_count(sp, counts_init.get_count(sp));
        }
    }
    
    // Fill up the rest of the counts
    while (t_st_next <= t_max) {
        
        if (verbose) {
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
