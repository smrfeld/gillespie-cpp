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
    _propensities = other._propensities;
    _no_times_rxns_occur_in_tau = other._no_times_rxns_occur_in_tau;
    _propensities_sigma = other._propensities_sigma;
    
    _aux_mu = other._aux_mu;
    _aux_var = other._aux_var;
    _highest_order_event = other._highest_order_event;
    _epsilon_cao = other._epsilon_cao;
};
void TauLeaping::_move(TauLeaping& other) {
    _propensities = other._propensities;
    _no_times_rxns_occur_in_tau = other._no_times_rxns_occur_in_tau;
    _propensities_sigma = other._propensities_sigma;

    _aux_mu = other._aux_mu;
    _aux_var = other._aux_var;
    _highest_order_event = other._highest_order_event;
    _epsilon_cao = other._epsilon_cao;
    
    // Reset the other
    other._propensities.clear();
    other._no_times_rxns_occur_in_tau.clear();
    other._aux_var.clear();
    other._aux_mu.clear();
    other._highest_order_event.clear();
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
// MARK: - Tau step size selection
// ***************

// According to
//  Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151.

double TauLeaping::calculate_tau_step_size_cao(const std::vector<Rxn> &rxn_list, const Counts &counts) {
    
    // Calculate aux variables and find highest order event for every species
    
    // Reset them
    for (auto species: counts.get_species()) {
        _aux_mu[species] = 0.0;
        _aux_var[species] = 0.0;
        _highest_order_event[species] = 0.0;
    }
    
    // Go through rxns
    for (auto i=0; i<rxn_list.size(); i++) {
        for (auto const &pr: rxn_list.at(i).get_reactant_multiplicity()) {
            _aux_mu[pr.first] += -1 * pr.second * _propensities.at(i);
            _aux_var[pr.first] += pow(-1 * pr.second,2) * _propensities.at(i);
            
            // Highest order event
            if (_propensities.at(i) > _propensities.at(_highest_order_event.at(pr.first))) {
                _highest_order_event[pr.first] = i;
            }
        }
        for (auto const &pr: rxn_list.at(i).get_product_multiplicity()) {
            _aux_mu[pr.first] += pr.second * _propensities.at(i);
            _aux_var[pr.first] += pow(pr.second,2) * _propensities.at(i);
            
            // Highest order event
            if (_propensities.at(i) > _propensities.at(_highest_order_event.at(pr.first))) {
                _highest_order_event[pr.first] = i;
            }
        }
    }
    
    // Calculate time step tau
    double tau_min = 100000000000.0;
    for (auto species: counts.get_species()) {
        double count = counts.get_count(species);
        double highest_order_propensity = _propensities.at(_highest_order_event.at(species));
        double mu = _aux_mu.at(species);
        double var = _aux_var.at(species);
        
        double num = std::max(_epsilon_cao * count / highest_order_propensity, 1.0);
        double tau_1 = num / abs(mu);
        double tau_2 = pow(num, 2) / abs(var);
        
        tau_min = std::min(tau_min, std::min(tau_1, tau_2));
    }
    
    return tau_min;
}

// ***************
// MARK: - Calculate no times each reaction occurs in tau
// ***************

std::pair<bool,double> TauLeaping::calculate_no_times_reaction_occurs_in_tau(const std::vector<Rxn> &rxn_list, const Counts &counts) {

    if (rxn_list.size() != _propensities.size()) {
        _propensities = std::vector<double>(rxn_list.size(), 0);
    }
    for (auto i=0; i<_propensities.size(); i++) {
        if (rxn_list.at(i).get_no_reactants() != 0) {
            _propensities[i] = _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_reactant_multiplicity(), counts);
        } else {
            _propensities[i] = _get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_product_multiplicity(), counts);
        }
    }
    
    // Ensure there is any propensity
    double rates_total = std::accumulate(_propensities.begin(),
                                         _propensities.end(),
                                         decltype(_propensities)::value_type(0));
    if (rates_total < _propensities_sigma) {
        return std::make_pair(false,0.0);
    }
    
    // Compute tau
    double tau = calculate_tau_step_size_cao(rxn_list, counts);
    // std::cout << "Tau computed: " << tau << std::endl;
    
    // Number of times each event occurs in tau interval
    if (rxn_list.size() != _no_times_rxns_occur_in_tau.size()) {
        _no_times_rxns_occur_in_tau = std::vector<int>(rxn_list.size(), 0);
    }
    for (auto i=0; i<_no_times_rxns_occur_in_tau.size(); i++) {
        _no_times_rxns_occur_in_tau[i] = sample_poisson(_propensities.at(i) * tau);
    }
    
    return std::make_pair(true,tau);
}

// ***************
// MARK: - Do the reaction
// ***************

void TauLeaping::update_states_in_tau(const std::vector<Rxn> &rxn_list, Counts &counts) const {
    // Update the states
    for (auto i=0; i<rxn_list.size(); i++) {
        for (auto const &pr: rxn_list.at(i).get_reactant_multiplicity()) {
            counts.increment_count(pr.first, -1 * pr.second * _no_times_rxns_occur_in_tau.at(i), false);
        }
        for (auto const &pr: rxn_list.at(i).get_product_multiplicity()) {
            counts.increment_count(pr.first, pr.second * _no_times_rxns_occur_in_tau.at(i), false);
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
