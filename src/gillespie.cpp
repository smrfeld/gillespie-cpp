#include "../include/gilsp_bits/gillespie.hpp"
#include "../include/gilsp_bits/helpers.hpp"

#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <numeric>      // std::accumulate

namespace gilsp {

double get_prop(double rate, std::map<std::string,int> species_mult, Counts &counts) {
    double prop = rate;
    
    for (auto const &pr: species_mult) {
        int count = counts.get_count(pr.first);
        prop *= binomial_coeff(count, pr.second);
    }
    
    return prop;
}

std::pair<Rxn*, double> choose_next_rxn(std::vector<Rxn> &rxn_list, Counts &counts) {

    std::vector<double> props(rxn_list.size(), 0);
    for (auto i=0; i<rxn_list.size(); i++) {
        double prev = 0.0;
        if (i != 0) {
            prev = props[i-1];
        }
        
        if (rxn_list.at(i).get_reactants().size() != 0) {
            props[i] = prev + get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_reactant_multiplicity(), counts);
        } else {
            props[i] = prev + get_prop(rxn_list.at(i).get_kr(), rxn_list.at(i).get_product_multiplicity(), counts);
        }
    }
    
    double props_cum = props.back();
    if (props_cum < 1.0e-8) {
        return std::make_pair(nullptr, 10000000.0);
    }
    
    // Choose the reaction
    double r = props_cum * (double) rand() / (RAND_MAX);
    
    // Find the lower bound
    auto low = std::lower_bound (props.begin(), props.end(), r);
    int idx = low - props.begin();
    
    // Find time
    double t = (double) rand() / (RAND_MAX);
    double dt_next = log(1.0 / t) / props_cum;
    
    return std::make_pair(&rxn_list[idx], dt_next);
}

void do_rxn(Rxn* rxn, Counts &counts) {
    for (auto const &pr: rxn->get_reactant_multiplicity()) {
        counts.increment_count(pr.first, -1*pr.second);
    }
    for (auto const &pr: rxn->get_product_multiplicity()) {
        counts.increment_count(pr.first, pr.second);
    }
}

CountsHist run_gillespie(std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, bool verbose, std::vector<std::string> conserved_species, bool write_as_we_go, std::string write_dir) {
    
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
    
        // Pick next rxn
        auto pr = choose_next_rxn(rxn_list, counts);
        
        // Check if exists
        if (pr.first == nullptr) {
            break;
        }
        
        // Store
        while (t_st_next < t_curr + pr.second) {
            
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
        t_curr += pr.second;
        
        // Do the reaction
        do_rxn(pr.first, counts);
        
        // Conserve species
        for (auto const &sp: conserved_species) {
            counts.set_count(sp, counts_init.get_count(sp));
        }
    }
    
    return counts_hist;
}

};
