#include <string>
#include <vector>

#include "rxn.hpp"
#include "counts_hist.hpp"

#ifndef TAU_LEAPING_H
#define TAU_LEAPING_H

namespace gilsp {

// ***************
// MARK: - TauLeaping
// ***************

class TauLeaping {
    
private:
            
    // Rates
    std::vector<double> _rates;
    
    // No times events occur in tau interval
    std::vector<int> _no_times_rxns_occur_in_tau;
    
    // Propensity sigma
    double _rates_sigma = 1.0e-8;
    
    // Get propensity
    double _get_prop(double rate, const std::vector<std::pair<std::string,int>> &species_mult, const Counts &counts);
        
    // Internal copy func/clean up
    void _clean_up();
    void _copy(const TauLeaping& other);
    void _move(TauLeaping &other);
    
public:
            
    // ***************
    // MARK: - Constructor
    // ***************
    
    TauLeaping();
    TauLeaping(const TauLeaping& other);
    TauLeaping& operator=(const TauLeaping& other);
    TauLeaping(TauLeaping&& other);
    TauLeaping& operator=(TauLeaping&& other);
    ~TauLeaping();
    
    // ***************
    // MARK: - Choose next reaction
    // ***************
  
    std::pair<bool,double> calculate_no_times_reaction_occurs_in_tau(const std::vector<Rxn> &rxn_list, const Counts &counts);

    // ***************
    // MARK: - Do the reaction
    // ***************
    
    void update_states_in_tau(const std::vector<Rxn> &rxn_list, Counts &counts);

    // ***************
    // MARK: - Run and return the counts history
    // ***************
    
    CountsHist run(const std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, bool verbose=false, std::vector<std::string> conserved_species=std::vector<std::string>(), bool write_as_we_go=false, std::string write_dir="");
    };
};

#endif
