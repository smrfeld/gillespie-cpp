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

struct OptionsTauLeaping {
    bool verbose=false;
    std::vector<std::string> conserved_species=std::vector<std::string>();
    bool write_as_we_go=false;
    std::string write_dir="";
    bool with_fixed_tau=false;
    double tau_fixed=0.1;
};

class TauLeaping {
    
private:
    
    // Rates
    std::vector<double> _propensities;
    
    // No times events occur in tau interval
    std::vector<int> _no_times_rxns_occur_in_tau;
        
    // Propensity sigma
    double _propensities_sigma = 1.0e-8;
    
    // Aux variables
    std::map<std::string, double> _aux_mu, _aux_var;
    
    // Highest order rate for every species
    std::map<std::string, int> _highest_order_event;
    
    // Epsilon
    double _epsilon_cao = 0.03;
    
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
    // MARK: - Tau step size selection
    // ***************

    // According to
    //  Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151.

    double calculate_tau_step_size_cao(const std::vector<Rxn> &rxn_list, const Counts &counts);

    // ***************
    // MARK: - Calculate no times each reaction occurs in tau
    // ***************
  
    std::pair<bool,double> calculate_no_times_reaction_occurs_in_tau(const std::vector<Rxn> &rxn_list, const Counts &counts, bool with_fixed_tau, double tau_fixed=0.1);

    // ***************
    // MARK: - Do the reaction
    // ***************
    
    void update_states_in_tau(const std::vector<Rxn> &rxn_list, Counts &counts) const;

    // ***************
    // MARK: - Run and return the counts history
    // ***************
    
    CountsHist run(const std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, OptionsTauLeaping options=OptionsTauLeaping());
    };
};

#endif
