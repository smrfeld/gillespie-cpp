#include <string>
#include <vector>

#include "rxn.hpp"
#include "counts_hist.hpp"

namespace gilsp {

// ***************
// MARK: - Gillespie
// ***************

class Gillespie {
    
private:
            
    // Propensity
    std::vector<double> _props;
    
    // Propensity sigma
    double _props_sigma = 1.0e-8;
    
    // Get propensity
    double _get_prop(double rate, const std::map<std::string,int> &species_mult, const Counts &counts);
    
    // Internal copy func/clean up
    void _clean_up();
    void _copy(const Gillespie& other);
    void _move(Gillespie &other);
    
public:
            
    // ***************
    // MARK: - Constructor
    // ***************
    
    Gillespie();
    Gillespie(const Gillespie& other);
    Gillespie& operator=(const Gillespie& other);
    Gillespie(Gillespie&& other);
    Gillespie& operator=(Gillespie&& other);
    ~Gillespie();
    
    // ***************
    // MARK: - Choose next reaction
    // ***************
  
    std::pair<Rxn const*, double> choose_next_rxn(const std::vector<Rxn> &rxn_list, const Counts &counts);

    // ***************
    // MARK: - Do the reaction
    // ***************
    
    void do_rxn(Rxn const* rxn, Counts &counts);

    // ***************
    // MARK: - Run and return the counts history
    // ***************
    
    CountsHist run(const std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, bool verbose=false, std::vector<std::string> conserved_species=std::vector<std::string>(), bool write_as_we_go=false, std::string write_dir="");
    };
};
