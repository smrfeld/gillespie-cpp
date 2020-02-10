#include <string>
#include <vector>

#include "rxn.hpp"
#include "counts_hist.hpp"

namespace gilsp {

std::pair<Rxn*, double> choose_next_rxn(std::vector<Rxn> &rxn_list, Counts &counts);

void do_rxn(Rxn* rxn, Counts &counts);

CountsHist run_gillespie(std::vector<Rxn> &rxn_list, Counts &counts, double dt_st_every, double t_max, bool verbose=false, std::vector<std::string> conserved_species=std::vector<std::string>(), bool write_as_we_go=false, std::string write_dir="");
};
