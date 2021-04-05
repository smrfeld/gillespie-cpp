#include <iostream>
#include <vector>

#include <gilsp>

using namespace gilsp;
using namespace std;

int main() {
    
    // Random seed
    srand (time(NULL));

    // Make the reaction A->0 with rate 3
    Rxn rxn = Rxn("rxn", 3.0, {"A"}, {});
    
    // Initial counts of particles
    Counts counts = Counts();
    counts.set_count("A", 100);
 
    // Run
    TauLeaping g;
    vector<Rxn> rxn_list({rxn});
    double dt_st_every = 0.1;
    double t_max = 100.0;
    CountsHist counts_hist = g.run(rxn_list, counts, dt_st_every, t_max);

    // Print it
    counts_hist.print();
    
    // Write it
    counts_hist.write_count_hist_all_species("../data_tau_leaping");
    
    return 0;
}
