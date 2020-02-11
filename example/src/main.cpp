#include <iostream>
#include <vector>

#include <gilsp>

using namespace gilsp;
using namespace std;

int main() {
    
    // Random seed
    srand (time(NULL));

    // Make the reaction A->0 with rate 3
    vector<string> reactants({"A"});
    vector<string> products;
    Rxn rxn = Rxn("rxn", 3.0, reactants, products);
    
    // Initial counts of particles
    Counts counts = Counts();
    counts.set_count("A", 100);
 
    // Run
    vector<Rxn> rxn_list({rxn});
    double dt_st_every = 0.1;
    double t_max = 100.0;
    CountsHist counts_hist = run_gillespie(rxn_list, counts, dt_st_every, t_max);
    
    // Print it
    counts_hist.print();
    
    // Write it
    counts_hist.write_all_count_hist("../data");
    
    return 0;
}
