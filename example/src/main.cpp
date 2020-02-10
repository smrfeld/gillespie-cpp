#include <iostream>
#include <vector>

#include <gilsp>

using namespace gilsp;
using namespace std;

int main() {
    
    /* initialize random seed: */
    srand (time(NULL));

    vector<string> reactants({"A"});
    vector<string> products;
    Rxn rxn = Rxn("rxn", 3.0, reactants, products);
    
    Counts counts = Counts();
    counts.set_count("A", 100);
 
    vector<Rxn> rxn_list({rxn});
    double dt_st_every = 0.1;
    double t_max = 100.0;
    auto counts_hist = run_gillespie(rxn_list, counts, dt_st_every, t_max);
    
    counts_hist.print();
    
    return 0;
}
