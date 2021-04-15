#include <iostream>
#include <vector>
#include <filesystem>
#include <random>

#include <gilsp>

using namespace gilsp;
using namespace std;

int main() {

    // Initial count
    int no_init_prey = 241;
    int no_init_hunter = 602;
        
    // Reactions
    vector<Rxn> rxn_list;
    rxn_list.push_back(Rxn("1", 1.0, {"P"}, {"P","P"}));
    rxn_list.push_back(Rxn("2", 0.2, {"H"}, {}));
    rxn_list.push_back(Rxn("3", 0.0016, {"H", "P"}, {"H","H"}));

    // Run
    std::string dir_name = "../data_bimol_gillespie/";
    std::filesystem::create_directories(dir_name);
    Gillespie g;
    for (auto seed=0; seed<10; seed++) {
        
        // Random seed
        srand(seed);
        
        // Reset initial counts
        Counts counts;
        counts.set_count("H",no_init_hunter);
        counts.set_count("P",no_init_prey);
        
        double dt_st_every = 0.25;
        double t_max = 30.0;
        CountsHist counts_hist = g.run(rxn_list, counts, dt_st_every, t_max);

        // Write
        std::string fname = dir_name + pad_str(seed, 3) + ".txt";
        counts_hist.write_count_hist_all_species_single_file(fname);
        
    };
    
    return 0;
}
