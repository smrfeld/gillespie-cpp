#include <string>
#include <map>
#include <vector>

#include "counts.hpp"

namespace gilsp {

    // ***************
    // MARK: - CountsHist History
    // ***************
    
    class CountsHist {
        
    private:
                
        std::map<std::string, std::vector<int>> _counts_hist;
        std::vector<double> _t_hist;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const CountsHist& other);
        void _move(CountsHist &other);
        
    public:
                
        // ***************
        // MARK: - Constructor
        // ***************
        
        CountsHist();
        CountsHist(const CountsHist& other);
        CountsHist& operator=(const CountsHist& other);
        CountsHist(CountsHist&& other);
        CountsHist& operator=(CountsHist&& other);
        ~CountsHist();
        
        // ***************
        // MARK: - Get/set
        // ***************
        
        void store_counts(double time, const Counts &counts);
        
        std::pair<std::vector<double>, std::vector<int>> get_count_hist(std::string species) const;
        
        void write_count_hist(std::string species, std::string fname) const;
        
        void write_all_count_hist(std::string dir_name) const;
        
        void read_all_count_hists(std::string dir_name);
        
        std::map<std::string, double> get_average_counts(double time_start) const;
        
        double get_t_max() const;
        
        // ***************
        // MARK: - Print
        // ***************
        
        void print() const;
    };
};
