#include <string>
#include <map>
#include <vector>

#ifndef COUNTS_H
#define COUNTS_H

namespace gilsp {

    // ***************
    // MARK: - Counts
    // ***************
    
    class Counts {
        
    private:
                
        std::map<std::string, int> _counts;
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const Counts& other);
        void _move(Counts &other);
        
    public:
                
        // ***************
        // MARK: - Constructor
        // ***************
        
        Counts();
        Counts(const Counts& other);
        Counts& operator=(const Counts& other);
        Counts(Counts&& other);
        Counts& operator=(Counts&& other);
        ~Counts();
        
        // ***************
        // MARK: - Get/set
        // ***************
        
        void set_count(std::string species, int count);
        
        void increment_count(std::string species, int increment, bool require_non_negative);
        
        void increment_count(const std::map<std::string,int> &increments, bool require_non_negative);
        
        int get_count(std::string species) const;
        
        std::vector<std::string> get_species() const;
        
        // ***************
        // MARK: - Print
        // ***************
        
        void print() const;
    };
};

#endif
