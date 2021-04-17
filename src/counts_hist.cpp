#include "../include/gilsp_bits/counts_hist.hpp"
#include "../include/gilsp_bits/helpers.hpp"

#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>
#include <sstream>
#include "assert.h"
#include <filesystem>

namespace gilsp {

// ***************
// MARK: - Constructor
// ***************

CountsHist::CountsHist() {
};
CountsHist::CountsHist(const CountsHist& other) {
    _copy(other);
};
CountsHist::CountsHist(CountsHist&& other) {
    _move(other);
};
CountsHist& CountsHist::operator=(const CountsHist& other) {
    if (this != &other) {
        _clean_up();
        _copy(other);
    };
    return *this;
};
CountsHist& CountsHist::operator=(CountsHist&& other) {
    if (this != &other) {
        _clean_up();
        _move(other);
    };
    return *this;
};
CountsHist::~CountsHist()
{
    _clean_up();
};
void CountsHist::_clean_up() {
    // Nothing....
};
void CountsHist::_copy(const CountsHist& other) {
    _counts_hist = other._counts_hist;
    _t_hist = other._t_hist;
};
void CountsHist::_move(CountsHist& other) {
    _counts_hist = other._counts_hist;
    _t_hist = other._t_hist;

    // Reset the other
    other._counts_hist.clear();
    other._t_hist.clear();
};

// ***************
// MARK: - Get/set
// ***************

void CountsHist::store_counts(double time, const Counts &counts) {
    _t_hist.push_back(time);
    
    for (auto const &species: counts.get_species()) {
        auto count = counts.get_count(species);
        
        auto it = _counts_hist.find(species);
        if (it == _counts_hist.end()) {
            _counts_hist[species] = std::vector<int>();
        }
        
        _counts_hist[species].push_back(count);
        
        assert (_counts_hist[species].size() == _t_hist.size());
    }
}

std::pair<std::vector<double>, std::vector<int>> CountsHist::get_count_hist(std::string species) const {
    return std::make_pair(_t_hist, _counts_hist.at(species));
}

void CountsHist::write_count_hist_single_species(std::string species, std::string fname) const {
    std::ofstream f;
    f.open(fname);
    
    assert (f.is_open());

    // Go through time
    for (auto t=0; t<_t_hist.size(); t++) {
        f << _t_hist.at(t) << " " << _counts_hist.at(species).at(t) << "\n";
    }
    f.close();
}

void CountsHist::write_count_hist_all_species_seperate_files_by_species(std::string dir_name) const {
    for (auto const &pr: _counts_hist) {
        if (std::string(1,dir_name.back()) != "/") {
            write_count_hist_single_species(pr.first, dir_name + "/" + pr.first + ".txt");
        } else {
            write_count_hist_single_species(pr.first, dir_name + pr.first + ".txt");
        }
    }
}

void CountsHist::write_count_hist_all_species_single_file(std::string fname) const {
    std::ofstream f;
    f.open(fname);
    assert (f.is_open());
    
    // Write species labels - map is sorted by default
    f << "t";
    for (auto const &pr: _counts_hist) {
        f << " " << pr.first;
    }
    f << "\n";
    
    // Write values in time
    for (auto t=0; t<_t_hist.size(); t++) {
        f << _t_hist.at(t);
        for (auto const &pr: _counts_hist) {
            f << " " << pr.second.at(t);
        }
        f << "\n";
    }
    
    f.close();
}

void CountsHist::write_count_hist_all_species_seperate_files_by_timepoint(std::string dir_name) const {
    std::ofstream f;
    std::string dir_name_mod;
    if (std::string(1,dir_name.back()) != "/") {
        dir_name_mod = dir_name + "/";
    } else {
        dir_name_mod = dir_name;
    }
    
    for (auto t=0; t<_t_hist.size(); t++) {
        f.open(dir_name_mod + pad_str(t,5) + ".txt");
        
        assert (f.is_open());

        for (auto const &pr: _counts_hist) {
            f << pr.first << " " << pr.second.at(t) << "\n";
        }
        
        f.close();
    }
}


void CountsHist::read_count_hist_all_species(std::string dir_name) {
    std::ifstream f;
    
    std::string time="";
    std::string count="";
    std::istringstream iss;
    std::string line;
    std::string species;
    
    // Clear
    _t_hist.clear();
    _counts_hist.clear();

    bool read_time = true;
    for (const auto & fname : std::filesystem::directory_iterator(dir_name)) {

        // Check name
        if (fname.path().extension() != ".txt") {
            continue;
        }
        
        f.open(fname.path());
        
        species = fname.path().stem();
        _counts_hist[species] = std::vector<int>();
        
        assert (f.is_open());

        while (getline(f,line)) {
            
            if (line == "") { continue; };
            iss = std::istringstream(line);
            iss >> time >> count;
            
            // Append
            if (read_time) {
                _t_hist.push_back(atof(time.c_str()));
                // std::cout << "_t_hist: " << _t_hist.back() << std::endl;
            }
            _counts_hist[species].push_back(atoi(count.c_str()));
            // std::cout << "counts: " << species << " " << _counts_hist[species].back() << std::endl;

            // Reset
            time=""; count="";
        }
        
        f.close();
        
        // Only read time once
        read_time = false;
    }
    
    // Check lengths
    for (auto const &pr: _counts_hist) {
        assert (_t_hist.size() == pr.second.size());
    }
}

std::map<std::string, double> CountsHist::get_average_counts(double time_start) const {
    auto low = std::lower_bound (_t_hist.begin(), _t_hist.end(), time_start);
    int idx = low - _t_hist.begin();
    
    std::map<std::string, double> counts_average;
    for (auto const &pr: _counts_hist) {
        // Average
        double val = 0.0;
        for (auto i=idx; i<pr.second.size(); i++) {
            val += pr.second[i];
        }
        val /= pr.second.size() - idx;
        
        // Store
        counts_average[pr.first] = val;
    }
    
    return counts_average;
}

double CountsHist::get_t_max() const {
    return _t_hist.back();
}

// ***************
// MARK: - Print
// ***************

void CountsHist::print() const {
    for (auto const &pr: _counts_hist) {
        std::cout << pr.first << std::endl;
        for (auto const &x: pr.second) {
            std::cout << x << " ";
        }
        std::cout << "" << std::endl;
    }
}

};
