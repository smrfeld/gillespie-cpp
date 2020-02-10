#include "../include/gilsp_bits/counts_hist.hpp"

#include <iostream>
#include <fstream>

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

void CountsHist::store_counts(double time, Counts counts) {
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

void CountsHist::write_count_hist(std::string species, std::string fname) const {
    // ...
}

void CountsHist::write_all_count_hist(std::string dir_name) const {
    // ..
}

void CountsHist::read_all_count_hists(std::string dir_name) {
    //..
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
