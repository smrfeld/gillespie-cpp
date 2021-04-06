#include "../include/gilsp_bits/counts.hpp"

#include "assert.h"
#include <iostream>
#include <fstream>

namespace gilsp {

// ***************
// MARK: - Constructor
// ***************

Counts::Counts() {
};
Counts::Counts(const Counts& other) {
    _copy(other);
};
Counts::Counts(Counts&& other) {
    _move(other);
};
Counts& Counts::operator=(const Counts& other) {
    if (this != &other) {
        _clean_up();
        _copy(other);
    };
    return *this;
};
Counts& Counts::operator=(Counts&& other) {
    if (this != &other) {
        _clean_up();
        _move(other);
    };
    return *this;
};
Counts::~Counts()
{
    _clean_up();
};
void Counts::_clean_up() {
    // Nothing....
};
void Counts::_copy(const Counts& other) {
    _counts = other._counts;
};
void Counts::_move(Counts& other) {
    _counts = other._counts;
    
    // Reset the other
    other._counts.clear();
};
// ***************
// MARK: - Get/set
// ***************

void Counts::set_count(std::string species, int count) {
    _counts[species] = count;
    
    assert (_counts[species] >= 0);
}

void Counts::increment_count(std::string species, int increment, bool require_non_negative) {
    auto it = _counts.find(species);
    if (it != _counts.end()) {
        it->second += increment;
    } else {
        _counts[species] = increment;
    }
    
    if (require_non_negative) {
        assert (_counts[species] >= 0);
    } else {
        _counts[species] = std::max(0,_counts.at(species));
    }
}

void Counts::increment_count(const std::map<std::string,int> &increments, bool require_non_negative) {
    for (auto pr: increments) {
        auto it = _counts.find(pr.first);
        if (it != _counts.end()) {
            it->second += pr.second;
        } else {
            _counts[pr.first] = pr.second;
        }
        
        if (require_non_negative) {
            assert (_counts[pr.first] >= 0);
        } else {
            _counts[pr.first] = std::max(0,_counts.at(pr.first));
        }
    }
}

int Counts::get_count(std::string species) const {
    return _counts.at(species);
    /*
    auto it = _counts.find(species);
    if (it != _counts.end()) {
        return it->second;
    } else {
        return 0;
    }
     */
}

std::vector<std::string> Counts::get_species() const {
    std::vector<std::string> species;
    for (auto const &pr: _counts) {
        species.push_back(pr.first);
    }
    
    return species;
}

// ***************
// MARK: - Print
// ***************

void Counts::print() const {
    for (auto const &pr: _counts) {
        std::cout << pr.first << ": " << pr.second << std::endl;
    }
}

};
