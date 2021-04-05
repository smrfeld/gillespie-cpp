#include <string>
#include <map>
#include <vector>

#ifndef REACTION_H
#define REACTION_H

namespace gilsp {

    // ***************
    // MARK: - Reaction
    // ***************
    
    class Rxn {
        
    private:
                
        // Name
        std::string _name;
        
        // Rate
        double _kr;
        
        // Products and reactants
        std::vector<std::string> _r_list;
        std::vector<std::string> _p_list;

        // Species
        std::vector<std::string> _species;
        
        // Multiplicities
        std::vector<std::pair<std::string,int>> _species_r_multiplicity;
        std::vector<std::pair<std::string,int>> _species_p_multiplicity;

        // Internal update species and r/p multiplicity
        void _update_species_and_multiplicities();
        
        // Internal copy func/clean up
        void _clean_up();
        void _copy(const Rxn& other);
        void _move(Rxn &other);
        
    public:
                
        
        // ***************
        // MARK: - Constructor
        // ***************
        
        Rxn(std::string name, double kr, const std::vector<std::string> reactants, const std::vector<std::string> products);
        Rxn(const Rxn& other);
        Rxn& operator=(const Rxn& other);
        Rxn(Rxn&& other);
        Rxn& operator=(Rxn&& other);
        ~Rxn();
        
        // ***************
        // MARK: - Getters/setters
        // ***************
        
        std::string get_name() const;
        void set_name(std::string name);
        
        double get_kr() const;
        void set_kr(double kr);
        
        const std::vector<std::string>& get_reactants() const;
        void set_reactants(std::vector<std::string> reactants);
        int get_no_reactants() const;
        
        const std::vector<std::string>& get_products() const;
        void set_products(std::vector<std::string> products);
        int get_no_products() const;
        
        const std::vector<std::pair<std::string,int>>& get_reactant_multiplicity() const;
        const std::vector<std::pair<std::string,int>>& get_product_multiplicity() const;
        
        const std::vector<std::string>& get_species() const;
                
        // ***************
        // MARK: - Print
        // ***************
        
        void print() const;
    };
};

#endif
