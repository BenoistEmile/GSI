#pragma once

#include "Spectrum.h"

#include <unordered_map>

//__________________________________________________________________________________________________________

class Peptide {
private:
    const std::size_t id;
    const std::string sequence;
    std::vector<Pic*>* pics;
    std::unordered_map<std::size_t, std::vector<float>> proteins;
    std::vector<std::size_t> deactivated_proteins;
public:
    //constr/destr
    Peptide(const std::size_t id, const std::string sequence);
    ~Peptide();

    //print
    friend std::ostream& operator<<(std::ostream& os, const Peptide& peptide);

    // getters
    const std::size_t Get_Id() const;
    const std::string Get_Sequence() const;
    const std::unordered_map<std::size_t, std::vector<float>>& Get_Proteins() const;
    const std::vector<Pic*>* Get_Pics() const;
    const bool Get_Protein_Activation(std::size_t protein) const;

    // others
    void Add_Protein(const std::size_t protein);
    void Remove_Protein(const std::size_t protein);
    void Define_Probabilities(float value);
    void Define_Probabilities(float value, std::size_t protein);
    void Define_Probabilities(float value, std::size_t protein, std::size_t position);
    void Build_Spectrum(std::unordered_map<char, double> modifications);
};

