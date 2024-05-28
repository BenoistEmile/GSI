
#include "Peptide.h"
#include "Amino_Acids.h"
#include "Fonction.h"

//__________________________________________________________________________________________________________

Peptide::Peptide(const std::size_t id, const std::string sequence) : id(id), sequence(sequence) {
    proteins = {};
    pics = new std::vector<Pic*>;
}

Peptide::~Peptide() {
    for (Pic* pic : *pics) {
        delete pic;
    }
    delete pics;
}

//__________________________________________________________________________________________________________

const std::size_t Peptide::Get_Id() const {
    return id;
}

const std::string Peptide::Get_Sequence() const {
    return sequence;
}

const std::unordered_map<std::size_t, std::vector<float>>& Peptide::Get_Proteins() const {
    return proteins;
}

const std::vector<Pic*>* Peptide::Get_Pics() const {
    return pics;
}

//__________________________________________________________________________________________________________

std::ostream& operator<<(std::ostream& os, const Peptide& peptide)
{
    os << "Peptide " << peptide.id << " : " << peptide.sequence << "\n";
    for (auto iter = peptide.proteins.begin(); iter != peptide.proteins.end(); ++iter) {
        os << "   protein " << iter->first << " : (";
        for (std::size_t i = 0; i < iter->second.size(); ++i) {
            os << iter->second[i] << " ,";
        }
        os << ")\n";
    }
    return os;
}

//__________________________________________________________________________________________________________

void Peptide::Add_Protein(const std::size_t protein) {
    if (this->proteins.find(protein) == proteins.end()) {
        proteins[protein] = { 1.0f };
    }
    else {
        proteins[protein].push_back(1.0f);
    }
}

//_________________________

void Peptide::Define_Probabilities(float value) {
    for (auto iter = proteins.begin(); iter != proteins.end(); ++iter) {
        proteins[iter->first] = std::vector<float>(iter->second.size(), value);
    }
}

void Peptide::Define_Probabilities(float value, std::size_t protein) {
    proteins[protein] = std::vector<float>(proteins[protein].size(), value);
}

void Peptide::Define_Probabilities(float value, std::size_t protein, std::size_t position) {
    proteins[protein][position] = value;
}

//_________________________

void Peptide::Build_Spectrum(std::unordered_map<char, double> modifications) {
    //std::vector<Pic*>* pics = new std::vector<Pic*>;
    for (Pic* pic : *pics) {
        delete pic;
    }
    pics->clear();
    double sumMass = AminoAcid::Get_Y_Mass();
    for (std::size_t i = 1; i <= sequence.size(); ++i) {
        //sumMass = Round_Precision(Amino_Acids.at(sequence[sequence.size() - i]).Get_Mono_Isotopic_Mass() + sumMass, 3);
        sumMass += Amino_Acids.at(sequence[sequence.size() - i]).Get_Mono_Isotopic_Mass();
        if (modifications.contains(Amino_Acids.at(sequence[sequence.size() - i]).Get_Letter())) {
            sumMass += modifications[Amino_Acids.at(sequence[sequence.size() - i]).Get_Letter()];
        }
        (*pics).push_back(new Pic(sumMass, 100.f));
    }
    sumMass = AminoAcid::Get_B_Mass();
    auto j = (*pics).begin();
    for (std::size_t i = 0; i < sequence.size() - 1; ++i) {
        //sumMass = Round_Precision(Amino_Acids.at(sequence[i]).Get_Mono_Isotopic_Mass() + sumMass, 3);
        sumMass += Amino_Acids.at(sequence[i]).Get_Mono_Isotopic_Mass();
        if (modifications.contains(Amino_Acids.at(sequence[i]).Get_Letter())) {
            sumMass += modifications[Amino_Acids.at(sequence[i]).Get_Letter()];
        }
        while ((**j).mass < sumMass) {
            ++j;
        }
        if ((**j).mass != sumMass) {
            j = (*pics).insert(j, new Pic(sumMass, 100.0));
        }
        ++j;
    }
    //return pics;
}