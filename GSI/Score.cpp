
#include "Score.h"

//__________________________________________________________________________________________________________

Identification::Identification(const std::size_t spectrum) : peptides(), spectrum(spectrum) {}
Identification::Identification(const std::vector<std::size_t>, const std::size_t spectrum) : peptides(peptides), spectrum(spectrum) {}

std::ostream& operator<<(std::ostream& os, const Identification& identification)
{
    os << "Identification (";
    for (auto peptide: identification.peptides) {
        os << peptide << ", ";
    }
    os << identification.spectrum << ")";
    return os;
}

void Identification::Add_Peptide(std::size_t peptide) {
    peptides.push_back(peptide);
}

//__________________________________________________________________________________________________________

Score::Score(const std::vector<std::size_t> peptides, const std::size_t spectrum, double score) : Identification(peptides, spectrum) , score(score) {}

const Identification* Score::Get_Edge() const {
    return new Identification(peptides, spectrum);
}

std::ostream& operator<<(std::ostream& os, const Score& score)
{
    os << "Score (";
    for (auto& peptide: score.peptides) {
        os << peptide << ", ";
    }
    os << score.spectrum << ") = " << score.score;
    return os;
}