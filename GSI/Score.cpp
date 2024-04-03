
#include "Score.h"

//__________________________________________________________________________________________________________

Identification::Identification(const std::size_t peptide, const std::size_t spectrum) : peptide(peptide) ,spectrum(spectrum) {}

std::ostream& operator<<(std::ostream& os, const Identification& identification)
{
    os << "Identification (" << identification.peptide << " ," << identification.spectrum << ")";
    return os;
}

//__________________________________________________________________________________________________________

Score::Score(const std::size_t peptide, const std::size_t spectrum, double score) : Identification(peptide, spectrum) , score(score) {}

const Identification* Score::Get_Edge() const {
    return new Identification(peptide, spectrum);
}

std::ostream& operator<<(std::ostream& os, const Score& score)
{
    os << "Score (" << score.peptide << " ," << score.spectrum << ") = " << score.score;
    return os;
}