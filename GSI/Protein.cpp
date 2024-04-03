#include "Protein.h"


//__________________________________________________________________________________________________________
Protein::Protein(const std::size_t id, const std::string sequence) : id(id), sequence(sequence), peptides({}) {}

Protein::~Protein() {}

//__________________________________________________________________________________________________________

std::ostream& operator<<(std::ostream& os, const Protein& protein) {
    os << "Protein " << protein.id << " : " << protein.sequence << "\n";
    os << "    peptides : [";
    if (protein.peptides.size() > 0) {
        for (auto iter = protein.peptides.begin(); iter != protein.peptides.end() - 1; ++iter) {
            os << (*iter) << " ,";
        }
        os << protein.peptides.back();
    }
    os << "]";
    return os;
}

//__________________________________________________________________________________________________________

const std::size_t Protein::Get_Id() const {
    return id;
}

const std::string& Protein::Get_Sequence() const {
    return sequence;
}

const std::vector<std::size_t>& Protein::Get_Peptides() const {
    return peptides;
}

//__________________________________________________________________________________________________________

bool Protein::Is_Digested() const {
    return false;
}

//__________________________________________________________________________________________________________

void Protein::Add_Peptide(std::size_t peptide) {
    peptides.push_back(peptide);
}