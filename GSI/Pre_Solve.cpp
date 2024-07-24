#include "Model.h"

void Model::Pre_Solve(const float Pmin) {
    std::set<std::size_t> peptides_with_spectra;

    for (auto& score: scores) {
        if (not peptides_with_spectra.contains(score->peptide)) {
            peptides_with_spectra.emplace(score->peptide);
        }
    }

    for (auto& peptide: peptides) {
        if (peptides_with_spectra.contains(peptide->Get_Id())) {
            continue;
        }
        for (auto& protein: peptide->Get_Proteins()) {
            if ((this->Get_Protein(protein.first).Get_Removed_Edges() + 1) / this->Get_Protein(protein.first).Get_Peptides().size() <= Pmin) {
                
            }
        }
    }
}