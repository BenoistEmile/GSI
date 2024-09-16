#include "Model.h"

//__________________________________________________________________________________________________________

Model::Model() : proteins({}), peptides({}), peptides_sequences({}) {};

Model::~Model() {
    for (std::vector<Protein*>::iterator iter = proteins.begin(); iter != proteins.end(); ++iter) {
        delete* iter;
    }
    for (std::vector<Peptide*>::iterator iter = peptides.begin(); iter != peptides.end(); ++iter) {
        delete* iter;
    }
    for (std::vector<Spectrum*>::iterator iter = spectra.begin(); iter != spectra.end(); ++iter) {
        delete* iter;
    }
    for (std::vector<Score*>::iterator iter = scores.begin(); iter != scores.end(); ++iter) {
        delete* iter;
    }
    solution.Clear();
}

//__________________________________________________________________________________________________________

void Model::Clear(bool clear_proteins, bool clear_peptides, bool clear_spectra, bool clear_scores, bool clear_solution) {
    if (clear_proteins || clear_peptides || clear_spectra || clear_scores || clear_solution) {
        solution.Clear();
        if (clear_proteins || clear_peptides || clear_spectra || clear_scores) {
            for (std::vector<Score*>::iterator iter = scores.begin(); iter != scores.end(); ++iter) {
                delete* iter;
            }
            scores.clear();
            if (clear_spectra) {
                for (std::vector<Spectrum*>::iterator iter = spectra.begin(); iter != spectra.end(); ++iter) {
                    delete* iter;
                }
                spectra.clear();
            }
            if (clear_proteins || clear_peptides) {
                for (std::vector<Peptide*>::iterator iter = peptides.begin(); iter != peptides.end(); ++iter) {
                    delete* iter;
                }
                peptides.clear();
                peptides_sequences = {};
                if (clear_proteins) {
                    for (std::vector<Protein*>::iterator iter = proteins.begin(); iter != proteins.end(); ++iter) {
                        delete* iter;
                    }
                    proteins.clear();
                }
            }
        }
    }
}

void Model::Correct_Abundances() {
    std::unordered_map<std::size_t, unsigned int> expected_abundances;
    float total_detect;
    float protein_detect;
    for (auto& identification: solution.identifications) {
        auto pos = expected_abundances.find(identification->peptide);
        if (pos == expected_abundances.end()) {
            expected_abundances[identification->peptide] = 1;
        }
        else {
            expected_abundances[identification->peptide]++;
        }
    }
    for (auto& expected_abundance: expected_abundances) {
        total_detect = 0;
        for (auto& protein: this->Get_Peptide(expected_abundance.first).Get_Proteins()) {
            if (solution.abundances.find(protein.first) == solution.abundances.end()) {
                continue;
            }
            for (auto& iter: protein.second) {
                total_detect += iter;
            }
        }
        for (auto& protein: this->Get_Peptide(expected_abundance.first).Get_Proteins()) {
            if (solution.abundances.find(protein.first) == solution.abundances.end()) {
                continue;
            }
            protein_detect = 0;
            for (auto& iter: protein.second) {
                protein_detect += iter;
            }
            auto pos = solution.corr_abundances.find(protein.first);
            if (pos == solution.corr_abundances.end()) {
                solution.corr_abundances[protein.first] = std::tuple<float, float, float>(expected_abundance.second * protein_detect / total_detect, 0, 0);
            }
            else {
                std::get<0>(solution.corr_abundances[protein.first]) += expected_abundance.second * protein_detect / total_detect;
            }
        }
    }
    for (auto& corr_abundance: solution.corr_abundances) {
        total_detect = 0;
        std::get<1>(corr_abundance.second) = std::get<0>(corr_abundance.second);
        std::get<2>(corr_abundance.second) = std::get<0>(corr_abundance.second);
        std::get<1>(corr_abundance.second) /= this->Get_Protein(corr_abundance.first).Get_Peptides().size();
        for (std::size_t pep_id: this->Get_Protein(corr_abundance.first).Get_Peptides()) {
            for (auto& protein: this->Get_Peptide(pep_id).Get_Proteins()) {
                if (corr_abundance.first == protein.first) {
                    for (auto& detect: protein.second) {
                        total_detect += detect;
                    }
                }
            }
        }
        std::get<2>(corr_abundance.second) /= total_detect;
    }
}

//__________________________________________________________________________________________________________

const std::size_t Model::Number_Of_Proteins() const {
    return proteins.size();
}

const std::size_t Model::Number_Of_Peptides() const {
    return peptides.size();
}

const std::size_t Model::Number_Of_Spectra() const {
    return spectra.size();
}

const std::size_t Model::Number_Of_Scores() const {
    return scores.size();
}

const Protein& Model::Get_Protein(std::size_t protein) const {
    return (*proteins[protein]);
}

const Protein& Model::Get_Protein(std::string accession) const {
    return (*proteins[proteins_accession.at(accession)]);
}

const Peptide& Model::Get_Peptide(std::size_t peptide) const {
    return (*peptides[peptide]);
}

const Spectrum& Model::Get_Spectrum(std::size_t spectrum) const {
    return (*spectra[spectrum]);
}

const Score& Model::Get_Score(std::size_t score) const {
    return (*scores[score]);
}