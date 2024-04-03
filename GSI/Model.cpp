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

const Peptide& Model::Get_Peptide(std::size_t peptide) const {
    return (*peptides[peptide]);
}

const Spectrum& Model::Get_Spectrum(std::size_t spectrum) const {
    return (*spectra[spectrum]);
}

const Score& Model::Get_Score(std::size_t score) const {
    return (*scores[score]);
}