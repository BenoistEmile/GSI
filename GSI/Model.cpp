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
    std::unordered_map<std::size_t, unsigned int> expected_abundances; // Expected abundance for each peptide.
    float total_detect;
    float protein_detect;
    for (auto& identification: solution.identifications) { // Filling of expected_abundances
        auto pos = expected_abundances.find(identification->peptide);
        if (pos == expected_abundances.end()) {
            expected_abundances[identification->peptide] = 1;
        }
        else {
            expected_abundances[identification->peptide]++;
        }
    }
    for (auto& expected_abundance: expected_abundances) { // Repartition of expected abundance between identified proteins
        total_detect = 0;
        for (auto& protein: this->Get_Peptide(expected_abundance.first).Get_Proteins()) { // Computation of peptide total detectability
            if (solution.abundances.find(protein.first) == solution.abundances.end()) {
                continue;
            }
            for (auto& iter: protein.second) {
                total_detect += iter;
            }
        }
        for (auto& protein: this->Get_Peptide(expected_abundance.first).Get_Proteins()) { // Computation of raw corrected abundance for each identified protein
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
    for (auto& corr_abundance: solution.corr_abundances) { // Computation of 2 corrected abundances (/N_pep and /Total_detect)
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
    for (auto& protein : solution.abundances) {
        auto pos = solution.corr_abundances.find(protein.first);
        if (pos == solution.corr_abundances.end()) {
            std::cout << "Protein " << protein.first << " not corrected.";
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

void Model::Load_test() {
    proteins.push_back(new Protein(0, "AAAAAAAAAAAAAAAAAAAAA", "P1"));
    proteins.push_back(new Protein(1, "AAAAAAAAAAAAAAAAAAAAA", "P2"));
    proteins.push_back(new Protein(2, "AAAAAAAAAAAAAAAAAAAAA", "P3"));
    proteins.push_back(new Protein(3, "AAAAAAAAAAAAAAAAAAAAA", "P4"));
    peptides.push_back(new Peptide(0, "BBBBBBBBBBBB"));
    peptides.back()->Add_Protein(0);
    proteins[0]->Add_Peptide(0);
    peptides.back()->Add_Protein(1);
    proteins[1]->Add_Peptide(0);
    peptides.back()->Add_Protein(2);
    proteins[2]->Add_Peptide(0);
    peptides.back()->Define_Probabilities(1, 0, 0);
    peptides.back()->Define_Probabilities(0.5, 1, 0);
    peptides.back()->Define_Probabilities(1, 2, 0);
    peptides.push_back(new Peptide(1, "BBBBBBBBBBBB"));
    peptides.back()->Add_Protein(2);
    proteins[2]->Add_Peptide(1);
    peptides.back()->Add_Protein(3);
    proteins[3]->Add_Peptide(1);
    peptides.back()->Define_Probabilities(1, 2, 0);
    peptides.back()->Define_Probabilities(0.5, 3, 0);
    std::vector<Pic*>* pics = new std::vector<Pic*>;
    spectra.push_back(new Spectrum(0, pics));
    spectra.push_back(new Spectrum(1, pics));
    spectra.push_back(new Spectrum(2, pics));
    spectra.push_back(new Spectrum(3, pics));
    spectra.push_back(new Spectrum(4, pics));
    spectra.push_back(new Spectrum(5, pics));
    scores.push_back(new Score(0, 0, 0));
    scores.push_back(new Score(0, 1, 0));
    scores.push_back(new Score(0, 2, 0));
    scores.push_back(new Score(0, 3, 0));
    scores.push_back(new Score(1, 4, 0));
    scores.push_back(new Score(1, 5, 0));
}