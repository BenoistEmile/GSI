#include "Model.h"

#include <fstream>
#include <filesystem>

void Model::Load_Scores(const std::string file_name, std::vector<Score*>(parser)(std::ifstream& file)) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/scores/" + file_name);
    if (file) {
        std::vector<Score*> result = parser(file);
        scores.insert(scores.end(), result.begin(), result.end());
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}

void Model::Load_Scores_SpecOMS(const std::string file_name) {
    std::ifstream file(std::filesystem::current_path() / "data" / "scores" / file_name);
    if (file) {
        bool first_line = true;
        std::vector<std::string> row;
        std::string word, line;
        int index_spectrum, index_peptide, index_shared_masses;
        unsigned int scores_sum;
        std::vector<std::vector<std::tuple<std::size_t, int>>*> spectra_scores;
        std::string spectrumID;
        while (getline(file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ';')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_spectrum = std::find(row.begin(), row.end(), "spectrumID") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "peptide") - row.begin();
                index_shared_masses = std::find(row.begin(), row.end(), "sharedMassesAfterSpecFit") - row.begin();
                continue;
            }
            for (auto& iter : peptides) {
                if (iter->Get_Sequence() == row[index_peptide]) {
                    if (spectrumID != row[index_spectrum]) {
                        spectrumID = row[index_spectrum];
                        spectra.push_back(new Spectrum(spectra.size(), new std::vector<Pic*>));
                        spectra_scores.push_back(new std::vector<std::tuple<std::size_t, int>>);
                    }
                    spectra_scores.at(spectra_scores.size() - 1)->push_back({iter->Get_Id(), std::stoi(row[index_shared_masses])});
                }
            }
        }
        for (int i = 0; i < spectra_scores.size(); i++) {
            scores_sum = 0;
            for (auto& psm : *spectra_scores.at(i)) {
                scores_sum += std::get<1>(psm);
            }
            for (auto& psm : *spectra_scores.at(i)) {
                scores.push_back(new Score(std::get<0>(psm), i, 1.0 - ((double)std::get<1>(psm) / scores_sum)));
            }
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}