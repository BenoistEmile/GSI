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
        std::string word, line, sequence;
        int index_spectrum, index_peptide, index_shared_masses, shared_masses;
        unsigned int scores_sum, count;
        std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, int>>*> spectra_scores;
        std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, int>>*>::const_iterator spectrum_scores;
        std::string spectrumID;
        std::size_t spectrum_id;
        while (getline(file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ';')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_spectrum = std::find(row.begin(), row.end(), "spectrum index") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "sequence after specfit") - row.begin();
                index_shared_masses = std::find(row.begin(), row.end(), "score after specfit") - row.begin();
                continue;
            }
            sequence = row[index_peptide];
            for (auto iter = sequence.begin(); iter != sequence.end(); iter++) {
                if (*iter == '(') {
                    std::string::iterator iter2 = iter;
                    while (*iter2 != ')') {
                        iter2++;
                    }
                    iter2++;
                    sequence.erase(iter, iter2);
                }
            }
            spectrum_id = std::stoi(row[index_spectrum]);
            shared_masses = std::stoi(row[index_shared_masses]);
            for (auto& iter : peptides) {
                if (iter->Get_Sequence() == sequence) {
                    spectrum_scores = spectra_scores.find(spectrum_id);
                    if (spectrum_scores == spectra_scores.end()) {
                        spectra_scores[spectrum_id] = new std::vector<std::tuple<std::size_t, int>>;
                        spectra_scores.at(spectrum_id)->push_back({iter->Get_Id(), shared_masses});
                    }
                    else {
                        spectrum_scores->second->push_back({iter->Get_Id(), shared_masses});
                    }
                    break;
                    // if (spectrumID != row[index_spectrum]) {
                    //     spectrumID = row[index_spectrum];
                    //     spectra.push_back(new Spectrum(spectra.size(), new std::vector<Pic*>));
                    //     spectra_scores.push_back(new std::vector<std::tuple<std::size_t, int>>);
                    // }
                    // spectra_scores.at(spectra_scores.size() - 1)->push_back({iter->Get_Id(), std::stoi(row[index_shared_masses])});
                }
            }
        }
        for (auto& spectrum_scores : spectra_scores) {
            scores_sum = 0;
            for (auto psm = spectrum_scores.second->begin(); psm != spectrum_scores.second->end(); psm++) {
                scores_sum += std::get<1>(*psm);
            }
            for (auto psm = spectrum_scores.second->begin(); psm != spectrum_scores.second->end(); psm++) {
                scores.push_back(new Score(std::get<0>(*psm), spectrum_scores.first, 1.0 - ((double)std::get<1>(*psm) / scores_sum)));
            }
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}