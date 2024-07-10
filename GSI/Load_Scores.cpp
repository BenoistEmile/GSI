#include "Model.h"

#include <fstream>
#include <filesystem>

bool Same_Peptide(const std::string seq_1, const std::string seq_2, const int min_size, const int max_size, const int max_len_diff) {
    if (seq_1.size() < min_size || seq_2.size() < min_size || seq_1.size() > max_size || seq_2.size() > max_size) {
        return false;
    }
    // else if (abs(seq_1.size() - seq_2.size()) > max_len_diff) {
    //     return false;
    // }
    else if (seq_1.find(seq_2) != seq_1.npos || seq_2.find(seq_1) != seq_2.npos) {
        return true;
    }
    else {
        return false;
    }
}

int LevensteinDistance(const std::string& seq_1, const std::string& seq_2) {
    int m = seq_1.length();
    int n = seq_2.length();
    std::vector<std::vector<int>> D(m+1, std::vector<int>(n+1, 0));

    // Initialisation de la matrice D
    for (int i = 0; i <= m; i++) {
        D[i][0] = i;
    }
    for (int j = 0; j <= n; j++) {
        D[0][j] = j;
    }

    // Cout de substitution
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (seq_1[i - 1] != seq_2[j - 1]) {
                D[i][j] = 1;
            }
        }
    }

    for (int i =1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            D[i][j] = std::min(D[i-1][j] + 1, std::min(D[i][j-1] + 1, D[i-1][j-1] + D[i][j]));
        }
    }
    return D[m][n];
}

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
        unsigned int scores_sum, count, len_to_remove;
        std::unordered_map<std::size_t, std::unordered_map<std::size_t, int>*> spectra_scores;
        std::unordered_map<std::size_t, std::unordered_map<std::size_t, int>*>::const_iterator spectrum_scores;
        std::unordered_map<std::size_t, int>::iterator same_peptide;
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
                if (index_peptide == row.size()) {
                    index_peptide = std::find(row.begin(), row.end(), "peptide") - row.begin();
                    if (index_peptide == row.size()) {
                        throw "Can't find peptide sequences";
                    }
                }
                index_shared_masses = std::find(row.begin(), row.end(), "score after specfit") - row.begin();
                if (index_shared_masses == row.size()) {
                    index_shared_masses = std::find(row.begin(), row.end(), "sharedMassesAfterAlign") - row.begin();
                    if (index_peptide == row.size()) {
                        throw "Can't find shared masses";
                    }
                }
                continue;
            }
            sequence = row[index_peptide];
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence[i] == '(') {
                    int len_to_remove = 1;
                    int begin_remove = i;
                    while (sequence[i] != ')') {
                        i++;
                        len_to_remove++;
                    }
                    sequence.erase(begin_remove, len_to_remove);
                    i -= len_to_remove;
                }
            }
            spectrum_id = std::stoi(row[index_spectrum]);
            shared_masses = std::stoi(row[index_shared_masses]);
            bool found_peptide = false;
            for (auto& iter : peptides) {
                // if (iter->Get_Sequence() == sequence) {
                if (Same_Peptide(iter->Get_Sequence(), sequence, 7, 25, 3)) {
                    found_peptide = true;
                    spectrum_scores = spectra_scores.find(spectrum_id);
                    if (spectrum_scores == spectra_scores.end()) {
                        spectra_scores[spectrum_id] = new std::unordered_map<std::size_t, int>;
                        spectra_scores.at(spectrum_id)->emplace(iter->Get_Id(), shared_masses);
                    }
                    else {
                        same_peptide = spectrum_scores->second->find(iter->Get_Id());
                        if (same_peptide == spectrum_scores->second->end()) {
                            spectrum_scores->second->emplace(iter->Get_Id(), shared_masses);
                        }
                        else if (shared_masses > same_peptide->second) {
                            same_peptide->second = shared_masses;
                        }
                    }
                    break;
                }
            }
            // if (not found_peptide) {
            //     std::cout << spectrum_id << ", " << sequence << std::endl;
            // }
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

void Model::Load_Scores_Prospect(const std::string file_name, const int min_length, const int max_length, const int min_pics) {
    std::ifstream file(std::filesystem::current_path() / "data" / "scores" / file_name);
    if (file) {
        bool first_line = true;
        std::vector<std::string> row;
        std::string word, line, sequence;
        int index_spectrum, index_peptide, index_shared_masses, shared_masses;
        unsigned int scores_sum, count;
        std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, int>>*> spectra_scores;
        std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, int>>*>::const_iterator spectrum_scores;
        std::size_t spectrum_id;
        while (getline(file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ';')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_spectrum = std::find(row.begin(), row.end(), "spectrum id") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "preliminary: peptide") - row.begin();
                index_shared_masses = std::find(row.begin(), row.end(), "preliminary: number of common peaks") - row.begin();
                continue;
            }
            sequence = row[index_peptide];
            spectrum_id = std::stoi(row[index_spectrum]);
            shared_masses = std::stoi(row[index_shared_masses]);
            if (sequence.size() < min_length || sequence.size() > max_length || shared_masses < min_pics) {
                continue;
            }

            for (auto& iter : peptides) {
                // if (iter->Get_Sequence() == sequence) {
                if (Same_Peptide(iter->Get_Sequence(), sequence, 6, 25, 3)) {
                // if (LevensteinDistance(iter->Get_Sequence(), sequence) <= 2) {
                    spectrum_scores = spectra_scores.find(spectrum_id);
                    // std::cout << "Peptide : " << iter->Get_Sequence() << ", sequence : " << sequence << ", score : " << LevensteinDistance(iter->Get_Sequence(), sequence) << std::endl;
                    if (spectrum_scores == spectra_scores.end()) {
                        spectra_scores[spectrum_id] = new std::vector<std::tuple<std::size_t, int>>;
                        spectra_scores.at(spectrum_id)->push_back({iter->Get_Id(), shared_masses});
                        std::cout << iter->Get_Sequence() << ", " << sequence << ", " << shared_masses << std::endl;
                    }
                    else {
                        spectrum_scores->second->push_back({iter->Get_Id(), shared_masses});
                    }
                    break;
                }
            }

            // int score, best_score, best_peptide;
            // for (auto& iter : peptides) {
            //     score = LevensteinDistance(iter->Get_Sequence(), sequence);
            //     if (score < best_score) {
            //         best_score = score;
            //         best_peptide = iter->Get_Id();
            //     }
            // }
            // spectra_scores.find(spectrum_id);
            // std::cout << "Peptide : " << iter->Get_Sequence() << ", sequence : " << sequence << ", score : " << LevensteinDistance(iter->Get_Sequence(), sequence) << std::endl;
            // if (spectrum_scores == spectra_scores.end()) {
            //     spectra_scores[spectrum_id] = new std::vector<std::tuple<std::size_t, int>>;
            //     spectra_scores.at(spectrum_id)->push_back({iter->Get_Id(), shared_masses});
            // }
            // else {
            //     spectrum_scores->second->push_back({iter->Get_Id(), shared_masses});
            // }
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

void Model::Load_Scores_XTandem(const std::string file_name) {
    std::ifstream file(std::filesystem::current_path() / "data" / "scores" / file_name);
    if (file) {
        bool first_line = true;
        std::vector<std::string> row;
        std::string word, line, sequence;
        int index_spectrum, index_peptide, index_shared_masses, shared_masses;
        unsigned int scores_sum, count, len_to_remove;
        std::unordered_map<std::size_t, std::unordered_map<std::size_t, int>*> spectra_scores;
        std::unordered_map<std::size_t, std::unordered_map<std::size_t, int>*>::const_iterator spectrum_scores;
        std::unordered_map<std::size_t, int>::iterator same_peptide;
        std::size_t spectrum_id;
        while (getline(file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ';')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_spectrum = std::find(row.begin(), row.end(), "Scan") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "Sequence") - row.begin();
                index_shared_masses = std::find(row.begin(), row.end(), "hyperscore") - row.begin();
                continue;
            }
            sequence = row[index_peptide];
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence[i] == '(') {
                    int len_to_remove = 1;
                    int begin_remove = i;
                    while (sequence[i] != ')') {
                        i++;
                        len_to_remove++;
                    }
                    sequence.erase(begin_remove, len_to_remove);
                    i -= len_to_remove;
                }
            }
            spectrum_id = std::stoi(row[index_spectrum]);
            shared_masses = std::stoi(row[index_shared_masses]);
            bool found_peptide = false;
            for (auto& iter : peptides) {
                // if (iter->Get_Sequence() == sequence) {
                if (Same_Peptide(iter->Get_Sequence(), sequence, 0, 50, 3)) {
                    found_peptide = true;
                    spectrum_scores = spectra_scores.find(spectrum_id);
                    if (spectrum_scores == spectra_scores.end()) {
                        spectra_scores[spectrum_id] = new std::unordered_map<std::size_t, int>;
                        spectra_scores.at(spectrum_id)->emplace(iter->Get_Id(), shared_masses);
                    }
                    else {
                        same_peptide = spectrum_scores->second->find(iter->Get_Id());
                        if (same_peptide == spectrum_scores->second->end()) {
                            spectrum_scores->second->emplace(iter->Get_Id(), shared_masses);
                        }
                        else if (shared_masses > same_peptide->second) {
                            same_peptide->second = shared_masses;
                        }
                    }
                    break;
                }
            }
            if (not found_peptide) {
                std::cout << spectrum_id << ", " << sequence << std::endl;
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