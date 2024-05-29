#include "Model.h"

#include <stdlib.h>
#include <time.h>

//__________________________________________________________________________________________________________

void Model::Define_Probabilities(float value) {
    for (unsigned int i = 0; i < peptides.size(); ++i) {
        (*peptides[i]).Define_Probabilities(value);
    }
}

void Model::Define_Probabilities(bool per_protein) {
    srand((unsigned int)time(NULL));
    if (per_protein) {
        float value;
        for (unsigned int i = 0; i < proteins.size(); ++i) {
            value = rand() / (float)RAND_MAX;
            for (unsigned int j = 0; j < (*proteins[i]).Get_Peptides().size(); ++j) {
                (*(peptides[(*proteins[i]).Get_Peptides()[j]])).Define_Probabilities(value, i);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < peptides.size(); ++i) {
            auto proteins = (*peptides[i]).Get_Proteins();
            for (auto iter = proteins.begin(); iter != proteins.end(); ++iter) {
                for (unsigned int j = 0; j < iter->second.size(); ++j) {
                    (*peptides[i]).Define_Probabilities(rand() / (float)RAND_MAX, iter->first, j);
                }
            }
        }
    }
}

void Model::Define_Probabilities(void (define_probabilities_function)(std::vector<Protein*> proteins, std::vector<Peptide*> peptides)) {
    define_probabilities_function(proteins, peptides);
}

void Model::Define_Probabilities(const std::string file_name) {
    std::ifstream proba_file(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name);
    if (proba_file) {
        bool first_line = true;
        std::size_t previous_protein = 0;
        std::unordered_map<std::size_t, std::size_t> positions;
        std::vector<std::string> row;
        std::string word, line;
        std::unordered_map<std::size_t, std::size_t>::iterator position;
        std::size_t protein, peptide;
        int index_protein, index_peptide, index_proba;
        while (getline(proba_file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ',')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_protein = std::find(row.begin(), row.end(), "protein_id") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "peptide_id") - row.begin();
                index_proba = std::find(row.begin(), row.end(), "Prob") - row.begin();
                continue;
            }
            protein = std::stoi(row[index_protein]);
            peptide = std::stoi(row[index_peptide]);
            if (protein != previous_protein) {
                previous_protein = protein;
                positions.clear();
            }
            position = positions.find(peptide);
            if (position == positions.end()) {
                positions[peptide] = 0;
                (*peptides[peptide]).Define_Probabilities(std::stof(row[index_proba]), protein, 0);
            }
            else {
                position->second++;
                (*peptides[peptide]).Define_Probabilities(std::stof(row[index_proba]), protein, position->second);
            }
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
    proba_file.close();
}

void Model::Define_Probabilities(const std::string file_name, std::ofstream& output_file) {
    std::ifstream proba_file(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name);
    if (proba_file) {
        bool first_line = true;
        std::size_t previous_protein = 0;
        std::unordered_map<std::size_t, std::size_t> positions;
        std::vector<std::string> row;
        std::string word, line;
        std::unordered_map<std::size_t, std::size_t>::iterator position;
        std::size_t protein, peptide;
        int index_protein, index_peptide, index_proba;
        unsigned int count = 0;
        while (getline(proba_file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ',')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_protein = std::find(row.begin(), row.end(), "protein_id") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "peptide_id") - row.begin();
                index_proba = std::find(row.begin(), row.end(), "Prob") - row.begin();
                continue;
            }
            protein = std::stoi(row[index_protein]);
            peptide = std::stoi(row[index_peptide]);
            if (protein != previous_protein) {
                previous_protein = protein;
                positions.clear();
            }
            position = positions.find(peptide);
            if (position == positions.end()) {
                positions[peptide] = 0;
                (*peptides[peptide]).Define_Probabilities(std::stof(row[index_proba]), protein, 0);
                count++;
            }
            else {
                position->second++;
                (*peptides[peptide]).Define_Probabilities(std::stof(row[index_proba]), protein, position->second);
                count++;
            }
        }
        output_file << "Loaded " << count << " scores from file " << file_name << std::endl << std::endl;
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
    proba_file.close();
}

void Model::Define_Probabilities(const std::string file_name, float min_proba) {
    std::ifstream proba_file(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name);
    if (proba_file) {
        bool first_line = true;
        std::size_t previous_protein = 0;
        std::unordered_map<std::size_t, std::size_t> positions;
        std::vector<std::string> row;
        std::string word, line;
        std::unordered_map<std::size_t, std::size_t>::iterator position;
        std::size_t protein, peptide;
        float proba;
        int index_protein, index_peptide, index_proba;
        while (getline(proba_file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ',')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_protein = std::find(row.begin(), row.end(), "protein_id") - row.begin();
                index_peptide = std::find(row.begin(), row.end(), "peptide_id") - row.begin();
                index_proba = std::find(row.begin(), row.end(), "Prob") - row.begin();
                continue;
            }
            protein = std::stoi(row[index_protein]);
            peptide = std::stoi(row[index_peptide]);
            proba = std::stof(row[index_proba]);
            if (proba >= min_proba) {
                if (protein != previous_protein) {
                    previous_protein = protein;
                    positions.clear();
                }
                position = positions.find(peptide);
                if (position == positions.end()) {
                    positions[peptide] = 0;
                    (*peptides[peptide]).Define_Probabilities(proba, protein, 0);
                }
                else {
                    position->second++;
                    (*peptides[peptide]).Define_Probabilities(proba, protein, position->second);
                }
            }
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
    proba_file.close();
}

void Model::Define_Probabilities_2(const std::string file_name, const float min_proba) {
    std::ifstream proba_file(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name);
    if (proba_file) {
        bool first_line = true;
        std::size_t previous_protein = 0;
        std::unordered_map<std::size_t, std::size_t> positions;
        std::vector<std::string> row;
        std::string word, line, sequence;
        std::unordered_map<std::size_t, std::size_t>::iterator position;
        std::size_t protein, peptide;
        float proba;
        int index_protein, index_peptide, index_proba, index_seq;
        while (getline(proba_file, line)) {
            row.clear();
            std::stringstream s(line);
            while (getline(s, word, ',')) {
                row.push_back(word);
            }
            if (first_line) {
                first_line = false;
                index_protein = std::find(row.begin(), row.end(), "protein_id") - row.begin();
                index_proba = std::find(row.begin(), row.end(), "Prob") - row.begin();
                index_seq = std::find(row.begin(), row.end(), "peptide") - row.begin();
                continue;
            }
            protein = std::stoi(row[index_protein]);
            proba = std::stof(row[index_proba]);
            sequence = row[index_seq];
            if (proba >= min_proba) {
                if (protein != previous_protein) {
                    previous_protein = protein;
                    positions.clear();
                }
                if (peptide >= peptides.size()) {
                    peptides.push_back(new Peptide(peptides.size(), sequence));
                    if (peptides.at(peptide)->Get_Id() != peptide) {
                        throw "The last built peptide id does not match the expected id.";
                    }
                }
                peptides.at(peptide)->Add_Protein(protein);
                (*proteins.at(protein)).Add_Peptide(peptide);
                position = positions.find(peptide);
                if (position == positions.end()) {
                    positions[peptide] = 0;
                    (*peptides[peptide]).Define_Probabilities(proba, protein, 0);
                }
                else {
                    position->second++;
                    (*peptides[peptide]).Define_Probabilities(proba, protein, position->second);
                }
            }
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
    proba_file.close();
}