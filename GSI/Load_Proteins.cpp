#include "Model.h"

#include <fstream>
#include <filesystem>

//__________________________________________________________________________________________________________

void Model::Load_Proteins(const std::string file_name) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/proteins/" + file_name);
    if (file) {
        std::string line;
        std::string sequence = "";
        while (getline(file, line)) {
            if (line[0] == '>') {
                if (Is_Valid_Sequence(sequence)) {
                    proteins.push_back(new Protein(proteins.size(), sequence));
                }
                sequence = "";
            }
            else {
                sequence += line;
            }
        }
        if (Is_Valid_Sequence(sequence)) {
            proteins.push_back(new Protein(proteins.size(), sequence));
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}

void Model::Load_Proteins(const std::string file_name, std::vector<Protein*>(parser)(std::ifstream& file)) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/proteins/" + file_name);
    if (file) {
        std::vector<Protein*> result = parser(file);
        proteins.insert(proteins.end(), result.begin(), result.end());
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}

void Model::Load_Proteins_Accession(const std::string file_name, std::ofstream& output_file) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/proteins/" + file_name);
    if (file) {
        std::string line, accession;
        std::string sequence = "";
        std::size_t index1, index2;
        while (getline(file, line)) {
            if (line[0] == '>') {
                if (Is_Valid_Sequence(sequence)) {
                    proteins.push_back(new Protein(proteins.size(), sequence, accession));
                }
                sequence = "";
                if (line.substr(0, 4) == ">sp|") {
                    index1 = 4;
                }
                else {
                    index1 = 1;
                }
                index2 = line.find("|", index1);
                if (index2 != line.npos) {
                    accession = line.substr(index1, index2 - index1);
                }
                else {
                    accession = "NA";
                }
            }
            else {
                sequence += line;
            }
        }
        if (Is_Valid_Sequence(sequence)) {
            proteins.push_back(new Protein(proteins.size(), sequence, accession));
        }
        output_file << "Loaded " << proteins.size() << " proteins from file " << file_name << std::endl << std::endl;
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}