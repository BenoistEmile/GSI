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