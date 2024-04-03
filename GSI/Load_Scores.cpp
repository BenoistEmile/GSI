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