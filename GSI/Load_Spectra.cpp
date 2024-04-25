#include "Model.h"

#include <fstream>
#include <filesystem>

//__________________________________________________________________________________________________________

void Model::Load_Spectra(const std::string file_name) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/spectra/" + file_name);
    if (file) {
        std::string line;
        std::vector<Pic*>* pics = new std::vector<Pic*>;
        std::size_t space_position;
        double mass;
        double intensity;
        while (getline(file, line)) {
            if (line[0] < '0' || line[0] > '9') {
                if ((*pics).size() != 0) {
                    try {
                        spectra.push_back(new Spectrum(spectra.size(), pics));
                    }
                    catch (std::string& str_err) {
                        std::cout << str_err << std::endl;
                        std::cout << "Skipping current spectrum" << std::endl;
                        pics = new std::vector<Pic*>;
                        continue;
                    }
                    pics = new std::vector<Pic*>;
                }
            }
            else {
                mass = std::stod(line, &space_position);
                intensity = std::stod(line.substr(space_position));
                pics->push_back(new Pic(mass ,intensity));
            }
        }
        if ((*pics).size() != 0) {
            try {
                spectra.push_back(new Spectrum(spectra.size(), pics));
            }
            catch (std::string& str_err) {
                std::cout << str_err << std::endl;
                std::cout << "Skipping current spectrum" << std::endl;
                pics = new std::vector<Pic*>;
            }
        }
        else {
            delete pics;
        }
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}

void Model::Load_Spectra(const std::string file_name, std::vector<Spectrum*>(parser)(std::ifstream& file)) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/spectra/" + file_name);
    if (file) {
        std::vector<Spectrum*> result = parser(file);
        spectra.insert(spectra.end(), result.begin(), result.end());
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}

void Model::Load_Spectra(const std::string file_name, std::ofstream& output_file) {
    std::ifstream file(std::filesystem::current_path().generic_string() + "/data/spectra/" + file_name);
    if (file) {
        std::string line;
        std::vector<Pic*>* pics = new std::vector<Pic*>;
        std::size_t space_position;
        double mass;
        double intensity;
        while (getline(file, line)) {
            if (line[0] < '0' || line[0] > '9') {
                if ((*pics).size() != 0) {
                    try {
                        spectra.push_back(new Spectrum(spectra.size(), pics));
                    }
                    catch (std::string& str_err) {
                        std::cout << str_err << std::endl;
                        std::cout << "Skipping current spectrum" << std::endl;
                        pics = new std::vector<Pic*>;
                        continue;
                    }
                    pics = new std::vector<Pic*>;
                }
            }
            else {
                try {
                    mass = std::stod(line, &space_position);
                    intensity = std::stod(line.substr(space_position));
                    pics->push_back(new Pic(mass ,intensity));
                }
                catch (std::invalid_argument& error) {
                    std::cout << error.what() << std::endl;
                    std::cout << "Invalid pic, skipping to next pic" << std::endl;
                }
            }
        }
        if ((*pics).size() != 0) {
            try {
                spectra.push_back(new Spectrum(spectra.size(), pics));
            }
            catch (std::string& str_err) {
                std::cout << str_err << std::endl;
                std::cout << "Skipping current spectrum" << std::endl;
                pics = new std::vector<Pic*>;
            }
        }
        else {
            delete pics;
        }
        output_file << "Loaded " << spectra.size() << " spectra from file " << file_name << std::endl << std::endl;
    }
    else {
        std::cout << "ERROR : Impossible to open the file named : " << file_name << std::endl;
    }
}