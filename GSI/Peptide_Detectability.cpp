#include "Model.h"
#include <filesystem>

//__________________________________________________________________________________________________________

void Model::AP3_Fasta() {
    std::ofstream temp_fasta(std::filesystem::current_path() / "local" / "temp_fasta.fasta");
    for (auto& protein : proteins) {
        temp_fasta << ">" << protein->Get_Accession() << "|" << std::endl;
        for (unsigned int i = 0; i < protein->Get_Sequence().length(); i += 60) {
            temp_fasta << protein->Get_Sequence().substr(i, 60) << std::endl;
        }
    }
    temp_fasta.close();
}

//__________________________________________________________________________________________________________

void Model::Peptide_Detectability(std::string env_name, std::string digestion_file_name) const {
    std::filesystem::path curr_path = std::filesystem::current_path();
    std::filesystem::path script_path = std::filesystem::path("~/stage/code/DbyDeep-main/dbydeep_model.py");
    std::filesystem::path data_path = curr_path / "data" / "digestion" / digestion_file_name;
    data_path += ".csv";
    std::string command = "conda run -n " + env_name
    + " python " + script_path.generic_string()
    + " --data-path " + data_path.generic_string()
    + " --save-path " + (curr_path / "data" / "digestion").generic_string()
    + " --job-name " + digestion_file_name + "_result"
    + " --model-path " + (script_path.parent_path() / "data" / "DbyDeep.h5").generic_string();
    // std::string command = ("conda run -n " + env_name + " --cwd ~/stage/code/DbyDeep-main pwd");
    // std::cout << "Command : " << command << std::endl;
    int return_code = system(command.c_str());
    if (return_code != 0) {
        std::cout << "An error occured during peptide detectability prediction" << std::endl;
    }
}

void Model::Peptide_Detectability(int detectability_model, float min_detect, int minimum_number_of_amino_acids, int maximum_number_of_amino_acids) {
    int return_code;
    this->In_Silico_Digestion_2("digestion_file", minimum_number_of_amino_acids, maximum_number_of_amino_acids);
    switch (detectability_model) {
        case 1: {
            std::string command = "conda run -n detectability python ./python/DbyDeep_script.py";
            return_code = system(command.c_str());
            if (return_code != 0) {
                throw "An error occured during peptide detectability prediction : " + return_code;
            }
            break;
        }
        case 2: {
            AP3_Fasta();
            std::ofstream param_file(std::filesystem::current_path() / "local" / "ap3_params.param");
            param_file << "InputFilePath=\"" << (std::filesystem::current_path() / "local" / "temp_fasta.fasta").string() << "\"" << std::endl;
            param_file << "Model=\"Saccharomyces cerevisiae\"\nIdentifierParsingRule=\">(.*?)\\|\"\n";
            param_file << "ResultPath=\"" << (std::filesystem::current_path() / "local" / "ap3_results").string() << "\"" << std::endl;
            param_file << "AllowMinPeptideLength=\"" << minimum_number_of_amino_acids << "\"" << std::endl;
            param_file << "AllowMaxPeptideLength=\"" << maximum_number_of_amino_acids << "\"" << std::endl;
            param_file << "AllowMissingCutNumber=\"0\"" << std::endl;
            param_file.close();
            std::string ap3_loc = "AP3.exe";
            std::string command_ap3 = ap3_loc + " " + (std::filesystem::current_path() / "local" / "ap3_params.param").string();
            // return_code = system(command_ap3.c_str());
            // if (return_code != 0) {
            //     throw "An error occured during peptide detectability prediction : " + return_code;
            // }
            std::string command = "conda run -n detectability python ./python/ap3_script.py";
            return_code = system(command.c_str());
            if (return_code != 0) {
                throw "An error occured during detectability results processing : " + return_code;
            }
            break;
        }
    }
    this->Define_Probabilities_2("output_file", min_detect);
}