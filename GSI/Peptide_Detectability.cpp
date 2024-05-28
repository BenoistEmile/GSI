#include "Model.h"
#include <filesystem>

//__________________________________________________________________________________________________________

void Model::AP3_Fasta() {
    std::ofstream temp_fasta(std::filesystem::current_path() / "local" / "temp_fasta.fasta");
    for (auto& protein : proteins) {
        temp_fasta << ">" << protein->Get_Accession() << std::endl;
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

void Model::Peptide_Detectability(int detectability_model, int minimum_number_of_amino_acids, int maximum_number_of_amino_acids) {
    In_Silico_Digestion_2("digestion_file", minimum_number_of_amino_acids, maximum_number_of_amino_acids);
    switch (detectability_model) {
        case 1: {
            std::string command = "conda run -n detectability python ./python/DbyDeep_script.py";
            int return_code = system(command.c_str());
            if (return_code != 0) {
                std::cout << "An error occured during peptide detectability prediction" << std::endl;
            }
            
        }
        break;
        case 2: {
            AP3_Fasta();
        }
        break;
    }
}