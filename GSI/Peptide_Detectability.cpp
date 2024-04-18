#include "Model.h"

//__________________________________________________________________________________________________________

void Model::Peptide_detectability(std::string env_name, std::string digestion_file_name) const {
    std::string command = "conda run -n " + env_name + 
    " --cwd ~/stage/code/DbyDeep-main \
    python dbydeep_model.py \
    --data-path ../GSI/GSI/data/digestion/" + digestion_file_name + ".csv \
    --save-path ../GSI/GSI/data/digestion/ \
    --job-name " + digestion_file_name + "_result";
    // std::string command = ("conda run -n " + env_name + " --cwd ~/stage/code/DbyDeep-main pwd");
    // std::cout << "Command : " << command << std::endl;
    int return_code = system(command.c_str());
    if (return_code != 0) {
        std::cout << "An error occured during peptide detectability prediction" << std::endl;
    }
}