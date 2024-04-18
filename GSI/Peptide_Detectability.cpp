#include "Model.h"
#include <filesystem>

//__________________________________________________________________________________________________________

void Model::Peptide_detectability(std::string env_name, std::string digestion_file_name) const {
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