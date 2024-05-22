#include "Model.h"
#define FMT_HEADER_ONLY
#include "fmt/format.h"

//__________________________________________________________________________________________________________

void Model::Analyse_Solution(std::unordered_map<std::size_t, unsigned int> sample, std::string file_name) const {
    std::ofstream output_file;
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".csv";
    if (!fileExists(file_path)) {
        std::cout << "Proteins in sample :" << std::endl;
        for (auto iter = sample.begin(); iter != sample.end(); iter++) {
            if (solution.Is_In_Solution(iter->first)) {
                std::cout << "- " << iter->first << " : Quantity in sample : " << iter->second << ", Computed abundance : " << solution.Get_Abundance(iter->first) << std::endl;
            }
            else {
                std::cout << "- " << iter->first << " : Was not identified by the model." << std::endl;
            }
        }
        std::cout << "Proteins identified out of the sample :" << std::endl;
        for (auto iter = solution.abundances.begin(); iter != solution.abundances.end(); iter++) {
            if (!sample.contains(iter->first)) {
                std::cout << "- " << iter->first << " : Quantity in sample : 0, Computed abundance : " << iter->second << std::endl;
            }
        }
    }
    else {
        std::cout << "ERROR : There already is a file named : " << std::endl;
    }
}

//__________________________________________________________________________________________________________

std::unordered_map<std::size_t, unsigned int> Model::Random_Sample(unsigned int min_proteins, unsigned int max_proteins, unsigned int min_abundance, unsigned int max_abundance) {
    std::unordered_map<std::size_t, unsigned int> sample;
    std::vector<int> proteins_ids(proteins.size());
    std::vector<size_t> sampled_proteins;
    std::iota(proteins_ids.begin(), proteins_ids.end(), 0);
    if (max_proteins == 0) {
        max_proteins = proteins.size();
    }
    std::ranges::sample(proteins_ids, std::back_inserter(sampled_proteins), round(rand()/(float)RAND_MAX*(max_proteins - min_proteins)) + min_proteins, std::mt19937{std::random_device{}()});
    for (auto iter = sampled_proteins.begin(); iter != sampled_proteins.end(); iter++) {
        sample[*iter] = round(rand()/(float)RAND_MAX*(max_abundance - min_abundance)) + min_abundance;
        std::cout << *iter << " : " << sample.at(*iter) << std::endl;
    }
    return sample;
}

//__________________________________________________________________________________________________________

void Model::Run_Tests_Sample_Data(unsigned int n_tests, unsigned int min_proteins, unsigned int max_proteins, unsigned int min_abundance, unsigned int max_abundance, std::string file_name) {
    std::ofstream output_file;
    std::vector<int> times;
    std::vector<float> errors;
    int min_time, max_time;
    float mean_error, med_error, min_error, max_error, mean_time, med_time;
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".csv";
    for (int i = 0; i < n_tests; i++) {
        std::unordered_map<std::size_t, unsigned int> sample = this->Random_Sample(min_proteins, max_proteins, min_abundance, max_abundance);
        this->Simulated_Sample(sample);
        this->Compute_Score(1);
        times.push_back(this->Solve(0.5f, 0.5f));
        for (auto iter = sample.begin(); iter != sample.end(); iter++) {
            if (this->solution.Is_In_Solution(iter->first)) {
                errors.push_back(abs(iter->second - this->solution.Get_Abundance(iter->first)));
            }
            else {
                errors.push_back(iter->second);
            }
        }
        for (auto iter = this->solution.abundances.begin();iter != this->solution.abundances.end(); iter++) {
            if (!sample.contains(iter->first)) {
                errors.push_back(iter->second);
            }
        }
        this->Clear(false, false, true, true, true);
    }
    mean_error = std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();
    std::size_t n = errors.size()/2;
    std::nth_element(errors.begin(), errors.begin(), errors.end());
    med_error = errors[n];
    min_error = *std::min_element(errors.begin(), errors.end());
    max_error = *std::max_element(errors.begin(), errors.end());
    mean_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    n = times.size()/2;
    std::nth_element(times.begin(), times.begin(), times.end());
    med_time = times[n];
    min_time = *std::min_element(times.begin(), times.end());
    max_time = *std::max_element(times.begin(), times.end());
    if (!fileExists(file_path)) {
        output_file.open(file_path);
        output_file << "loaded proteins,peptides digested,n_tests,min_proteins,max_proteins,min_abundance,max_abundance,mean_error,med_error,min_error,max_error,mean_time,med_time,min_time,max_time" << std::endl;
    }
    else {
        output_file.open(file_path, std::ios_base::app);
    }
    output_file << this->Number_Of_Proteins() << "," << this->Number_Of_Peptides() << "," << n_tests << "," << min_proteins << "," << max_proteins << "," << min_abundance << "," << max_abundance // Parameters
    << "," << mean_error << "," << med_error << "," << min_error << "," << max_error << "," << mean_time << "," << med_time << "," << min_time << "," << max_time << std::endl; // Results
    output_file.close();
}

//__________________________________________________________________________________________________________

void Model::Run_Tests_Real_Data(unsigned int n_tests, std::string file_name) {
    std::ofstream output_file;
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".csv";
    Solution sol;
    std::unordered_map<std::size_t, float> abundances;
    std::unordered_map<std::size_t, std::vector<float>> identifications;
    std::vector<int> times;
    int min_time, max_time;
    float mean_time, std_time = 0;
    float mean_abundance, min_abundance, max_abundance, std_abundance;
    for (int i = 0; i < n_tests; i++) {
        times.push_back(this->Solve());
        sol = this->Get_Solution();
        abundances = sol.Get_abundances();
        for (auto iter_id = identifications.begin(); iter_id != identifications.end(); iter_id++) {
            if (!abundances.contains(iter_id->first)) {
                iter_id->second.push_back(0);
            }
        }
        for (auto iter_ab = abundances.begin(); iter_ab != abundances.end(); iter_ab++) {
            if (identifications.contains(iter_ab->first)) {
                identifications.at(iter_ab->first).push_back(iter_ab->second);
            }
            else {
                identifications[iter_ab->first] = {iter_ab->second};
            }
        }
        this->Clear(false, false, false, false, true);
    }
    if (!fileExists(file_path)) {
        output_file.open(file_path);
        output_file << "Protein_ID,";
        for (int i = 0; i < n_tests; i++) {
            output_file << "run_" << i << ",";
        }
        output_file << "mean,std_deviation,min,max" << std::endl;
        output_file << "time,";
        mean_time = std::accumulate(times.begin(), times.end(), 0.0) / n_tests;
        min_time = *std::min_element(times.begin(), times.end());
        max_time =*std::max_element(times.begin(), times.end());
        for (auto iter_time = times.begin(); iter_time != times.end(); iter_time++) {
            output_file << *iter_time << ",";
            std_time += pow(*iter_time - mean_time, 2);
        }
        std_time /= n_tests;
        output_file << mean_time << "," << std_time << "," << min_time << "," << max_time << std::endl;
        for (auto iter = identifications.begin(); iter != identifications.end(); iter++) {
            output_file << iter->first << ",";
            mean_abundance = std::accumulate(iter->second.begin(), iter->second.end(), 0.0) / n_tests;
            min_abundance = *std::min_element(iter->second.begin(), iter->second.end());
            max_abundance =*std::max_element(iter->second.begin(), iter->second.end());
            std_abundance = 0;
            for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
                output_file << *iter2 << ",";
                std_abundance += pow(*iter2 - mean_abundance, 2);
            }
            std_abundance /= n_tests;
            output_file << mean_abundance << "," << std_abundance << "," << min_abundance << "," << max_abundance << std::endl;
        }
        output_file.close();
    }
    else {
        std::cout << "The file " << file_path << "already exists." << std::endl;
    }
}

//__________________________________________________________________________________________________________

std::ofstream Model::Open_Output_File(std::string file_name) {
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".txt";
    if (fileExists(file_path)) {
        unsigned int i = 1;
        std::string new_file_name = file_name + "_" + std::to_string(i);
        while (fileExists(file_path.parent_path() / (new_file_name + ".txt"))) {
            i++;
            new_file_name = file_name + "_" + std::to_string(i);
        }
        file_name = new_file_name;
        file_path = std::filesystem::current_path() / "solution" / file_name;
        file_path += ".txt";
    }
    std::ofstream output_file;
    output_file.open(file_path);
    return output_file;
}

//__________________________________________________________________________________________________________

void Model::Test_Psi_Values(std::set<std::tuple<float, float>> psi_values, std::ofstream& output_file, std::string job_name) {
    std::string file_name;
    for (auto& iter : psi_values) {
        this->Solve(std::get<0>(iter), std::get<1>(iter));
        this->Save_Solution(output_file, false, false, false);
        file_name = job_name + "_" + std::to_string(std::get<0>(iter)) + "_" + std::to_string(std::get<1>(iter));
        this->Save_Solution(file_name, true, true, false, true);
        this->Clear(false, false, false, false, true);
    }
}

//__________________________________________________________________________________________________________

void Model::Run_Test(std::string prefix, float psi1, float psi2, unsigned int threshold, unsigned int max_edges, float min_detect, bool compute_detect) {
    std::string file_name;
    std::ofstream lower_edges_file;
    if (prefix != "") {
        file_name = prefix + "_" + std::to_string(threshold) + "_" + std::to_string(max_edges) + "_" + std::to_string((int)std::round(psi1)) + "_" + std::to_string((int)std::round(psi2)) + "_" + fmt::format("{0:.2f}", min_detect);
    }
    else {
        file_name = std::to_string(threshold) + "_" + std::to_string(max_edges) + "_" + std::to_string((int)std::round(psi1)) + "_" + std::to_string((int)std::round(psi2)) + "_" + fmt::format("{0:.2f}", min_detect);
    }
    std::filesystem::path file_path = std::filesystem::current_path() / "models" / ("lower_edges_" + file_name + ".csv");
    lower_edges_file.open(file_path);
    lower_edges_file << "Peptide,Spectrum,Score" << std::endl;
    this->In_Silico_Digestion(prefix);
    if (compute_detect) {
        this->Peptide_detectability("Dby_Deep", prefix);
    }
    this->Define_Probabilities(std::string(prefix + "_result.csv"), min_detect);
    this->Build_Theoretical_Spectra();
    this->Compute_Score_SpecOMS(0U, 99999U, 2, 0U, threshold, max_edges);
    std::cout << this->Number_Of_Scores() << std::endl;
    for (std::size_t iter_score = 0; iter_score < this->Number_Of_Scores(); iter_score++) {
        Score score = this->Get_Score(iter_score);
        lower_edges_file << score.peptide << "," << score.spectrum << "," << score.score << std::endl;
    }
    lower_edges_file.close();
    this->Solve(psi1, psi2);
    this->Save_Solution("results_" + file_name, true, true, false, true);
}

void Model::Run_Multiple_Tests(std::set<std::tuple<float, float, unsigned int, unsigned int, float>> parameters, std::string prefix) {
    for (auto& iter : parameters) {
        this->Run_Test(prefix, std::get<0>(iter), std::get<1>(iter), std::get<2>(iter), std::get<3>(iter), std::get<4>(iter), false);
        this->Clear(false, true, false, true, true);
    }
}

//__________________________________________________________________________________________________________

