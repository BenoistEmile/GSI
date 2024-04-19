#include "Model.h"

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