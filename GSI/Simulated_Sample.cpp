#include "Model.h"

//__________________________________________________________________________________________________________

void Model::Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, float error_rate) {
    srand((unsigned int)time(NULL));
    std::unordered_map<std::size_t, std::size_t> positions;
    float probability;
    float final_probability;
    for (auto iter = sample.begin(); iter != sample.end(); ++iter) {
        positions = {};
        for (auto iter2 = (*proteins[iter->first]).Get_Peptides().begin(); iter2 != (*proteins[iter->first]).Get_Peptides().end(); ++iter2) {
            if (positions.contains(*iter2)){
                positions.at(*iter2) += 1;
            }
            else {
                positions[*iter2] = 0;
            }
            const std::vector<Pic*>* pics = (*peptides[*iter2]).Get_Pics(); //.Build_Spectrum();
            probability = (*peptides[*iter2]).Get_Proteins().at(iter->first)[positions[*iter2]];
            final_probability = error_rate * ((1 - round(probability)) - probability) + probability;
            const Origin* origin = new Origin(*iter2, iter->first, positions[*iter2]);
            for (unsigned int i = 0; i < iter->second; ++i) {
                if (rand() / (float)RAND_MAX <= final_probability) {
                    spectra.push_back(new Spectrum(spectra.size() ,pics ,origin ,true));
                }
            }
        }
    }
}

void Model::Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, void (simulated_sample_function)(std::unordered_map<std::size_t, unsigned int> sample, std::vector<Protein*> proteins, std::vector<Peptide*> peptides, std::vector<Spectrum*> spectra)) {
    simulated_sample_function(sample, proteins, peptides, spectra);
}