#include "Model.h"

#include <stdlib.h>
#include <time.h>

//__________________________________________________________________________________________________________

void Model::Define_Probabilities(float value) {
    for (unsigned int i = 0; i < peptides.size(); ++i) {
        (*peptides[i]).Define_Probabilities(value);
    }
}

void Model::Define_Probabilities(bool per_protein) {
    srand((unsigned int)time(NULL));
    if (per_protein) {
        float value;
        for (unsigned int i = 0; i < proteins.size(); ++i) {
            value = rand() / (float)RAND_MAX;
            for (unsigned int j = 0; j < (*proteins[i]).Get_Peptides().size(); ++j) {
                (*(peptides[(*proteins[i]).Get_Peptides()[j]])).Define_Probabilities(value, i);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < peptides.size(); ++i) {
            auto proteins = (*peptides[i]).Get_Proteins();
            for (auto iter = proteins.begin(); iter != proteins.end(); ++iter) {
                for (unsigned int j = 0; j < iter->second.size(); ++j) {
                    (*peptides[i]).Define_Probabilities(rand() / (float)RAND_MAX, iter->first, j);
                }
            }
        }
    }
}

void Model::Define_Probabilities(void (define_probabilities_function)(std::vector<Protein*> proteins, std::vector<Peptide*> peptides)) {
    define_probabilities_function(proteins, peptides);
}