
#pragma once

#include "Score.h"
#include "Peptide.h"
#include "Protein.h"
#include "Fonction.h"

#include <vector>
#include <string>

struct Solution {
    std::unordered_map<std::size_t, float> abundances;
    std::vector<const Identification*> identifications; 

    ~Solution() {
        for (const Identification* identification : identifications) {
            delete identification;
        }
    }
    /*
    void operator+=(const Score& score) {
        identifications.push_back(score.Get_Edge());
    };
    */

    void Add_Protein(std::size_t protein ,float abundance) {
        abundances[protein] = abundance;
    };

    void Add_Score(const Score* score) {
        identifications.push_back(score->Get_Edge());
    };

    void Clear() {
        abundances.clear();
        for (const Identification* identification : identifications) {
            delete identification;
        }
        identifications.clear();
    };

    void Print(const std::vector<Peptide*> peptides , std::vector<Spectrum*> spectra) const {
        std::cout << "\nSelected proteins :" << std::endl;
        for (auto& couple : abundances) {
            std::cout << "   - " << couple.first << " : " << couple.second << std::endl;
        }
        for (const Identification* identification : identifications) {
            std::cout << *identification << std::endl;
        }
        if (spectra.size() && spectra[0]->Is_Simulated()) {
            std::cout << "\nNumber of wrong selected edges : ";
            unsigned int compteur = 0;
            for (const Identification* identification : identifications) {
                if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                    compteur++;
                }
            }
            std::cout << compteur << std::endl;
        }
        /*
        std::cout << "\nSelected edges :" << std::endl;
        for (const Identification* identification : identifications) {
            std::cout << "   - (p" << identification->peptide << " ,s" << identification->spectrum << ")" << std::endl;
        }
        */
    }
};

//__________________________________________________________________________________________________________

class Model {
private:
    std::vector<Protein*> proteins;
    std::vector<Peptide*> peptides;
    std::unordered_map<std::string, std::size_t> peptides_sequences;
    std::vector<Spectrum*> spectra;
    std::vector<Score*> scores;
    Solution solution;
public:
    Model();
    ~Model();

    void Load_Proteins(const std::string file_name);
    void Load_Proteins(const std::string file_name, std::vector<Protein*>(parser)(std::ifstream& file));

    void In_Silico_Digestion(int minimum_number_of_amino_acids = 7 ,int maximum_number_of_amino_acids = 25);
    void In_Silico_Digestion(std::vector<std::string>(digest_function)(const std::string& sequence) ,int minimum_number_of_amino_acids = 7,int maximum_number_of_amino_acids = 7);

    void Build_Theoretical_Spectra(std::unordered_map<char, double> modifications = {});

    void Define_Probabilities(float value);
    void Define_Probabilities(bool per_protein = false);
    void Define_Probabilities(void (define_probabilities_function)(std::vector<Protein*> proteins, std::vector<Peptide*> peptides));

    void Load_Spectra(const std::string file_name);
    void Load_Spectra(const std::string file_name, std::vector<Spectrum*>(parser)(std::ifstream& file));
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, float error_rate = 0.f);
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, void (simulated_sample_function)(std::unordered_map<std::size_t, unsigned int> sample, std::vector<Protein*> proteins, std::vector<Peptide*> peptides, std::vector<Spectrum*> spectra));

    //void Load_Scores(const std::string file_name);
    void Load_Scores(const std::string file_name , std::vector<Score*>(parser)(std::ifstream& file));
    void Compute_Score(void (compute_score_function)(std::vector<Peptide*> &peptides, std::vector<Spectrum*> &spectra, std::vector<Score*> &scores));
    void Compute_Score(unsigned int randoms = 0);
    void Compute_Score2(unsigned int minimum_number_of_masses = 0, unsigned int maximum_number_of_masses = 99999, int accuracy = 2, unsigned int number_of_copies = 2, unsigned int threshold = 0, unsigned int maximum_number_of_edges = 0);

    void Solve(const float psi1 = 0.5f ,const float psi2 = 0.5f);

    void Print_Solution() const {
        solution.Print(peptides ,spectra);
    };

    //void Evaluate_Solution();

    //void Save() const;
    //void Load();

    void Clear(bool clear_proteins = true, bool clear_peptides = true, bool clear_spectra = true, bool clear_scores = true, bool clear_solution = true);

    const std::size_t Number_Of_Proteins() const;
    const std::size_t Number_Of_Peptides() const;
    const std::size_t Number_Of_Spectra() const;
    const std::size_t Number_Of_Scores() const;
    const Protein& Get_Protein(std::size_t protein) const;
    const Peptide& Get_Peptide(std::size_t peptide) const;
    const Spectrum& Get_Spectrum(std::size_t spectrum) const;
    const Score& Get_Score(std::size_t score) const;
};
