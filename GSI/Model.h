
#pragma once

#include "Score.h"
#include "Peptide.h"
#include "Protein.h"
#include "Fonction.h"

#include <vector>
#include <string>

/*
* Ce fichier contient la définition du type Model ainsi que du type Solution.
*/

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

    /*
    * Modifie la quantité associé à une protéine dans la solution
    */
    void Add_Protein(std::size_t protein ,float abundance) {
        abundances[protein] = abundance;
    };

    /*
    * Permet d'ajouter une arête sélectionnée dans la solution
    */
    void Add_Score(const Score* score) {
        identifications.push_back(score->Get_Edge());
    };

    /*
    * Vide l'entièreté de la solution
    */
    void Clear() {
        abundances.clear();
        for (const Identification* identification : identifications) {
            delete identification;
        }
        identifications.clear();
    };

    /*
    * Permet d'afficher la solution.
    * On retrouve en premier la liste des protéines sélectionnées avec leurs abondances.
    * Puis la liste des arêtes sélectionnées (peptide, spectre)
    * Puis le nombre d'arêtes sélectionnées à tort
    */
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

    /*
    * Charge le fichier de protéines en paramètre. Il doit être au format fasta
    */
    void Load_Proteins(const std::string file_name);
    /*
    * Charge le fichier de protéines en paramètre avec le parser donné. Le parser est une fonction prenant un std::ifstream en paramètre et renvoyant un std::vector<Protein*>.
    */
    void Load_Proteins(const std::string file_name, std::vector<Protein*>(parser)(std::ifstream& file));

    /*
    * Effectue la digestion in-silico des protéines avec de la trypsine. Seuls les peptides dont la taille respecte les limites en paramètre sont générés.
    */
    void In_Silico_Digestion(int minimum_number_of_amino_acids = 7 ,int maximum_number_of_amino_acids = 25);
    /*
    * Effectue la digestion in-silico des protéines avec la fonction de digestion en paramètre. Seuls les peptides dont la taille respecte les limites en paramètre sont générés
    * La fonction de digestion prend un std::string& (séquence de protéine) en paramètre et renvoie un std::vector<std::string> (vecteur de séquences de peptide)
    */
    void In_Silico_Digestion(std::vector<std::string>(digest_function)(const std::string& sequence) ,int minimum_number_of_amino_acids = 7,int maximum_number_of_amino_acids = 7);

    /*
    * Génère les spectres théoriques à partir des spectres théoriques avec la possibilité d'ajouter des modifications.
    * Le paramètre modifications associe des caractères (acides aminés) à des masses (modification).
    */
    void Build_Theoretical_Spectra(std::unordered_map<char, double> modifications = {});

    /*
    * Définie les probabilités sur les arêtes protéines-peptides en affectant la valeur en paramètre à chaque arête.
    */
    void Define_Probabilities(float value);
    /*
    * Définie les probabilités aléatoires (loi uniforme) sur les arêtes protéine-peptide. Si le paramètre per_protein est à vrai, alors la probabilité des arêtes adjacentes à une même protéine obtiennent la même probabilité.
    */
    void Define_Probabilities(bool per_protein = false);
    /*
    * Définie les probabilités sur les arêtes protéines-peptides selon la fonction en paramètre. Cette fonction prend en paramètre un std::vector<Protein*> (liste des protéines) et un std::vector<Peptide*> (liste des peptides).
    * Voir la classe Peptide pour avoir des détails sur la façon dont on peut attribuer les probabilités.
    */
    void Define_Probabilities(void (define_probabilities_function)(std::vector<Protein*> proteins, std::vector<Peptide*> peptides));

    /*
    * Charge le fichier de spectres en paramètre au format ms2
    */
    void Load_Spectra(const std::string file_name);
    /*
    * Charge le fichier de spectres en paramètre avec le parser donné. Le parser prend en paramètre un std::ifstream& et renvoie un std::vector<Spectrum*>.
    */
    void Load_Spectra(const std::string file_name, std::vector<Spectrum*>(parser)(std::ifstream& file));

    /*
    * Génère un ensemble de spectres simulés à partir des protéines, des peptides théoriques et des probabilités.
    * Pour cela, il faut fournir un échantillon sous la forme d'une map associant des numéros de protéines à des abondances.
    * Le paramètre error_rate permet de s'éloigner d'un "monde parfait" lorsqu'on le rapproche de 1.
    */
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, float error_rate = 0.f);
    /*
    * Génère un ensemble de spectres simulés à partir des protéines, des peptides théoriques et des probabilités.
    * Pour cela, il faut fournir un échantillon sous la forme d'une map associant des numéros de protéines à des abondances.
    * La façon dont les spectres sont sélectionnés est directement dépendante de la fonction à fournir en paramètre.
    * Cette fonction prend en paramètre un std::unordered_map<std::size_t, unsigned int> (correspond au sample fourni), un std::vector<Protein*>, un std::vector<Peptide*> et un std::vector<Spectrum*> (liste de spectre à remplir).
    */
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, void (simulated_sample_function)(std::unordered_map<std::size_t, unsigned int> sample, std::vector<Protein*> proteins, std::vector<Peptide*> peptides, std::vector<Spectrum*> spectra));

    //void Load_Scores(const std::string file_name);
    /*
    * Permet d'affecter des scores aux arêtes spectres-peptides à partir d'un fichier à indiquer en paramètre. Le parser doit être fourni.
    * Il prend en paramètre un std::ifstream& et renvoie un std::vector<Score*> (liste des scores).
    */
    void Load_Scores(const std::string file_name , std::vector<Score*>(parser)(std::ifstream& file));
    /*
    * Permet d'affecter des scores aux arêtes spectres-peptides à partir de la fonction en paramètre.
    * Celle-ci doit prendre en paramètre un std::vector<Peptide*>, un std::vector<Spectrum*> et un std::vector<Score*> (liste de scores à remplir).
    */
    void Compute_Score(void (compute_score_function)(std::vector<Peptide*> &peptides, std::vector<Spectrum*> &spectra, std::vector<Score*> &scores));
    /*
    * Permet d'affecter des scores aux arêtes spectres-peptides dans le cas d'un échantillon simulé en ajoutant une arête pour chaque PSM correct
    * De plus, pour chaque spectre, plusieurs arêtes sont ajoutées avec une position aléatoire. Ce nombre d'arêtes correspond au paramètre "randoms".
    * Que l'arête soit l'arête correcte ou aléatoire, leur score est toujours de 0 (le meilleur).
    */
    void Compute_Score(unsigned int randoms = 0);
    /*
    * Score calculés à partir de l'algorithme SpecOMS. Les paramètres correspondent aux paramètres de SpecOMS. Si le paramètre maximum_number_of_edges est fixé à 0, il est considéré à +infiny.
    */
    void Compute_Score_SpecOMS(unsigned int minimum_number_of_masses = 0, unsigned int maximum_number_of_masses = 99999, int accuracy = 2, unsigned int number_of_copies = 2, unsigned int threshold = 0, unsigned int maximum_number_of_edges = 0);

    /*
    * Calcule une solution pour le modèle courant. psi1 correspond au coefficient de l'objectif sur les Deltas, psi2 correspond au coefficient de l'objectif sur les arêtes spectre-peptide.
    */
    void Solve(const float psi1 = 0.5f ,const float psi2 = 0.5f);

    /*
    * Affiche la solution courante.
    */
    void Print_Solution() const {
        solution.Print(peptides ,spectra);
    };

    //void Evaluate_Solution();

    //void Save() const;
    //void Load();

    /*
    * Permet de réinitialiser certaines parties de l'instance courante.
    * A noter que la réinitialisation de certaines parties peut entraîner la réinitialisation d'autres sans qu'elles ne soient demandées explicitement.
    */
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
