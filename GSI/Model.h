
#pragma once

#include "Score.h"
#include "Peptide.h"
#include "Protein.h"
#include "Fonction.h"

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <random>

/*
* Ce fichier contient la d�finition du type Model ainsi que du type Solution.
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
    * Modifie la quantit� associ� � une prot�ine dans la solution
    */
    void Add_Protein(std::size_t protein ,float abundance) {
        abundances[protein] = abundance;
    };

    /*
    * Permet d'ajouter une ar�te s�lectionn�e dans la solution
    */
    void Add_Score(const Score* score) {
        identifications.push_back(score->Get_Edge());
    };

    /*
    * Vide l'enti�ret� de la solution
    */
    void Clear() {
        abundances.clear();
        for (const Identification* identification : identifications) {
            delete identification;
        }
        identifications.clear();
    };

    /*
    * Indique si une protéine est présente dans la solution.
    */
    const bool Is_In_Solution(std::size_t protein) const {
        return abundances.contains(protein);
    }

    /*
    * Accès aux abondances des protéines dans la solution.
    */
    const float Get_Abundance(std::size_t protein) const {
        return abundances.at(protein);
    }

    const std::unordered_map<std::size_t, float>& Get_abundances() const {
        return abundances;
    }

    /*
    * Permet d'afficher la solution.
    * On retrouve en premier la liste des prot�ines s�lectionn�es avec leurs abondances.
    * Puis la liste des ar�tes s�lectionn�es (peptide, spectre)
    * Puis le nombre d'ar�tes s�lectionn�es � tort
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

    /*
    Enregistre la solution dans un fichier texte.
    */
    void Save(const std::vector<Peptide*> peptides , std::vector<Spectrum*> spectra, const std::string file_name, bool overwrite = false) const {
        std::ofstream output_file;
        std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
        file_path += ".txt";
        if (!fileExists(file_path) || overwrite) {
            output_file.open(file_path);
            output_file << "Selected proteins :" << std::endl;
            for (auto& couple : abundances) {
                output_file << "   - " << couple.first << " : " << couple.second << std::endl;
            }
            for (const Identification* identification : identifications) {
                output_file << *identification << std::endl;
            }
            if (spectra.size() && spectra[0]->Is_Simulated()) {
                output_file << "\nNumber of wrong selected edges : ";
                unsigned int compteur = 0;
                for (const Identification* identification : identifications) {
                    if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                        compteur++;
                    }
                }
                output_file << compteur << std::endl;
            }
            output_file.close();
            std::cout << "Saved solution to " << file_name << std::endl;
        }
        else {
            std::cout << "ERROR : There already is a file named : " << file_name << std::endl;
        }
    }

    /*
    * Enregistre la solution dans le fichier de log fourni.
    */
    void Save(std::vector<Spectrum*> spectra, std::ofstream& output_file, bool proteins = true, bool ident = true) const {
        output_file << "SOLUTION :" << std::endl;
        output_file << "Selected " << abundances.size() << " proteins." << std::endl;
        if (proteins) {
            output_file << "Selected proteins :" << std::endl;
            for (auto& couple : abundances) {
                output_file << "   - " << couple.first << " : " << couple.second << std::endl;
            }
        }
        if (ident) {
            for (const Identification* identification : identifications) {
                output_file << *identification << std::endl;
            }
        }
        if (spectra.size() && spectra[0]->Is_Simulated()) {
            output_file << "\nNumber of wrong selected edges : ";
            unsigned int compteur = 0;
            for (const Identification* identification : identifications) {
                if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                    compteur++;
                }
            }
            output_file << compteur << std::endl;
        }
        output_file << std::endl << std::endl;
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
    * Renvoie un fichier dans lequel les paramètres et résultats de l'exécution du programme pourront être enregistrés.
    */
    std::ofstream Open_Output_File(std::string file_name);

    /*
    * Charge le fichier de prot�ines en param�tre. Il doit �tre au format fasta
    */
    void Load_Proteins(const std::string file_name);
    /*
    * Charge le fichier de prot�ines en param�tre avec le parser donn�. Le parser est une fonction prenant un std::ifstream en param�tre et renvoyant un std::vector<Protein*>.
    */
    void Load_Proteins(const std::string file_name, std::vector<Protein*>(parser)(std::ifstream& file));
    /*
    * Charge le fichier de protéines en paramètre. Il doit être au format fasta.
    * Sauvegarde des informations sur le jeu de données chargé (nom du fichier, nombre de protéines) dans le fichier fourni.
    */
    void Load_Proteins(const std::string file_name, std::ofstream& output_file);

    /*
    * Effectue la digestion in-silico des prot�ines avec de la trypsine. Seuls les peptides dont la taille respecte les limites en param�tre sont g�n�r�s.
    */
    void In_Silico_Digestion(int minimum_number_of_amino_acids = 7 ,int maximum_number_of_amino_acids = 25);
    /*
    * Effectue la digestion in-silico des prot�ines avec la fonction de digestion en param�tre. Seuls les peptides dont la taille respecte les limites en param�tre sont g�n�r�s
    * La fonction de digestion prend un std::string& (s�quence de prot�ine) en param�tre et renvoie un std::vector<std::string> (vecteur de s�quences de peptide)
    */
    void In_Silico_Digestion(std::vector<std::string>(digest_function)(const std::string& sequence) ,int minimum_number_of_amino_acids = 7,int maximum_number_of_amino_acids = 7);
    /*
    * Effectue la digestion in-silico des protéines avec de la trypsine. Seuls les peptides dont la taille respecte les limites en paramètre sont générés.
    * Génère un fichier contenant les peptides générés (au format utilisé par DbyDeep).
    */
    void In_Silico_Digestion(std:: string file_name, int minimum_number_of_amino_acids = 7 ,int maximum_number_of_amino_acids = 25);
    /*
    * Effectue la digestion in-silico des protéines avec de la trypsine. Seuls les peptides dont la taille respecte les limites en paramètre sont générés.
    * Génère un fichier contenant les peptides générés (au format utilisé par DbyDeep).
    * Sauvegarde le nombre de peptides générés dans le fichier fourni.
    */
    void In_Silico_Digestion(std:: string file_name, std::ofstream& log_file, int minimum_number_of_amino_acids = 7 ,int maximum_number_of_amino_acids = 25);

    /*
    * 
    */
    void Peptide_detectability(std::string env_name = "Dby_Deep", std::string digestion_file_name = "digestion") const;

    /*
    * G�n�re les spectres th�oriques � partir des spectres th�oriques avec la possibilit� d'ajouter des modifications.
    * Le param�tre modifications associe des caract�res (acides amin�s) � des masses (modification).
    */
    void Build_Theoretical_Spectra(std::unordered_map<char, double> modifications = {});

    /*
    * D�finie les probabilit�s sur les ar�tes prot�ines-peptides en affectant la valeur en param�tre � chaque ar�te.
    */
    void Define_Probabilities(float value);
    /*
    * D�finie les probabilit�s al�atoires (loi uniforme) sur les ar�tes prot�ine-peptide. Si le param�tre per_protein est � vrai, alors la probabilit� des ar�tes adjacentes � une m�me prot�ine obtiennent la m�me probabilit�.
    */
    void Define_Probabilities(bool per_protein = false);
    /*
    * D�finie les probabilit�s sur les ar�tes prot�ines-peptides selon la fonction en param�tre. Cette fonction prend en param�tre un std::vector<Protein*> (liste des prot�ines) et un std::vector<Peptide*> (liste des peptides).
    * Voir la classe Peptide pour avoir des d�tails sur la fa�on dont on peut attribuer les probabilit�s.
    */
    void Define_Probabilities(void (define_probabilities_function)(std::vector<Protein*> proteins, std::vector<Peptide*> peptides));
    /*
    * Définit les probabilités sur les arêtes protéines-peptides à partir d'un fichier csv en paramètre.
    * Le fichier csv doit contenir les colonnes suivantes : protein_id, peptide_id, Prob.
    * La digestion in-silico doit avoir été réalisée en demandant la génération d'un fichier.
    */
    void Define_Probabilities(const std::string file_name);
    /*
    * Définit les probabilités sur les arêtes protéines-peptides à partir d'un fichier csv en paramètre.
    * Le fichier csv doit contenir les colonnes suivantes : protein_id, peptide_id, Prob.
    * La digestion in-silico doit avoir été réalisée en demandant la génération d'un fichier.
    * Sauvegarde des informations (méthode de génération des probabilités) dans le fichier fourni.
    */
    void Define_Probabilities(const std::string file_name, std::ofstream& output_file);

    /*
    * Charge le fichier de spectres en param�tre au format ms2
    */
    void Load_Spectra(const std::string file_name);
    /*
    * Charge le fichier de spectres en param�tre avec le parser donn�. Le parser prend en param�tre un std::ifstream& et renvoie un std::vector<Spectrum*>.
    */
    void Load_Spectra(const std::string file_name, std::vector<Spectrum*>(parser)(std::ifstream& file));
    /*
    * Charge le fichier de spectres en paramètre au format ms2
    * Sauvegarde des informations sur les spectres (fichier d'origine, nombre de spectres) dans le fichier fourni.
    */
    void Load_Spectra(const std::string file_name, std::ofstream& output_file);

    /*
    * G�n�re un ensemble de spectres simul�s � partir des prot�ines, des peptides th�oriques et des probabilit�s.
    * Pour cela, il faut fournir un �chantillon sous la forme d'une map associant des num�ros de prot�ines � des abondances.
    * Le param�tre error_rate permet de s'�loigner d'un "monde parfait" lorsqu'on le rapproche de 1.
    */
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, float error_rate = 0.f);
    /*
    * G�n�re un ensemble de spectres simul�s � partir des prot�ines, des peptides th�oriques et des probabilit�s.
    * Pour cela, il faut fournir un �chantillon sous la forme d'une map associant des num�ros de prot�ines � des abondances.
    * La fa�on dont les spectres sont s�lectionn�s est directement d�pendante de la fonction � fournir en param�tre.
    * Cette fonction prend en param�tre un std::unordered_map<std::size_t, unsigned int> (correspond au sample fourni), un std::vector<Protein*>, un std::vector<Peptide*> et un std::vector<Spectrum*> (liste de spectre � remplir).
    */
    void Simulated_Sample(std::unordered_map<std::size_t, unsigned int> sample, void (simulated_sample_function)(std::unordered_map<std::size_t, unsigned int> sample, std::vector<Protein*> proteins, std::vector<Peptide*> peptides, std::vector<Spectrum*> spectra));

    //void Load_Scores(const std::string file_name);
    /*
    * Permet d'affecter des scores aux ar�tes spectres-peptides � partir d'un fichier � indiquer en param�tre. Le parser doit �tre fourni.
    * Il prend en param�tre un std::ifstream& et renvoie un std::vector<Score*> (liste des scores).
    */
    void Load_Scores(const std::string file_name , std::vector<Score*>(parser)(std::ifstream& file));
    /*
    * Permet d'affecter des scores aux ar�tes spectres-peptides � partir de la fonction en param�tre.
    * Celle-ci doit prendre en param�tre un std::vector<Peptide*>, un std::vector<Spectrum*> et un std::vector<Score*> (liste de scores � remplir).
    */
    void Compute_Score(void (compute_score_function)(std::vector<Peptide*> &peptides, std::vector<Spectrum*> &spectra, std::vector<Score*> &scores));
    /*
    * Permet d'affecter des scores aux ar�tes spectres-peptides dans le cas d'un �chantillon simul� en ajoutant une ar�te pour chaque PSM correct
    * De plus, pour chaque spectre, plusieurs ar�tes sont ajout�es avec une position al�atoire. Ce nombre d'ar�tes correspond au param�tre "randoms".
    * Que l'ar�te soit l'ar�te correcte ou al�atoire, leur score est toujours de 0 (le meilleur).
    */
    void Compute_Score(unsigned int randoms = 0);
    /*
    * Score calcul�s � partir de l'algorithme SpecOMS. Les param�tres correspondent aux param�tres de SpecOMS. Si le param�tre maximum_number_of_edges est fix� � 0, il est consid�r� � +infiny.
    */
    void Compute_Score_SpecOMS(unsigned int minimum_number_of_masses = 0, unsigned int maximum_number_of_masses = 99999, int accuracy = 2, unsigned int number_of_copies = 2, unsigned int threshold = 0, unsigned int maximum_number_of_edges = 0);
    /*
    * Score calcul�s � partir de l'algorithme SpecOMS. Les param�tres correspondent aux param�tres de SpecOMS. Si le param�tre maximum_number_of_edges est fix� � 0, il est consid�r� � +infiny.
    * Sauvegarde des informations (paramètres, nombre de scores calculés, méthode utilisée, temps d'exécution) dans le fichier fourni.
    */
    void Compute_Score_SpecOMS(std::ofstream& output_file, unsigned int minimum_number_of_masses = 0, unsigned int maximum_number_of_masses = 99999, int accuracy = 2, unsigned int number_of_copies = 2, unsigned int threshold = 0, unsigned int maximum_number_of_edges = 0);

    /*
    * Calcule une solution pour le mod�le courant. psi1 correspond au coefficient de l'objectif sur les Deltas, psi2 correspond au coefficient de l'objectif sur les ar�tes spectre-peptide.
    */
    int Solve(const float psi1 = 0.5f ,const float psi2 = 0.5f);
    /*
    * Calcule une solution pour le mod�le courant. psi1 correspond au coefficient de l'objectif sur les Deltas, psi2 correspond au coefficient de l'objectif sur les ar�tes spectre-peptide.
    * Sauvegarde des informations (paramètres, durée d'exécution) dans le fichier fourni.
    */
    int Solve(std::ofstream& output_file, float psi1 = 0.5f ,const float psi2 = 0.5f);

    /*
    * Affiche la solution courante.
    */
    void Print_Solution() const {
        solution.Print(peptides ,spectra);
    };

    /*
    * Enregistre la solution
    */
    void Save_Solution(std::string file_name, bool overwrite = false) const {
        solution.Save(spectra, file_name, overwrite);
    }
    /*
    * Enregistre la solution dans le fichier de log fourni.
    */
    void Save_Solution(std::ofstream& output_file, bool proteins = true, bool ident = true) const {
        solution.Save(spectra, output_file, proteins, ident);
    }

    /*
    * Renvoie la solution.
    */
    const Solution& Get_Solution() const {
        return solution;
    }

    /*
    * Affiche et enregistre (solution/file_name.csv) les résultats obtenus ainsi que l'échantillon synthétique fourni au modèle.
    */
    void Analyse_Solution(std::unordered_map<std::size_t, unsigned int> sample, std::string file_name) const;

    //void Evaluate_Solution();

    //void Save() const;
    //void Load();

    /*
    * Permet de r�initialiser certaines parties de l'instance courante.
    * A noter que la r�initialisation de certaines parties peut entra�ner la r�initialisation d'autres sans qu'elles ne soient demand�es explicitement.
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

    void Add_peptide(std::string sequence) {
        peptides.push_back(new Peptide(peptides.size(), sequence));
    }

    /*
    * Génère un échantillon aléatoire de protéines pour les fonctions SImulated_Sample.
    * 
    */
    std::unordered_map<std::size_t, unsigned int> Random_Sample(unsigned int min_proteins = 1, unsigned int max_proteins = 0, unsigned int min_abundance = 1, unsigned int max_abundance = 10);

    /*
    * Teste le modèle sur des données synthétiques. Le modèle doit avoir été préparé avec des protéines, des peptides, des arêtes et scores protéine-peptide.
    */
   void Run_Tests_Sample_Data(unsigned int n_tests, unsigned int min_proteins = 1, unsigned int max_proteins = 0, unsigned int min_abundance = 1, unsigned int max_abundance = 10, std::string file_name = "test_results");

   /*
   * Teste le modèle sur des données réelles. => à modifier car on obtient le même résultat à chaque run
   * Le modèle doit avoir été préparé avec les éléments suivants :
   * - Protéines (Load_Proteins)
   * - Peptides (In_Silico_Digestion)
   * - Arêtes Peptides - Protéines (Define_Probabilities)
   * - Scores Peptides - Spectres (Load_Spectra + Compute_Score_SpecOMS)
   */
   void Run_Tests_Real_Data(unsigned int n_tests, std::string file_name);
};