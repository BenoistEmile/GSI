
#include "Model.h"
#include <unordered_map>
#include <set>

//__________________________________________________________________________________________________________

/*
* Méthode temporaire pour générer des scores sur les arêtes spectres-peptides. Elle m'a servi à réaliser certains tests sur le modèles
*/
void compute_score_function(std::vector<Peptide*> &peptides, std::vector<Spectrum*> &spectra, std::vector<Score*> &scores) {
	srand((unsigned int)time(NULL));
	for (Spectrum* spectrum : spectra) {
		std::size_t peptide_number;
		std::set<std::size_t> peptides_numbers = {};
		float scores_[3] = {0.26f ,0.75f ,0.99f};
		//float scores_[3] = { 0.4f ,0.7f ,0.9f };
		//float scores_[2] = { 0.2f ,0.8f };
		for (unsigned int i = 0; i < 3; ++i) {
			do {
				peptide_number = (std::size_t)(rand() % peptides.size());
			} while (peptides_numbers.contains(peptide_number));
			peptides_numbers.insert(peptide_number);
			scores.push_back(new Score(peptide_number, spectrum->Get_Id(), scores_[i]));
		}
	}
}

/*
* Il faut commencer par créer une variable de type Model pour pouvoir exécuter les tâches suivantes dessus :
* - Charger les protéines (Load_Proteins)
* - Construire les peptides théoriques (In_Silico_Digestion)
* - Construire les spectres théoriques (Build_Theoretical_Spectra)
* - Définir les probabilités sur les arêtes protéines-spectres (Define_Probabilities)
* - Charger des spectres expérimentaux ou générer des spectres simulés (Load_Spectra / Simulated_Sample)
* - Générer des scores sur les arêtes spectres-peptides (Load_Scores / Compute_Score / Compute_Score_SpecOMS)
* - Résoudre l'instance courante (Solve)
* - Afficher la solution (Print_Solution)
* 
* Chacune de ces étapes possède leur propre fichier source contenant leur documentation
*/
int main() {

	/*
	* Un exemple
	*/
	Model model;

	model.Load_Proteins("uniprot_humain_5_prot.fasta");
	std::cout << "proteins loaded : " << model.Number_Of_Proteins() << std::endl;

	model.In_Silico_Digestion("test_digestion", 0 ,40);
	std::cout << "proteins digested : " << model.Number_Of_Peptides() << std::endl;

	model.Build_Theoretical_Spectra();
	std::cout << "theoretical spectra built : " << model.Number_Of_Peptides() << std::endl;

	model.Define_Probabilities(std::string("test_result.csv"));
	// model.Define_Probabilities(1.0f);
	std::cout << "probabilities defined" << std::endl;
	// std::cout << model.Get_Peptide(207160).Get_Proteins().at(6113)[0] << std::endl;

	// model.Run_Tests_Sample_Data(10, 5, 5, 2, 2, "tests_results");
	
	// /*

	std::unordered_map<std::size_t, unsigned int> sample = model.Random_Sample(5, 5, 2, 2);
	// std::unordered_map<std::size_t, unsigned int> sample = {{1, 1}, {3, 3}, {4, 4}};
	model.Simulated_Sample(sample);
	// model.Load_Spectra("110616_yeast_ups_10fmol.ms2");
	std::cout << "spectra loaded : " << model.Number_Of_Spectra() << std::endl;

	// model.Compute_Score(1);
	model.Compute_Score_SpecOMS();
	std::cout << "scores computed : " << model.Number_Of_Scores() << std::endl;

	int time = model.Solve(0.5f, 0.5f);
	std::cout << "model solved in " << time << " milliseconds" << std::endl;

	// model.Print_Solution();

	// model.Save_Solution("test_save", true);

	model.Analyse_Solution(sample, "analyse");

	/**/

	return 0;

}
