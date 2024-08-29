#include "Model.h"
#include <unordered_map>
#include <set>
#include "fmt/format.h"

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

	// /*

	// std::ofstream output_file = model.Open_Output_File("result");

	// model.Load_Proteins_Accession("Sprot_Gallus_gallus_2021_04_13.fasta");
	model.Load_Proteins_Accession("Sprot_2024-02-05.fasta");
	// model.Load_Proteins_Accession("Sprot_Tax9606_human_2023_03_22.fasta");
	std::cout << "proteins loaded : " << model.Number_Of_Proteins() << std::endl;

	model.Peptide_Detectability(2, 0.00, 7, 25, false); // if SpecOMS is used, L2I must be true

	std::cout << "Peptide detectability computed" << std::endl;
	std::cout << "Peptides digested : " << model.Number_Of_Peptides() << std::endl;

	// model.Load_Spectra("QX001127_OVA.mgf", 50);
	model.Load_Spectra("QX002755_HeLa.mgf", 50);
	std::cout << "spectra loaded : " << model.Number_Of_Spectra() << std::endl;

	// model.Load_Scores_SpecOMS("specoms_output_HeLa.csv");
	// model.Load_Scores_XTandem("QX002755_Hela-WithAccess-b.csv", 0.0001, "QX002755_Hela-WithAccess-b_proteins.csv");
	model.Load_Scores_XTandem("QX002755_Hela-classical-Evalue-param_Sprot-2024-02-05.csv", 0.0, "QX002755_Hela-classical-Evalue-param_Sprot-2024-02-05_proteins.csv");
	std::cout << "scores computed : " << model.Number_Of_Scores() << std::endl;

	// model.Pre_Solve(0.9);

	std::set<std::pair<float, float>> params_set = {{1, 1}};
	for (auto& params : params_set) {

		std::cout << std::get<0>(params) << ", " << std::get<1>(params) << std::endl;

	// unsigned int count;
	// for (float i = 0.0; i < 1.001; i += 0.1) {
		auto start = std::chrono::high_resolution_clock::now();
		model.Solve(std::get<0>(params), std::get<1>(params));
		auto end = std::chrono::high_resolution_clock::now();
		auto duration_tot = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << duration_tot << std::endl;
	// model.Solve(1, 10);
		std::cout << "model solved" << std::endl;

		model.Print_Solution();

	// std::string file_name = "ova_SpecOMS_test_2_2_8_0_1.0_10.0_0.00";
		// std::string file_name = "OVA_test_2_8_0_" + fmt::format("{0:.1f}", params.first) + "_" + fmt::format("{0:.1f}", params.second) + "_0.00";
		std::string file_name = "HeLa_full_optimist_2_8_0_" + fmt::format("{0:.1f}", std::get<0>(params)) + "_" + fmt::format("{0:.1f}", std::get<1>(params)) + "_0.00";

		model.Save_Solution("results_" + file_name, true, true, true, true);

		std::filesystem::path file_path = std::filesystem::current_path() / "models" / ("upper_edges_" + file_name + ".csv");
		std::ofstream upper_edges_file(file_path);
		upper_edges_file << "accession,protein_id,peptide_id,Prob,rank" << std::endl;
		for (int i = 0; i < model.Number_Of_Peptides(); i++) {
			for (auto& protein : model.Get_Peptide(i).Get_Proteins()) {
				// for (auto& edge : std::get<1>(protein)) {
				for (int j = 0; j < std::get<1>(protein).size(); j++) {
					upper_edges_file << model.Get_Protein(std::get<0>(protein)).Get_Accession() << "," << model.Get_Protein(std::get<0>(protein)).Get_Id() << "," << model.Get_Peptide(i).Get_Id() << "," << std::get<1>(protein).at(j) << "," << j << std::endl;
				}
			}
		}
		upper_edges_file.close();
		file_path = std::filesystem::current_path() / "models" / ("lower_edges_" + file_name + ".csv");
		std::ofstream lower_edges_file(file_path);
		lower_edges_file << "Peptide,Spectrum,Score" << std::endl;
		for (std::size_t iter_score = 0; iter_score < model.Number_Of_Scores(); iter_score++) {
			Score score = model.Get_Score(iter_score);
			lower_edges_file << score.peptide << "," << score.spectrum << "," << score.score << std::endl;
		}
		lower_edges_file.close();

		// model.Save_Solution(output_file, false, false, true);
		// model.Save_Solution("results_" + file_name, true, true, true, true);
		model.Clear(false, false, false, false, true);
		std::cout << "here" << std::endl;
	}

	return 0;

}