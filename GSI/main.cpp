
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

	// model.Load_Proteins("uniprot_humain_moitie.fasta");
	// model.Load_Proteins("ups1-ups2-sequences.fasta");
	// model.Load_Proteins("c_albicans+ups.fasta");
	// model.Load_Proteins_Accession("yeast+ups1.fasta");
	model.Load_Proteins_Accession("Sprot_2024-02-05.fasta");
	std::cout << "proteins loaded : " << model.Number_Of_Proteins() << std::endl;

	model.Peptide_Detectability(2, 0.00, 7, 25, true); // if SpecOMS is used, L2I must be true
	
	// std::string digestion_file_name = "yeast_10fmol";
	// std::cout << std::string(digestion_file_name + "_result.csv") << std::endl;
	// model.In_Silico_Digestion(digestion_file_name);
	// model.In_Silico_Digestion_2(digestion_file_name);
	
	// model.Peptide_detectability("Dby_Deep", digestion_file_name);
	std::cout << "Peptide detectability computed" << std::endl;
	std::cout << "Peptides digested : " << model.Number_Of_Peptides() << std::endl;

	// std::ofstream peptides_fasta(std::filesystem::current_path().generic_string() + "/peptides_yeast+ups1.fasta");
	// for (int i = 0; i < model.Number_Of_Peptides(); i++) {
	// 	peptides_fasta << ">" << i << std::endl << model.Get_Peptide(i).Get_Sequence() << std::endl << std::endl;
	// }
	// peptides_fasta.close();

	// model.Build_Theoretical_Spectra();
	// std::cout << "theoretical spectra built : " << model.Number_Of_Peptides() << std::endl;

	// model.Define_Probabilities(std::string(digestion_file_name + "_result.csv"), 0.0);
	// model.Define_Probabilities_2(std::string(digestion_file_name + "_result.csv"), 0.0);
	// model.Define_Probabilities(0.2f);
	// std::cout << "probabilities defined" << std::endl;
	// std::cout << model.Get_Peptide(207160).Get_Proteins().at(6113)[0] << std::endl;

	// model.Run_Tests_Sample_Data(10, 5, 5, 2, 2, "tests_results");

	// std::unordered_map<std::size_t, unsigned int> sample = model.Random_Sample(100, 100, 1, 50);
	// std::unordered_map<std::size_t, unsigned int> sample = {{1, 1}, {3, 3}, {4, 4}};
	// model.Simulated_Sample(sample);
	// model.Load_Spectra("110616_yeast_ups_10fmol.ms2", 60);
	model.Load_Spectra("QX001127_OVA.mgf", 50);
	std::cout << "spectra loaded : " << model.Number_Of_Spectra() << std::endl;

	// model.Compute_Score(1);
	// model.Compute_Score_SpecOMS(0U, 99999U, 2, 2U, 7U, 2U);
	model.Load_Scores_SpecOMS("specoms_output_OVA.csv");
	// model.Load_Scores_Prospect("110618_yeast_ups_50fmol_r1_peptides.csv", 6, 25, 4);
	// model.Load_Scores_XTandem("QX002755_Hela.csv");
	std::cout << "scores computed : " << model.Number_Of_Scores() << std::endl;

	// std::set<std::tuple<float, float>> psi_values = {{1, 1}, {1, 10}, {1, 100}};
	// model.Test_Psi_Values(psi_values, output_file, "results_yeast+ups1");

	std::set<std::pair<float, float>> params_set = {{1, 1}, {1, 10}, {1, 100}, {1, 1000}, {10, 1}, {100, 1}, {1000, 1}, {0, 1}, {0.2, 0.8}, {0.4, 0.6}, {0.5, 0.5}, {0.6, 0.4}, {0.8, 0.2}, {1, 0}};
	for (auto& params : params_set) {

		std::cout << params.first << ", " << params.second << std::endl;

	// unsigned int count;
	// for (float i = 0.0; i < 1.001; i += 0.1) {
		model.Solve(params.first, params.second);
	// model.Solve(1, 10);
		std::cout << "model solved" << std::endl;

		model.Print_Solution();

		// std::string file_name = "HeLa_SpecOMS_2_8_0_1_10_0.00";
		std::string file_name = "OVA_SpecOMS_1_8_0_" + fmt::format("{0:.1f}", params.first) + "_" + fmt::format("{0:.1f}", params.second) + "_0.00";

		std::filesystem::path file_path = std::filesystem::current_path() / "models" / ("upper_edges_" + file_name + ".csv");
		std::ofstream upper_edges_file(file_path);
		upper_edges_file << "accession,protein_id,peptide_id,Prob" << std::endl;
		for (int i = 0; i < model.Number_Of_Peptides(); i++) {
			for (auto& protein : model.Get_Peptide(i).Get_Proteins()) {
				for (auto& edge : std::get<1>(protein)) {
					upper_edges_file << model.Get_Protein(std::get<0>(protein)).Get_Accession() << "," << model.Get_Protein(std::get<0>(protein)).Get_Id() << "," << model.Get_Peptide(i).Get_Id() << "," << edge << std::endl;
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
		model.Save_Solution("results_" + file_name, true, true, true, true);
		model.Clear(false, false, false, false, true);
	}
	// }
	// count = 0;
	// for (auto iter = model.Get_Solution().Get_abundances().begin(); iter != model.Get_Solution().Get_abundances().end(); iter++) {
	// 	if (iter->first <= 47) {
	// 		count++;
	// 	}
	// }
	// output_file << count << "/48 proteins from ups identified." << std::endl << std::endl;
	// model.Clear(false, false, false, false, true);
	// }

	// model.Print_Solution(false);
	
	// output_file.close();

	// model.Run_Test("yeast_10fmol_nonoise", 1, 10, 13, 4, 0, false);

	// for (int iter = 1; iter <= 5; iter++) {
	// 	for (auto& iter2 : {0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95}) {
	// 		Model model;
	// 		std::cout << std::to_string(iter) << std::endl;
	// 		model.Run_Test_Synthetic_Data("test_synth" + std::to_string(iter), "yeast+ups1.fasta", 1, 0.99, iter2, 1, 1, 10);
	// 		model.Clear();
	// 	}
	// }

	// std::set<std::tuple<float, float, unsigned int, unsigned int, float>> parameters = {{1, 10, 900, 4, 0.0}, {1, 10, 900, 10, 0.0}};
	// model.Run_Multiple_Tests(parameters, "yeast_10fmol");

	/**/

	return 0;

}
