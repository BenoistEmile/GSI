
#include "Model.h"
#include <unordered_map>
#include <set>

//__________________________________________________________________________________________________________

/*
* M�thode temporaire pour g�n�rer des scores sur les ar�tes spectres-peptides. Elle m'a servi � r�aliser certains tests sur le mod�les
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
* Il faut commencer par cr�er une variable de type Model pour pouvoir ex�cuter les t�ches suivantes dessus :
* - Charger les prot�ines (Load_Proteins)
* - Construire les peptides th�oriques (In_Silico_Digestion)
* - Construire les spectres th�oriques (Build_Theoretical_Spectra)
* - D�finir les probabilit�s sur les ar�tes prot�ines-spectres (Define_Probabilities)
* - Charger des spectres exp�rimentaux ou g�n�rer des spectres simul�s (Load_Spectra / Simulated_Sample)
* - G�n�rer des scores sur les ar�tes spectres-peptides (Load_Scores / Compute_Score / Compute_Score_SpecOMS)
* - R�soudre l'instance courante (Solve)
* - Afficher la solution (Print_Solution)
* 
* Chacune de ces �tapes poss�de leur propre fichier source contenant leur documentation
*/
int main() {

	/*
	* Un exemple
	*/
	Model model;

	model.Load_Proteins("uniprot_humain_moitier.fasta");
	std::cout << "proteins loaded : " << model.Number_Of_Proteins() << std::endl;

	model.In_Silico_Digestion(0 ,0);
	std::cout << "proteins digested : " << model.Number_Of_Peptides() << std::endl;

	model.Build_Theoretical_Spectra();
	std::cout << "theoretical spectra built : " << model.Number_Of_Peptides() << std::endl;

	model.Define_Probabilities(0.2f);
	std::cout << "probabilities defined" << std::endl;

	std::unordered_map<std::size_t, unsigned int> sample = {{2 , 5} ,{3 , 2}};
	model.Simulated_Sample(sample);
	std::cout << "spectra loaded : " << model.Number_Of_Spectra() << std::endl;

	model.Compute_Score(1);
	std::cout << "scores computed : " << model.Number_Of_Scores() << std::endl;

	model.Solve(0.5f, 0.5f);
	std::cout << "model solved" << std::endl;

	model.Print_Solution();

	return 0;

}