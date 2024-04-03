
#include "Model.h"

//__________________________________________________________________________________________________________

void Model::Build_Theoretical_Spectra(std::unordered_map<char, double> modifications) {
	for (Peptide* peptide : peptides) {
		peptide->Build_Spectrum(modifications);
	}
}