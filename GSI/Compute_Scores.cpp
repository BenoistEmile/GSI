
#include "Model.h"

#include <algorithm>
#include <set>
#include <ctime>
//#include <chrono>
//#include <thread>
//std::this_thread::sleep_for(std::chrono::milliseconds(2000));

//__________________________________________________________________________________________________________

struct Noeud {
	Noeud* relative;
	Noeud* prec;
	std::size_t id;
	unsigned int counter;

	Noeud(Noeud* relative, Noeud* prec, std::size_t id) : relative(relative), prec(prec), id(id), counter(1) {};

	~Noeud() {};

	Noeud& operator++() {
		++counter;
		if (relative != nullptr) {
			++*relative;
		}
		return *this;
	}
};

//__________________________________________________________________________________________________________

void Model::Compute_Score(void (compute_score_function)(std::vector<Peptide*> &peptides, std::vector<Spectrum*> &spectra, std::vector<Score*> &scores)) {
	compute_score_function(peptides, spectra, scores);
}

void Model::Compute_Score(unsigned int randoms) {
	for (auto iter = spectra.begin(); iter != spectra.end(); ++iter) {
		//srand((unsigned int)time(NULL));
		if ((*iter)->Get_Origin() != nullptr) {
			scores.push_back(new Score((*iter)->Get_Origin()->peptide, (*iter)->Get_Id(), 0.0));
			std::size_t peptide_number = (*iter)->Get_Origin()->peptide;
			std::set<std::size_t> peptides_numbers = { peptide_number };
			for (unsigned int _ = 0; _ < randoms; ++_) {
				do {
					peptide_number = (std::size_t)(rand() % peptides.size());
				} while (peptides_numbers.contains(peptide_number));
				peptides_numbers.insert(peptide_number);
				scores.push_back(new Score(peptide_number, (*iter)->Get_Id(), 0.0));
			}
		}
		else {
			std::cout << "ERROR : this function can only be used on simulated spectrum." << std::endl;
		}
	}
}

void Model::Compute_Score_SpecOMS(unsigned int minimum_number_of_masses, unsigned int maximum_number_of_masses, int accuracy, unsigned int number_of_copies, unsigned int threshold, unsigned int maximum_number_of_edges) {

	const std::size_t number_of_spectra = spectra.size();
	const std::size_t number_of_peptides = peptides.size();
	std::unordered_map<double, std::vector<std::size_t>*> buckets_cluster;

	//std::cout << "etape 1" << std::endl;
	/*
	* Filling the bucket cluster with theoretical spectra
	*/
	const std::vector<Pic*>* pics;
	for (auto iter = peptides.begin(); iter != peptides.end(); ++iter) {
		pics = (*iter)->Get_Pics();//(*iter)->Build_Spectrum();
		for (auto iter2 = (*pics).begin(); iter2 != (*pics).end(); ++iter2) {
			double rounded_pic = Round_Precision((*iter2)->mass, accuracy);
			if (buckets_cluster.contains(rounded_pic)) {
				buckets_cluster[rounded_pic]->push_back((*iter)->Get_Id());
			}
			else {
				buckets_cluster[rounded_pic] = new std::vector<std::size_t>{ (*iter)->Get_Id() };
			}
			delete* iter2;
		}
		delete pics;
	}

	//std::cout << "etape 2" << std::endl;
	/*
	* Filling the bucket cluster with experimental spectra
	*/
	const double min_accuracy = pow(10, -accuracy);
	for (Spectrum* spectrum : spectra) {
		std::vector<Pic*>* selected_pics = spectrum->Filter1(minimum_number_of_masses, maximum_number_of_masses, accuracy ,number_of_copies);
		if (spectrum->Is_Simulated()) {
			for (Pic* pic : *selected_pics) {
				if (buckets_cluster.contains(pic->mass)) {
					buckets_cluster[pic->mass]->push_back(number_of_peptides + spectrum->Get_Id());
				}
				delete pic;
			}
		}
		else {
			for (Pic* pic : *selected_pics) {
				std::cout << pic->mass << std::endl;
				if (buckets_cluster.contains(pic->mass)) {
					buckets_cluster[pic->mass]->push_back(number_of_peptides + spectrum->Get_Id());
				}
				double space = 0;
				for (unsigned int _ = 1; _ <= number_of_copies; ++_) {
					space += min_accuracy;
					double mass1 = Round_Precision(pic->mass + space, accuracy);
					double mass2 = Round_Precision(pic->mass - space, accuracy);
					if (buckets_cluster.contains(mass1)) {
						buckets_cluster[mass1]->push_back(number_of_peptides + spectrum->Get_Id());
					}
					if (buckets_cluster.contains(mass2)) {
						buckets_cluster[mass2]->push_back(number_of_peptides + spectrum->Get_Id());
					}
				}
				delete pic;
			}
		}
		delete selected_pics;
	}

	//std::cout << "etape 3" << std::endl;
	/*
	* Sorting and cleaning the bucket cluster
	*/
	std::vector<std::vector<std::size_t>*> sorted_buckets;
	for (auto iter = buckets_cluster.begin(); iter != buckets_cluster.end(); ++iter) {
		if (iter->second->back() >= peptides.size()) {
			sorted_buckets.insert(std::upper_bound(sorted_buckets.begin() ,sorted_buckets.end(), iter->second, [](std::vector<std::size_t>* a, std::vector<std::size_t>* b) -> bool {return (*a) < (*b); }), iter->second);
		}
		else {
			delete iter->second;
		}
	}

	//std::cout << "etape 4 : " << sorted_buckets.size() << std::endl;
	/*
	* Creation of the specTree
	*/
	Noeud** last = new Noeud * [number_of_spectra];
	for (unsigned int i = 0; i < number_of_spectra; ++i) {
		last[i] = nullptr;
	}
	Noeud* deletion = nullptr;
	Noeud* prec = nullptr;
	Noeud* curr = nullptr;
	std::size_t i = 0;
	if (sorted_buckets.size() > 0){
		while (sorted_buckets[0]->at(i) < number_of_peptides){
			curr = new Noeud(prec ,deletion, sorted_buckets[0]->at(i));
			prec = curr;
			deletion = curr;
			++i;
		}
		for (; i < sorted_buckets[0]->size(); ++i) {
			curr = new Noeud(prec, last[sorted_buckets[0]->at(i) - number_of_peptides], sorted_buckets[0]->at(i));
			prec = curr;
			last[sorted_buckets[0]->at(i) - number_of_peptides] = curr;
		}
		for (i = 1; i < sorted_buckets.size(); ++i) {
			unsigned int j = 0;
			while (j < sorted_buckets[i - 1]->size() && sorted_buckets[i - 1]->at(j) == sorted_buckets[i]->at(j)) {
				++j;
			}
			if (j == 0) {
				prec = nullptr;
			}
			else {
				Noeud* start;
				if (j == sorted_buckets[i - 1]->size()) {
					start = last[sorted_buckets[i - 1]->back() - number_of_peptides];
				}
				else if (sorted_buckets[i - 1]->at(j) < number_of_peptides) {
					unsigned int k = j + 1;
					while (sorted_buckets[i - 1]->at(k) < number_of_peptides) {
						++k;
					}
					start = last[sorted_buckets[i - 1]->at(k) - number_of_peptides];
					for (; k >= j && start != nullptr; --k, start = start->relative);
				}
				else {
					start = last[sorted_buckets[i - 1]->at(j) - number_of_peptides];
					if (start != nullptr) { // pour eviter le warning
						start = start->relative;
					}
				}
				++*start; // faire gaffe avec la recurcivite !
				prec = start;
			}
			if (j != sorted_buckets[i]->size()) {
				while (sorted_buckets[i]->at(j) < number_of_peptides) {
					curr = new Noeud(prec, deletion, sorted_buckets[i]->at(j));
					prec = curr;
					deletion = curr;
					++j;
				}
				for (; j < sorted_buckets[i]->size(); ++j) {
					curr = new Noeud(prec, last[sorted_buckets[i]->at(j) - number_of_peptides], sorted_buckets[i]->at(j));
					prec = curr;
					last[sorted_buckets[i]->at(j) - number_of_peptides] = curr;
				}
			}
		}

		//std::cout << "etape 5" << std::endl;
		/*
		* Extraction of pairs
		*/
		std::size_t k;
		unsigned int scores_sum;
		std::vector<unsigned int> total_scores;
		std::vector<std::size_t> pairs;
		total_scores.resize(peptides.size());
		pairs.reserve(peptides.size());
		Noeud* curr2;
		for (std::size_t i = 0; i < number_of_spectra; ++i) {
			curr = last[i];
			while (curr != nullptr) {
				curr2 = curr->relative;
				while (curr2->id >= number_of_peptides) {
					curr2 = curr2->relative;
				}
				while (curr2 != nullptr) {
					if (total_scores[curr2->id] == 0) {
						pairs.push_back(curr2->id);
					}
					total_scores[curr2->id] += curr->counter;
					curr2 = curr2->relative;
				}
				curr = curr->prec;
			}
			sort(pairs.begin(), pairs.end(), [total_scores](std::size_t i, std::size_t j)-> bool {return total_scores[i] > total_scores[j]; });
			k = 0;
			scores_sum = 0;
			while (k < pairs.size() && (k < maximum_number_of_edges || maximum_number_of_edges == 0) && total_scores[pairs[k]] >= threshold) {
				scores_sum += total_scores[pairs[k]];
				++k;
			}
		
			for (std::size_t l = k; l < pairs.size(); l++) {
				total_scores[pairs[l]] = 0;
			}
		
			if (k != 0) {
				--k;
				for (; k > 0; --k) {
					scores.push_back(new Score(pairs[k], i, 1.0 - ((double)total_scores[pairs[k]] / scores_sum)));
					total_scores[pairs[k]] = 0;
				}
				scores.push_back(new Score(pairs[0], i, 1.0 - ((double)total_scores[pairs[0]] / scores_sum)));
				total_scores[pairs[0]] = 0;
			}
			pairs.clear();
		}

		//std::cout << "etape 6" << std::endl;
		/*
		* Cleaning step
		*/
		prec = deletion;
		curr = deletion;
		while (curr != nullptr) {
			curr = prec->prec;
			delete prec;
			prec = curr;
		}
		for (std::size_t i = 0; i < number_of_spectra; ++i) {
			prec = last[i];
			curr = last[i];
			while (curr != nullptr) {
				curr = prec->prec;
				delete prec;
				prec = curr;
			}
		}
		for (auto iter = sorted_buckets.begin(); iter != sorted_buckets.end(); ++iter) {
			delete* iter;
		}
	}
	delete[] last;

	/*
	std::cout << "___________" << std::endl;
	
	for (auto iter = sorted_buckets.begin(); iter != sorted_buckets.end(); ++iter) {
		std::cout << "[";
		for (auto iter2 = (*iter)->begin(); iter2 != (*iter)->end(); ++iter2) {
			std::cout << (*iter2) << " ,";
		}
		std::cout << "]" << std::endl;
		delete *iter;
	}
	*/


}