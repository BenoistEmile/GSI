
#include "Model.h"

#include "ilcplex/ilocplex.h"
#include <unordered_map>
#include <set>

//__________________________________________________________________________________________________________

int Model::Solve(const float psi1, const float psi2) {

	std::clock_t c_start = std::clock();

	solution.Clear();

	const size_t nb_proteins = proteins.size();
	const size_t nb_peptides = peptides.size();
	const size_t nb_spectra = spectra.size();
	const size_t nb_scores = scores.size();

	std::vector<std::vector<std::size_t>*> peptides_spectra; // arêtes à calculer pour chaque peptide utile
	peptides_spectra.reserve(nb_peptides);
	std::vector<std::vector<std::size_t>*> spectra_peptides; // même chose pour spectres
	spectra_peptides.reserve(nb_spectra);
	std::vector<unsigned int> peptides_sure;
	peptides_sure.reserve(nb_peptides); // nombre d'arêtes sûres pour chaque peptide utile
	std::vector<const Score*> useful_scores;
	useful_scores.reserve(nb_scores); // garde que arêtes utiles
	std::unordered_map<std::size_t, unsigned int> edge_per_spectrum; // nombre d'arêtes pour  chaque spectre

	for (auto edge : scores) {
		if (edge_per_spectrum.contains(edge->spectrum)) {
			edge_per_spectrum[edge->spectrum] += 1;
		}
		else {
			edge_per_spectrum[edge->spectrum] = 1;
		}
	}

	std::unordered_map<std::size_t, std::size_t> peptides_index; // ancien id, nouvel id
	std::unordered_map<std::size_t, std::size_t> spectra_index; // mm chose

	for (auto edge : scores) {
		if (edge_per_spectrum[edge->spectrum] > 1) {
			if (!spectra_index.contains(edge->spectrum)) {
				spectra_index[edge->spectrum] = spectra_peptides.size();
				spectra_peptides.push_back(new std::vector<std::size_t>{ useful_scores.size() });
			}
			else {
				spectra_peptides[spectra_index[edge->spectrum]]->push_back(useful_scores.size());
			}
			if (peptides_index.contains(edge->peptide)) {
				peptides_spectra[peptides_index[edge->peptide]]->push_back(useful_scores.size()); //
			}
			else {
				peptides_index[edge->peptide] = peptides_spectra.size();
				peptides_spectra.push_back(new std::vector<std::size_t>{ useful_scores.size() }); //
				peptides_sure.push_back(0);
			}
			useful_scores.push_back(edge);
		}
		else {
			if (peptides_index.contains(edge->peptide)) {
				peptides_sure[peptides_index[edge->peptide]] += 1;
				solution.Add_Score(edge);
			}
			else {
				peptides_index[edge->peptide] = peptides_spectra.size();
				peptides_spectra.push_back(new std::vector<std::size_t>);
				peptides_sure.push_back(1);
				solution.Add_Score(edge);
			}
		}
	}

	std::unordered_map<std::size_t, std::size_t> proteins_index; // en faire un vecteur ?
	std::unordered_map<std::size_t, std::size_t> proteins_index_new_old;
	std::vector<std::size_t> useful_proteins; // identifiants d'origine
	bool leave;

	for (Protein* protein : proteins) {
		std::size_t i = 0;
		leave = false;
		while (i < protein->Get_Peptides().size() && !leave) {
			if (peptides_index.contains(protein->Get_Peptides()[i])) {
				leave = true;
			}
			i++;
		}
		if (leave) {
			proteins_index[protein->Get_Id()] = useful_proteins.size();
			proteins_index_new_old[useful_proteins.size()] = protein->Get_Id();
			useful_proteins.push_back(protein->Get_Id());
		}
	}

	std::vector<std::vector<std::tuple<std::size_t, float>>*> peptides_proteins(peptides_spectra.size());
	for (std::size_t i = 0; i < peptides_spectra.size(); ++i) {
		peptides_proteins[i] = new std::vector<std::tuple<std::size_t, float>>;
	}
	std::vector<std::vector<std::tuple<std::size_t, float>>*> useless_peptides_proteins;
	useless_peptides_proteins.reserve(peptides.size() - peptides_spectra.size());
	std::size_t counter = 0;
	bool founded;
	for (Peptide* peptide : peptides) {
		if (peptides_index.contains(peptide->Get_Id())) {
			for (auto &edges : peptide->Get_Proteins()) {
				for (float prob : edges.second) {
					peptides_proteins[peptides_index[peptide->Get_Id()]]->push_back(std::tuple<std::size_t, float>(proteins_index[edges.first], prob));
				}
			}
		}
		else {
			founded = false;
			for (auto& edges : peptide->Get_Proteins()) {
				if (proteins_index.contains(edges.first)) {
					for (float prob : edges.second) {
						if (!founded) {
							useless_peptides_proteins.push_back(new std::vector<std::tuple<std::size_t, float>>);
							counter++;
							founded = true;
						}
						useless_peptides_proteins[counter - 1]->push_back(std::tuple<std::size_t, float>(proteins_index[edges.first], prob));
					}
				}
			}
		}
	}

	std::size_t n = useful_proteins.size();
	std::size_t m1 = peptides_proteins.size();
	std::size_t m2 = useless_peptides_proteins.size();
	std::size_t l = spectra_peptides.size();
	std::size_t o = useful_scores.size();

#pragma region variables and model

	IloEnv env;
	IloModel model(env);

	IloNumVarArray Q = IloNumVarArray(env, n, 0, IloInfinity);
	IloNumVarArray Delta = IloNumVarArray(env, m1, 0, IloInfinity);
	IloBoolVarArray X = IloBoolVarArray(env, o);

#pragma endregion

#pragma region constraints

	for (std::size_t k = 0; k < l; ++k) {
		IloExpr constraintX(env);
		for (std::size_t h : *(spectra_peptides[k])) {
			constraintX += X[h];
		}
		model.add(constraintX == 1);
	}

	for (std::size_t j = 0; j < m1; ++j) {
		IloExpr constraintDelta1(env);
		IloExpr constraintDelta2(env);
		for (std::tuple<std::size_t, float> edge : (*peptides_proteins[j])) {
			constraintDelta1 += Q[std::get<0>(edge)] * std::get<1>(edge);
			constraintDelta2 -= Q[std::get<0>(edge)] * std::get<1>(edge);
		}
		for (std::size_t h : (*peptides_spectra[j])) {
			constraintDelta1 -= X[h];
			constraintDelta2 += X[h];
		}
		constraintDelta1 -= Delta[j] + peptides_sure[j];
		constraintDelta2 -= Delta[j] - peptides_sure[j];
		model.add(constraintDelta1 <= 0);
		model.add(constraintDelta2 <= 0);
	}

#pragma endregion

#pragma region objective

	IloExpr objective(env);

	for (std::size_t j = 0; j < m1; ++j) {
		objective += psi1 * Delta[j];
	}
	for (std::size_t j = 0; j < m2; ++j) {
		for (std::tuple<std::size_t, float> edge : (*useless_peptides_proteins[j])) {
			objective += psi1 * std::get<1>(edge) * Q[std::get<0>(edge)];
		}
	}
	for (std::size_t h = 0; h < o; ++h) {
		objective += psi2 * X[h] * useful_scores[h]->score;
	}

	model.add(IloMinimize(env, objective));

#pragma endregion

	IloCplex cplex(model);
	auto start = std::chrono::high_resolution_clock::now();
	cplex.solve();
	auto end = std::chrono::high_resolution_clock::now();
	auto duration_tot = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

#pragma region generate solution

	IloNumArray valuesQ(env);
	cplex.getValues(Q, valuesQ);
	for (std::size_t i = 0; i < n; ++i) {
		if (valuesQ[i] > 0) {
			solution.Add_Protein(proteins_index_new_old[i], (float)valuesQ[i]); // protein_index[i] est inversé : on donne le nouvel indice (solution = créer dès le début une autre map inversée)
		}
	}

	IloNumArray valuesX(env);
	cplex.getValues(X, valuesX);
	for (std::size_t h = 0; h < o; ++h) {
		if (valuesX[h] > 0.5) {
			solution.Add_Score(scores[h]); // même pb qu'avant (créer une map pour retrouver le bon indice)
		}
	}
#pragma endregion

	IloNumArray valuesD(env);
	cplex.getValues(Delta, valuesD);

	env.out() << "Solution status = " << cplex.getStatus() << std::endl;
	env.out() << "Solution value = " << cplex.getObjValue() << std::endl;

#pragma region delete
	for (auto edges : peptides_spectra) {
		delete edges;
	}
	for (auto edges : spectra_peptides) {
		delete edges;
	}
	for (auto edges : peptides_proteins) {
		delete edges;
	}
	for (auto edges : useless_peptides_proteins) {
		delete edges;
	}
#pragma endregion

	std::cout << (std::clock() - c_start) / CLOCKS_PER_SEC << " secondes" << std::endl;
	return duration_tot;
}

int Model::Solve(std::ofstream& output_file, float psi1, const float psi2) {
	int duration = this->Solve(psi1, psi2);
	output_file << "Solved the model in " << duration << " milliseconds" << std::endl;
	output_file << "Parameters : psi1 = " << psi1 << ", psi2 = " << psi2 << std:: endl << std::endl;
	return duration;
}