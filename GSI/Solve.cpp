
#include "Model.h"

#include "ilcplex/ilocplex.h"
#include <unordered_map>
#include <set>

//__________________________________________________________________________________________________________

int Model::Solve(const float psi1, const float psi2, const float Pmin, const float psi3) {

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
	std::vector<unsigned int> peptides_sure; // nombre d'arêtes sûres pour chaque peptide utile
	peptides_sure.reserve(nb_peptides);
	std::vector<const Score*> useful_scores; // garde que arêtes utiles
	useful_scores.reserve(nb_scores);
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
	std::unordered_map<std::size_t, std::size_t> spectra_index; // ancien id, nouvel id

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

	std::unordered_map<std::size_t, std::size_t> useless_peptides_index; // ancien id, nouvel id
	std::unordered_map<std::size_t, std::size_t> proteins_index; // en faire un vecteur ?
	std::unordered_map<std::size_t, std::size_t> proteins_index_new_old;
	std::vector<std::size_t> useful_proteins; // identifiants d'origine
	bool leave;
	std::vector<std::vector<std::size_t>*> proteins_useful_detectabilities;
	std::vector<float> useful_detectabilities;
	std::vector<std::vector<std::tuple<std::size_t, std::size_t>>*> peptides_proteins(peptides_spectra.size()); // new id protein, id detectability
	for (std::size_t i = 0; i < peptides_spectra.size(); ++i) {
		peptides_proteins[i] = new std::vector<std::tuple<std::size_t, std::size_t>>;
	}
	std::vector<float> useless_detectabilities;
	std::vector<std::vector<std::size_t>*> proteins_useless_detectabilities;
	std::vector<std::vector<std::tuple<std::size_t, std::size_t>>*> useless_peptides_proteins;
	useless_peptides_proteins.reserve(peptides.size() - peptides_spectra.size());
	// std::size_t counter = 0;
	// bool founded;

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
			proteins_useless_detectabilities.push_back(new std::vector<std::size_t>);
			for (std::size_t peptide_id: protein->Get_Peptides()) {
				if (peptides_index.contains(peptide_id)) {
					for (auto& iter_pep: this->Get_Peptide(peptide_id).Get_Proteins()) {
						if (iter_pep.first == protein->Get_Id()) {
							for (auto& iter_detect : iter_pep.second) {
								if (proteins_useful_detectabilities.size() < useful_proteins.size()) {
									proteins_useful_detectabilities.push_back(new std::vector<std::size_t>{useful_detectabilities.size()});
								}
								else {
									proteins_useful_detectabilities[proteins_index[protein->Get_Id()]]->push_back(useful_detectabilities.size());
								}
								peptides_proteins[peptides_index[peptide_id]]->push_back(std::tuple<std::size_t, std::size_t>(proteins_index[protein->Get_Id()], useful_detectabilities.size()));
								useful_detectabilities.push_back(iter_detect);
							}
						}
					}
				}
				else {
					if (!useless_peptides_index.contains(peptide_id)) {
						useless_peptides_index[peptide_id] = useless_peptides_proteins.size();
						useless_peptides_proteins.push_back(new std::vector<std::tuple<std::size_t, std::size_t>>);
						// counter++;
					}
					for (auto& iter_pep: this->Get_Peptide(peptide_id).Get_Proteins()) {
						if (iter_pep.first == protein->Get_Id()) {
							for (auto& iter_detect: iter_pep.second) {
								proteins_useless_detectabilities[proteins_index[protein->Get_Id()]]->push_back(useless_detectabilities.size());
								useless_peptides_proteins[useless_peptides_index[peptide_id]]->push_back(std::tuple<std::size_t, std::size_t>(proteins_index[protein->Get_Id()], useless_detectabilities.size()));
								useless_detectabilities.push_back(iter_detect);
							}
						}
					}
				}
			}
		}
	}

	// for (Peptide* peptide : peptides) {
	// 	if (!(peptides_index.contains(peptide->Get_Id()))) {
	// 		founded = false;
	// 		for (auto& edges : peptide->Get_Proteins()) {
	// 			if (proteins_index.contains(edges.first)) {
	// 				for (float prob : edges.second) {
	// 					if (!founded) {
	// 						useless_peptides_proteins.push_back(new std::vector<std::tuple<std::size_t, std::size_t>>);
	// 						counter++;
	// 						founded = true;
	// 					}
	// 					useless_peptides_proteins[counter - 1]->push_back(std::tuple<std::size_t, float>(proteins_index[edges.first], prob));
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	std::size_t n = useful_proteins.size();
	std::size_t m1 = peptides_proteins.size();
	std::size_t m2 = useless_peptides_proteins.size();
	std::size_t l = spectra_peptides.size();
	std::size_t o = useful_scores.size();
	std::size_t p = useful_detectabilities.size();
	std::size_t q = useless_detectabilities.size();

#pragma region variables and model

	IloEnv env;
	IloModel model(env);

	IloNumVarArray Q = IloNumVarArray(env, n, 0, IloInfinity);
	IloNumVarArray Q1 = IloNumVarArray(env, p, 0, IloInfinity);
	IloNumVarArray Q2 = IloNumVarArray(env, q, 0, IloInfinity);
	IloNumVarArray Delta = IloNumVarArray(env, m1, 0, IloInfinity);
	IloBoolVarArray X = IloBoolVarArray(env, o);
	IloBoolVarArray Y1 = IloBoolVarArray(env, p);
	IloBoolVarArray Y2 = IloBoolVarArray(env, q);

#pragma endregion

#pragma region constraints

	int M = 100000;

	for (std::size_t k = 0; k < l; ++k) {
		IloExpr constraintX(env);
		for (std::size_t h : *(spectra_peptides[k])) {
			constraintX += X[h];
		}
		model.add(constraintX == 1);
	}

	for (std::size_t i = 0; i < n; i++) {
		std::size_t n_detect = proteins_useful_detectabilities[i]->size() + proteins_useless_detectabilities[i]->size();
		IloExpr constraintY(env);
		for (std::size_t l: *(proteins_useful_detectabilities[i])) {
			constraintY += Y1[l] / n_detect;
			// IloExpr constraintQ1(env);
			// IloExpr constraintQ1 = Q1[l];
			// constraintQ1 += Q1[l];
			model.add(Q1[l] <= (Q[i] * useful_detectabilities[l]) + (M * (1-Y1[l])));
			model.add(Q1[l] >= (Q[i] * useful_detectabilities[l]) - (M * (1-Y1[l])));
			model.add(Q1[l] <= M * Y1[l]);
			model.add(Q1[l] >= -(M * Y1[l]));
			// model.add(IloIfThen(env, Y1[l] == 0, Q1[l] == 0));
			// model.add(IloIfThen(env, Y1[l] == 1, Q1[l] == Q[i] * useful_detectabilities[l]));
			// model.add(Q1[l] == Q[i] * useful_detectabilities[l]);
		}
		for (std::size_t l: *(proteins_useless_detectabilities[i])) {
			constraintY += Y2[l] / n_detect;
			// IloExpr constraintQ2(env);
			// IloExpr constraintQ2 = Q2[l];
			// constraintQ2 += Q2[l];
			model.add(Q2[l] <= (Q[i] * useless_detectabilities[l]) + (M * (1-Y2[l])));
			model.add(Q2[l] >= (Q[i] * useless_detectabilities[l]) - (M * (1-Y2[l])));
			model.add(Q2[l] <= M * Y2[l]);
			model.add(Q2[l] >= -(M * Y2[l]));
			// model.add(IloIfThen(env, Y2[l] == 0, Q2[l] == 0));
			// model.add(IloIfThen(env, Y2[l] == 1, Q2[l] == Q[i] * useless_detectabilities[l]));
			// model.add(Q2[l] == Q[i] * useless_detectabilities[l]);
		}
		model.add(constraintY >= Pmin);
	}


	for (std::size_t j = 0; j < m1; ++j) {
		IloExpr constraintDelta1(env);
		IloExpr constraintDelta2(env);
		for (std::tuple<std::size_t, std::size_t> edge : (*peptides_proteins[j])) {
			// constraintDelta1 += Q[std::get<0>(edge)] * useful_detectabilities[std::get<1>(edge)] * Y1[std::get<1>(edge)];
			// constraintDelta2 -= Q[std::get<0>(edge)] * useful_detectabilities[std::get<1>(edge)] * Y1[std::get<1>(edge)];
			// constraintDelta1 += Q[std::get<0>(edge)] * useful_detectabilities[std::get<1>(edge)];
			// constraintDelta2 -= Q[std::get<0>(edge)] * useful_detectabilities[std::get<1>(edge)];
			constraintDelta1 += Q1[std::get<1>(edge)];
			constraintDelta2 -= Q1[std::get<1>(edge)];
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
		for (std::tuple<std::size_t, std::size_t> edge : (*useless_peptides_proteins[j])) {
			// objective += psi1 * Q[std::get<0>(edge)] * useless_detectabilities[std::get<1>(edge)];// * Y2[std::get<1>(edge)];
			objective += psi1 * Q2[std::get<1>(edge)];
		}
	}
	for (std::size_t h = 0; h < o; ++h) {
		objective += psi2 * X[h] * useful_scores[h]->score;
	}
	for (std::size_t l = 0; l < p; l++) {
		objective += psi3 * (1 - Y1[l]) * useful_detectabilities[l];
	}
	for (std::size_t l = 0; l < p; l++) {
		objective += psi3 * (1 - Y2[l]) * useless_detectabilities[l];
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
			solution.Add_Score(useful_scores[h]); // même pb qu'avant (créer une map pour retrouver le bon indice)
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