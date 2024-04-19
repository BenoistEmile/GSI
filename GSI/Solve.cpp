
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

	//std::vector<std::vector<std::size_t>*> peptides_spectra_first; //
	//peptides_spectra_first.reserve(nb_peptides); //
	std::vector<std::vector<std::size_t>*> peptides_spectra;
	peptides_spectra.reserve(nb_peptides);
	std::vector<std::vector<std::size_t>*> spectra_peptides;
	spectra_peptides.reserve(nb_spectra);
	std::vector<unsigned int> peptides_sure;
	peptides_sure.reserve(nb_peptides);
	std::vector<const Score*> useful_scores;
	useful_scores.reserve(nb_scores);
	std::unordered_map<std::size_t, unsigned int> edge_per_spectrum;

	for (auto edge : scores) {
		if (edge_per_spectrum.contains(edge->spectrum)) {
			edge_per_spectrum[edge->spectrum] += 1;
		}
		else {
			edge_per_spectrum[edge->spectrum] = 1;
		}
	}

	std::unordered_map<std::size_t, std::size_t> peptides_index;
	std::unordered_map<std::size_t, std::size_t> spectra_index;

	/*
	for (auto edge : scores) {
		if (edge_per_spectrum[edge->spectrum] > 1) {
			useful_scores.push_back(*edge);
			if (peptides_index.contains(edge->peptide)) {
				peptides_spectra[peptides_index[edge->peptide]]->push_back(counter);
			}
			else {
				peptides_index[edge->peptide] = peptides_spectra.size();
				peptides_spectra.push_back(new std::vector<std::size_t>{ counter });
				peptides_sure.push_back(0);
			}
			if (spectra_index.contains(edge->spectrum)) {
				spectra_peptides[spectra_index[edge->spectrum]]->push_back(counter);
			}
			else {
				spectra_index[edge->spectrum] = spectra_peptides.size();
				spectra_peptides.push_back(new std::vector<std::size_t>{ counter });
			}
			counter++;
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
	*/
	bool first; // a enlever
	//unsigned int number_of_first_edge = 0; //
	for (auto edge : scores) {
		if (edge_per_spectrum[edge->spectrum] > 1) {
			first = !spectra_index.contains(edge->spectrum);
			if (first) {
				spectra_index[edge->spectrum] = spectra_peptides.size();
				spectra_peptides.push_back(new std::vector<std::size_t>{ useful_scores.size() });
				//++number_of_first_edge; // 
			}
			else {
				spectra_peptides[spectra_index[edge->spectrum]]->push_back(useful_scores.size());
			}
			if (peptides_index.contains(edge->peptide)) {
				/*
				if (first) { //
					peptides_spectra_first[peptides_index[edge->peptide]]->push_back(useful_scores.size());
				}
				else {
					peptides_spectra[peptides_index[edge->peptide]]->push_back(useful_scores.size());
				}
				*/
				peptides_spectra[peptides_index[edge->peptide]]->push_back(useful_scores.size()); //
			}
			else {
				peptides_index[edge->peptide] = peptides_spectra.size();
				/*
				if (first) { //
					peptides_spectra_first.push_back(new std::vector<std::size_t>{ useful_scores.size() });
					peptides_spectra.push_back(new std::vector<std::size_t>);
				}
				else {
					peptides_spectra.push_back(new std::vector<std::size_t>{ useful_scores.size() });
					peptides_spectra_first.push_back(new std::vector<std::size_t>);
				}
				*/
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
				//peptides_spectra_first.push_back(new std::vector<std::size_t>); //
				peptides_spectra.push_back(new std::vector<std::size_t>);
				peptides_sure.push_back(1);
				solution.Add_Score(edge);
			}
		}
	}

	std::unordered_map<std::size_t, std::size_t> proteins_index; // en faire un vecteur ?
	std::vector<std::size_t> useful_proteins;
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
			useful_proteins.push_back(protein->Get_Id());
		}
	}

	std::vector<std::vector<std::tuple<std::size_t, float>>*> peptides_proteins(peptides_spectra.size());
	for (std::size_t i = 0; i < peptides_spectra.size(); ++i) {
		peptides_proteins[i] = new std::vector<std::tuple<std::size_t, float>>;
	}
	std::vector<std::vector<std::tuple<std::size_t, float>>*> useless_peptides_proteins(peptides.size() - peptides_spectra.size(), new std::vector<std::tuple<std::size_t, float>>);
	for (std::size_t i = 0; i < peptides.size() - peptides_spectra.size(); ++i) {
		useless_peptides_proteins[i] = new std::vector<std::tuple<std::size_t, float>>;
	}
	std::size_t counter = 0;
	bool founded;
	for (Peptide* peptide : peptides) {
		if (peptides_index.contains(peptide->Get_Id())) {
			for (auto &edges : peptide->Get_Proteins()) {
				if (proteins_index.contains(edges.first)) { // Inutile ??
					for (float prob : edges.second) {
						peptides_proteins[peptides_index[peptide->Get_Id()]]->push_back(std::tuple<std::size_t, float>(proteins_index[edges.first], prob));
					}
				}
			}
		}
		else {
			founded = false;
			for (auto& edges : peptide->Get_Proteins()) {
				if (proteins_index.contains(edges.first)) {
					for (float prob : edges.second) {
						if (!founded) {
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

	/*
	//edge_per_spectrum : associe un spectre a un nombre d'aretes															  std::unordered_map<std::size_t, unsigned int>
	std::cout << "________edge_per_spectrum________" << std::endl;
	for (auto& couple : edge_per_spectrum) {
		std::cout << couple.first << " -> " << couple.second << std::endl;
	}
	//peptide_spectra : indic� par les peptides. B[peptide] = liste des indice d'aretes utiles li�es au peptide               std::vector<std::vector<std::size_t>*>
	std::cout << "________peptide_spectra________" << std::endl;
	for (auto iter : peptides_spectra) {
		PrintVector(*iter);
	}
	//spectra_peptides : indic� par les spectres. B[spectre] = liste des indices d'aretes utiles li�es au spectre			  std::vector<std::vector<std::size_t>*>
	std::cout << "________spectra_peptides________" << std::endl;
	for (auto iter : spectra_peptides) {
		PrintVector(*iter);
	}
	//peptides_sure : associe chaque peptide � un nombre d'ar�tes d�j� s�lectionn� puisqu'elles sont seules					  std::vector<unsigned int>
	std::cout << "________peptides_sure________" << std::endl;
	PrintVector(peptides_sure);
	//useful_scores : liste des ar�tes utiles																				  std::vector<Score>
	std::cout << "________useful_scores________" << std::endl;
	PrintVector(useful_scores);
	//peptides_index : r�atribution des indices des peptides																  std::unordered_map<std::size_t, std::size_t>
	std::cout << "________peptides_index________" << std::endl;
	for (auto& couple : peptides_index) {
		std::cout << couple.first << " -> " << couple.second << std::endl;
	}
	//spectra_index : r�atribution des indices des spectres																	  std::unordered_map<std::size_t, std::size_t>
	std::cout << "________spectra_index________" << std::endl;
	for (auto& couple : spectra_index) {
		std::cout << couple.first << " -> " << couple.second << std::endl;
	}
	//proteins_index : r�atribution des indices des proteines																  std::unordered_map<std::size_t, std::size_t>
	std::cout << "________proteins_index________" << std::endl;
	for (auto& couple : proteins_index) {
		std::cout << couple.first << " -> " << couple.second << std::endl;
	}
	//useful_proteins : liste des prot�ines utiles (indice de base)															  std::vector<std::size_t>
	std::cout << "________useful_proteins________" << std::endl;
	PrintVector(useful_proteins);
	//peptides_proteins : indice sur les peptides. A1[peptide] = liste de couple (proteine ,prob)							  std::vector<std::vector<std::tuple<std::size_t, float>>*>
	std::cout << "________peptides_proteins________" << std::endl;
	for (auto iter : peptides_proteins) {
		for (auto& tup : *iter) {
			std::cout << "(" << std::get<0>(tup) << " ," << std::get<1>(tup) << ")";
		}
		std::cout << std::endl;
	} 
	//unuseful_peptides_proteins : indice sur les peptides inutiles. A2[peptide inutile] = liste de couple (protein ,prob)	  std::vector<std::vector<std::tuple<std::size_t, float>>*>
	std::cout << "________unuseful_peptides_proteins________" << std::endl;
	for (auto iter : unuseful_peptides_proteins) {
		for (auto& tup : *iter) {
			std::cout << "(" << std::get<0>(tup) << " ," << std::get<1>(tup) << ")";
		}
		std::cout << std::endl;
	}
	*/

#pragma region variables and model

	IloEnv env;
	IloModel model(env);

	IloNumVarArray Q = IloNumVarArray(env, n, 0, IloInfinity);
	IloNumVarArray Delta = IloNumVarArray(env, m1, 0, IloInfinity);
	IloBoolVarArray X = IloBoolVarArray(env, o);

#pragma endregion

/*
#pragma region variables and model2

	IloEnv env;
	IloModel model(env);

	IloNumVarArray Q = IloNumVarArray(env, n, 0, IloInfinity);
	IloNumVarArray Delta = IloNumVarArray(env, m1, 0, IloInfinity);
	IloBoolVarArray X = IloBoolVarArray(env, o - number_of_first_edge);

#pragma endregion
*/

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

/*
#pragma region constraints2

	for (std::size_t k = 0; k < l; ++k) {
		if (spectra_peptides[k]->size() >= 2) {
			IloExpr constraintX(env);
			for (std::size_t h = 1; h < spectra_peptides[k]->size(); ++h) {	
				constraintX += X[spectra_peptides[k]->at(h) - k - 1];
			}
			model.add(constraintX <= 1);
		}
	}

	for (std::size_t j = 0; j < m1; ++j) {
		IloExpr constraintDelta1(env);
		IloExpr constraintDelta2(env);
		for (std::tuple<std::size_t, float> edge : (*peptides_proteins[j])) {
			constraintDelta1 += Q[std::get<0>(edge)] * std::get<1>(edge);
			constraintDelta2 -= Q[std::get<0>(edge)] * std::get<1>(edge);
		}
		for (std::size_t h : (*peptides_spectra_first[j])) {
			constraintDelta1 -= 1;
			constraintDelta2 += 1;
			for (std::size_t h2 = 1; h2 < spectra_peptides[spectra_index[useful_scores[h]->spectrum]]->size(); ++h2) {
				constraintDelta1 += X[h + h2 - spectra_index[useful_scores[h]->spectrum] - 1];
				constraintDelta2 -= X[h + h2 - spectra_index[useful_scores[h]->spectrum] - 1];
			}
		}
		for (std::size_t h : (*peptides_spectra[j])) {
			constraintDelta1 -= X[h - spectra_index[useful_scores[h]->spectrum] - 1];
			constraintDelta2 += X[h - spectra_index[useful_scores[h]->spectrum] - 1];
		}
		constraintDelta1 -= Delta[j] + peptides_sure[j];
		constraintDelta2 -= Delta[j] - peptides_sure[j];
		model.add(constraintDelta1 <= 0);
		model.add(constraintDelta2 <= 0);
	}

#pragma endregion
*/

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

/*
#pragma region objective2

	IloExpr objective(env);

	for (std::size_t j = 0; j < m1; ++j) {
		objective += psi1 * Delta[j];
	}
	for (std::size_t j = 0; j < m2; ++j) {
		for (std::tuple<std::size_t, float> edge : (*unuseful_peptides_proteins[j])) {
			objective += psi1 * std::get<1>(edge) * Q[std::get<0>(edge)];
		}
	}
	for (std::size_t k = 0; k < l; ++k) {
		objective += psi2 * useful_scores[spectra_peptides[k]->front()]->score;
		for (std::size_t h = 1; h < spectra_peptides[k]->size(); ++h) {
			objective -= psi2 * useful_scores[spectra_peptides[k]->front()]->score * X[spectra_peptides[k]->at(h) - k - 1];
			objective += psi2 * useful_scores[spectra_peptides[k]->at(h)]->score * X[spectra_peptides[k]->at(h) - k - 1];
		}
	}
	model.add(IloMinimize(env, objective));

#pragma endregion
*/


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
			solution.Add_Protein(proteins_index[i], (float)valuesQ[i]);
		}
	}

	IloNumArray valuesX(env);
	cplex.getValues(X, valuesX);
	for (std::size_t h = 0; h < o; ++h) {
		if (valuesX[h] > 0.5) {
			solution.Add_Score(scores[h]);
		}
	}
#pragma endregion

/*
#pragma region generate solution2

	IloNumArray valuesQ(env);
	cplex.getValues(Q, valuesQ);
	for (std::size_t i = 0; i < n; ++i) {
		if (valuesQ[i] > 0) {
			solution.Add_Protein(proteins_index[i], (float)valuesQ[i]);
		}
	}

	IloNumArray valuesX(env);
	cplex.getValues(X, valuesX);
	bool found;
	for (std::size_t k = 0; k < l; ++k) {
		found = false;
		for (std::size_t h = 1; h < spectra_peptides[k]->size(); ++h) {
			if (valuesX[spectra_peptides[k]->at(h) - k - 1] > 0.5) {
				solution.Add_Score(useful_scores[spectra_peptides[k]->at(h)]);
				found = true;
			}
		}
		if (!found) {
			solution.Add_Score(useful_scores[spectra_peptides[k]->front()]);
		}
	}
#pragma endregion
*/


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