
#include "Model.h"
#include <unordered_map>
#include <set>

/*TODO
faire la documentation
revoir l'organisation
*/
//__________________________________________________________________________________________________________

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

int main() {
	Model model;

	model.Load_Proteins("uniprot_humain_moitie.fasta");
	//model.Load_Proteins("prot.txt");
	std::cout << "proteins loaded : " << model.Number_Of_Proteins() << std::endl;
	model.In_Silico_Digestion(0 ,0);
	std::cout << "proteins digested : " << model.Number_Of_Peptides() << std::endl;
	model.Build_Theoretical_Spectra();
	std::cout << "theoretical spectra built : " << model.Number_Of_Peptides() << std::endl;
	model.Define_Probabilities(0.2f);
	std::cout << "probabilities defined" << std::endl;

	/*
	for (std::size_t i = 0; i < model.Number_Of_Peptides(); i++) {
		std::cout << model.Get_Peptide(i).Get_Sequence() << std::endl;
		PrintVector2<Pic*>(*(model.Get_Peptide(i).Build_Spectrum()));
	}
	*/


	/*
	std::unordered_map<size_t, unsigned int> sample;
	size_t j = 0;
	unsigned int compteur = 0;
	while (j < model.Number_Of_Proteins() && compteur < 20) {
		compteur += (unsigned int)model.Get_Protein(j).Get_Peptides().size();
		sample[j] = 1;
		j++;
	}
	model.Simulated_Sample({ sample });
	*/

	//model.Load_Spectra("save.ms2");
	//model.Load_Spectra("110616_yeast_ups_10fmol.ms2");
	std::unordered_map<std::size_t, unsigned int> sample = { {2 , 5} ,{3 , 2}	 };
	model.Simulated_Sample(sample);

	std::cout << "spectra loaded : " << model.Number_Of_Spectra() << std::endl;
	//model.Compute_Score(compute_score_function);
	model.Compute_Score(1);
	//model.Compute_Score2(0 ,0 ,2 ,0 ,2 ,2);
	std::cout << "scores computed : " << model.Number_Of_Scores() << std::endl;

	/*
	std::cout << "____________________" << std::endl;
	for (std::size_t i = 0; i < model.Number_Of_Scores(); ++i) {
		std::cout << model.Get_Score(i) << std::endl;
	}
	*/

	model.Solve(0.5f, 0.5f);
	std::cout << "model solved" << std::endl;

	model.Print_Solution();

	return 0;

}

/*
#include <iostream>
#include "ilcplex/ilocplex.h"
#include <math.h>

typedef IloArray<IloNumVarArray> NumVar2D;
typedef IloArray<NumVar2D> NumVar3D;

int main() {

#pragma region Define datas
	int nS = 4;
	int nD = 3;

	int* S = new int[nS] {10, 30, 40, 20};
	int* D = new int[nD] {20, 50, 30};

	int** C = new int* [nS];

	C[0] = new int[nD] {2, 3, 4};
	C[1] = new int[nD] {3, 2, 1};
	C[2] = new int[nD] {1, 4, 3};
	C[3] = new int[nD] {4, 5, 2};
#pragma endregion


#pragma region Define model
	IloEnv env;
	IloModel model(env);
#pragma endregion


#pragma region Define decision variables
	NumVar2D X{ env ,nS };
	for (int s = 0; s < nS; s++) {
		X[s] = IloNumVarArray(env, nD, 0, IloInfinity, ILOINT);
	}
#pragma endregion


#pragma region Define objective function
	IloExpr exp0(env);

	for (int s = 0; s < nS; s++) {
		for (int d = 0; d < nD; d++) {
			exp0 += C[s][d] * X[s][d];
		}
	}

	model.add(IloMinimize(env, exp0));
#pragma endregion


#pragma region Define constraints
	for (int s = 0 ;s < nS ;s++){
		IloExpr exp1(env);
		for (int d = 0; d < nD; d++) {
			exp1 += X[s][d];
		}
		model.add(exp1 <= S[s]);
	}

	for (int d = 0; d < nD; d++) {
		IloExpr exp2(env);
		for (int s = 0; s < nS; s++) {
			exp2 += X[s][d];
		}
		model.add(exp2 >= D[d]);
	}
#pragma endregion


#pragma region Solve
	IloCplex cplex(model);
	cplex.solve();
#pragma endregion
}
*/

/*
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
int main(int argc, char** argv){
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0.0, 40.0));
		vars.add(IloNumVar(env));
		vars.add(IloNumVar(env));
		model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
		model.add(-vars[0] + vars[1] + vars[2] <= 20);
		model.add(vars[0] - 3 * vars[1] + vars[2] <= 30);
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << endl;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	env.end();
	return 0;
}
*/

/*
#include <ilcplex/ilocplex.h>

int main() {
	IloEnv env;
	IloModel model(env);

	IloNumVar x1 =IloNumVar(env ,0 ,IloInfinity);
	IloNumVar x2 = IloNumVar(env, 0, IloInfinity);

	IloRange sum_to_one = IloRange(env, 1, 1);

	IloObjective obj = IloMaximize(env, 0);

	sum_to_one.setLinearCoef(x1, 1);
	sum_to_one.setLinearCoef(x2, 1);
	obj.setLinearCoef(x1, 2);
	obj.setLinearCoef(x2, 1);

	model.add(obj);
	model.add(sum_to_one);

	IloCplex solver(model);

	solver.solve();

	std::cout << solver.getObjValue() << std::endl;

	return 0;
}
*/