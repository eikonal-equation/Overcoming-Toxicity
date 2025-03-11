#include "CPolicyEvaluation.h"
#include "WriteToFile.h"

void CPolicyEvaluation::PolicyEvaluation(const array_2D_bool& aPolicySlice,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat,
	array_2D& aValfn) {
	int mSize = (fNx + 1) * (fNy + 1);

	VectorXd b(mSize);
	VectorXd v(mSize);
	SparseMatrix<double> A(mSize, mSize);
	vector<Triplet<double>> entries;

	int aXindex;
	int aYindex;
	int aYrho_index;
	double aXloc;
	double aYloc;
	double aYrho_loc;
	double x_pos_tr;
	double y_pos_tr;
	double yrho_pos_tr;


	for (int i = 0; i < fNx + 1; i++) {
		double xi = fx0 + double(i) * fDx;
		for (int j = 0; j < fNy + 1; j++) {

			entries.push_back(Triplet<double>(i * (fNx + 1) + j, i * (fNx + 1) + j, 1.0));

			//A_dense[i * (fN + 1) + j][i * (fN + 1) + j] = 1;

			if (j > 0 && i > 0 && i < fNx && aValfn[i][j] != 1 && aValfn[i][j] != 0)
			{
				int current_policy = int(aPolicySlice[i][j]);

				b(i * (fNx + 1) + j) = fProb_arrival*xi;


				aXindex = afullstate_ind_mat[i][j][current_policy][0];
				aYindex = afullstate_ind_mat[i][j][current_policy] [1];

				aXloc = afullstate_loc_mat[i][j][current_policy][0];
				aYloc = afullstate_loc_mat[i][j][current_policy][1];


				x_pos_tr = aXindex * fDx + fx0;
				y_pos_tr = aYindex * fDy + fy0;

				// without rho
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (x_pos_tr - aXloc) / fDx)); //bottom left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //bottom right
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (x_pos_tr - aXloc) / fDx)); //top left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //top right


			}
			else {
				b(i * (fNx + 1) + j) = aValfn[i][j];
			}
		}
	}

	A.setFromTriplets(entries.begin(), entries.end());

	//string afile = "A dense linear.dat";
	//io::writeToFile2D(afile, A_dense);

	//-------------- solving the linear system Av = b by LU factorization ----------------------
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;

	solver.compute(A);

	if (solver.info() != Success) {
		cout << "LU Decomposition Failed" << endl;
		myAssert(solver.info() == Success);
	}
	v = solver.solve(b);

	//------------ update the top slice of our value function-----------------
	for (int i = 0; i < fNx + 1; i++) {
		for (int j = 0; j < fNy + 1; j++) {
			aValfn[i][j] = v(i * (fNx + 1) + j);
		}
	}
}

void CPolicyEvaluation::ApproximatePolicyEvaluation(const array_2D_bool& aPolicySlice,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat,
	array_2D& aValfn, double tolerance) {

	int mSize = (fNx+ 1) * (fNy + 1);

	VectorXd b(mSize);
	VectorXd v(mSize);
	VectorXd InitialGuess(mSize);
	SparseMatrix<double> A(mSize, mSize);
	vector<Triplet<double>> entries;

	int aXindex;
	int aYindex;
	int aYrho_index;
	double aXloc;
	double aYloc;
	double aYrho_loc;
	double x_pos_tr;
	double y_pos_tr;
	double yrho_pos_tr;

	for (int i = 0; i < fNx + 1; i++) {
		double xi = fx0 + double(i) * fDx;
		for (int j = 0; j < fNy + 1; j++) {

			entries.push_back(Triplet<double>(i * (fNx + 1) + j, i * (fNx + 1) + j, 1.0));
			InitialGuess(i * (fNx + 1) + j) = aValfn[i][j];
			//A_dense[i * (fN + 1) + j][i * (fN + 1) + j] = 1;

			if (j > 0 && i > 0 && i < fNx && aValfn[i][j] != 1 && aValfn[i][j] != 0)
			{
				int current_policy = int(aPolicySlice[i][j]);

				b(i * (fNx + 1) + j) = fProb_arrival*xi;


				aXindex = afullstate_ind_mat[i][j][current_policy][0];
				aYindex = afullstate_ind_mat[i][j][current_policy][1];

				aXloc = afullstate_loc_mat[i][j][current_policy][0];
				aYloc = afullstate_loc_mat[i][j][current_policy][1];


				x_pos_tr = aXindex * fDx + fx0;
				y_pos_tr = aYindex * fDy + fy0;

				// without rho
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex - 1, -fProb_not_arrival*(y_pos_tr - aYloc) / fDy * (x_pos_tr - aXloc) / fDx)); //bottom left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //bottom right
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (x_pos_tr - aXloc) / fDx)); //top left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //top right


			}
			else {
				b(i * (fNx + 1) + j) = aValfn[i][j];
			}
		}
	}

	A.setFromTriplets(entries.begin(), entries.end());

	BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double> > approxSolver;
	//BiCGSTAB<SparseMatrix<double>, DiagonalPreconditioner<double> > approxSolver;


	approxSolver.compute(A);
	approxSolver.setTolerance(tolerance);
	v = approxSolver.solveWithGuess(b, InitialGuess);
	//v = approxSolver.solve(b);
	cout << "Approx Solver Iterations:" << approxSolver.iterations() << "\n";

	//------------ update the top slice of our value function-----------------
	for (int i = 0; i < fNx + 1; i++) {
		for (int j = 0; j < fNy + 1; j++) {
			aValfn[i][j] = v(i * (fNx + 1) + j);
		}
	}
}
