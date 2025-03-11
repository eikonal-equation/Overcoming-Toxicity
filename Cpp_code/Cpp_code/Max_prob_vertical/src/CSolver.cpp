#include "CSolver.h"
#include "WriteToFile.h"

array<double, 2> PDE_Solver::f_drift(const double aX, const double aY, const int aCtrl)
{
	array<double, 2> drift;
	drift[0] = aX * (1 - aX) * ((1 - aY) * (fRks * (1 - fEpsilon * double(aCtrl)) - 1) + double(aCtrl) * fGamma * aX * aY);
	drift[1] = aY * (1 - aY) * (1 + aX * (fRks * (1 - fEpsilon * double(aCtrl)) - 1)) - double(aCtrl) * fGamma * pow(aY, 2) * aX * (1 - aX);
	return drift;
}


void PDE_Solver::InitializeMat() {
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 0; j < fNy + 1; j++) {// y loop;

			if (j == 0) {//zero population
				fMyValfn[i][j] = 0;
			}
			else if (i < fDefeat_indx) {
				fMyValfn[i][j] = 0; //Defeat zone
			}
			else if (i > fVic_indx) {
				fMyValfn[i][j] = 1; //Victory zone
			}
			else {
				fMyValfn[i][j] = 1e-6;
			}

			fMyPolicy[i][j] = 0;
		}
	}
}

//4D: i -> x ; j -> y ; a (0,1) -> ctrl ; m (0,1,2) -> (x,y,y_rho)
void PDE_Solver::precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat)
{
	for (int i = fDefeat_indx; i < fVic_indx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			double xi = fx0 + double(i) * fDx;
			double yj = fy0 + double(j) * fDy;

			array<double, 3> future_state;
			array<double, 2> drift;
			for (int a = 0; a < 2; a++) {

				future_state[0] = xi;
				future_state[1] = yj;

				drift = f_drift(xi, yj, a);

				for (int ii = 0; ii < 2; ii++) {
					future_state[ii] += fTau * drift[ii];
				}
				future_state[2] = fRho * future_state[1]; // y_rho = rho*y



				int x_indx = find_index(future_state[0], fDx, fx0, fxmax, fNx);
				int y_indx = find_index(future_state[1], fDy, fy0, fymax, fNy);
				int yrho_indx = find_index(future_state[2], fDy, fy0, fymax, fNy);

				array<int, 3> future_indx{ x_indx, y_indx, yrho_indx };

				//if (i == 400 && j == 600 && a == 1) {
				//	cout << future_state[0] << endl;
				//	cout << future_state[1] << endl;
				//	cout << future_state[2] << endl;
				//	cout << future_indx[0] << endl;
				//	cout << future_indx[1] << endl;
				//	cout << future_indx[2] << endl;
				//}

				for (int jj = 0; jj < 3; jj++) {
					fullstate_loc_mat[i][j][a][jj] = future_state[jj];
					fullstate_ind_mat[i][j][a][jj] = future_indx[jj];
				}


			}
		}
	}

}


array<double,4> PDE_Solver::ValueIte_single_point(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	array<double, 4> returnVal;
	int counter = 0;
	for (int a = 0; a < 2; a++) {
		//extract loc, indx, and tau
		int x_indx = afullstate_ind_mat[aI][aJ][a][0];
		int y_indx = afullstate_ind_mat[aI][aJ][a][1];
		int yrho_indx = afullstate_ind_mat[aI][aJ][a][2];
		double x_loc = afullstate_loc_mat[aI][aJ][a][0];
		double y_loc = afullstate_loc_mat[aI][aJ][a][1];
		double yrho_loc = afullstate_loc_mat[aI][aJ][a][2];

		array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn, x_indx, y_indx);
		double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

		array<double, 4> aVal_sq_rho = aGrid.Stencil_for_Bilinear_Interp(fMyValfn, x_indx, yrho_indx);
		double aVal_rho = aGrid.Bilinear_Interp(aVal_sq_rho, x_indx, yrho_indx, x_loc, yrho_loc);

		returnVal[counter] = aVal;
		returnVal[counter + 1] = aVal_rho;
		counter = counter + 2;
	}
	return returnVal;
}

double PDE_Solver::ValueIte_single_step(const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	double new_err = 0;
	array<double, 4> InterpVals;
	for (int i = fDefeat_indx; i < fVic_indx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop

			array<double, 2> w_possible;
			int counter = 0;

			InterpVals
				= ValueIte_single_point(i, j, afullstate_loc_mat, afullstate_ind_mat, aGrid);

			for (int ii = 0; ii < 2; ii++) {
				w_possible[ii] = fProb_not_arrival * InterpVals[counter] + fProb_arrival * InterpVals[counter+1];
				counter = counter + 2;
			}

			//find the current minimum among 4 combos of drugs
			double max_val = 0;
			bool mypolicy = 0;
			if (w_possible[0] >= w_possible[1]) {
				max_val = w_possible[0];
				mypolicy = 0;
			}
			else {
				max_val = w_possible[1];
				mypolicy = 1;
			}

			//update the value function and the optimal policy
			if (max_val > fMyValfn[i][j]) {
				double this_change = max_val - fMyValfn[i][j];

				fMyValfn[i][j] = max_val;

				new_err = max(new_err, this_change);

				fMyPolicy[i][j] = mypolicy;

			}
		}
	}
	return new_err;
}


void PDE_Solver::Main_VPI_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);
	CPolicyEvaluation myPEclass(fNx, fNy, fx0, fy0, fDx, fDy, fDefeat_thres, fVictory_thres, fDefeat_indx, fVic_indx, fProb_not_arrival, fProb_arrival);

	InitializeMat();

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	//std::abort();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	double old_err = -1;
	double new_err = -1;
	int PE_counter = 0;
	//Main GS Value Iterations on the full state model
	for (long iter = 0; iter < fMaxIter; ++iter)
		//for (long iter = 0; iter < 1; ++iter)
	{
		new_err = ValueIte_single_step(myfull_loc_mat, myfull_ind_mat, myGrid);
		if (iter == 0) { old_err = new_err; }
		if (new_err < 0.1 && new_err < fDelta * old_err) {
			old_err = new_err;
			/*myPEclass.PolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn);*/
			myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-6);
			PE_counter++;
		}

		cout << "Iteration #: " << iter + 1 << " with err = " << new_err << "\n";
		cout << "# of PE: " << PE_counter << "\n";
		if (new_err < fTol) { break; }
	}

	//myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-6);
	string aStr = to_string(int(fLamb*1000+1e-12));
	string aStr2 = to_string(int(fRho*1000 + 1e-12));
	string filename1 = "output/Max_prob_vertical_valuefn_rho" + aStr2 + "_lamb" + aStr + ".dat";
	string filename2 = "output/Max_prob_vertical_policy_rho" + aStr2 + "_lamb" + aStr + ".dat";
	io::writeToFile2D<double>(filename1, fMyValfn);
	io::writeToFile2D<bool>(filename2, fMyPolicy);

}


double PDE_Solver::Solver_with_output() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);
	CPolicyEvaluation myPEclass(fNx, fNy, fx0, fy0, fDx, fDy, fDefeat_thres, fVictory_thres, fDefeat_indx, fVic_indx, fProb_not_arrival, fProb_arrival);

	InitializeMat();

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	//std::abort();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	double old_err = -1;
	double new_err = -1;
	int PE_counter = 0;
	//Main GS Value Iterations on the full state model
	for (long iter = 0; iter < fMaxIter; ++iter)
		//for (long iter = 0; iter < 1; ++iter)
	{
		new_err = ValueIte_single_step(myfull_loc_mat, myfull_ind_mat, myGrid);
		if (iter == 0) { old_err = new_err; }
		if (new_err < 0.1 && new_err < fDelta * old_err) {
			old_err = new_err;
			/*myPEclass.PolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn);*/
			myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-6);
			PE_counter++;
		}

		cout << "Iteration #: " << iter + 1 << " with err = " << new_err << "\n";
		cout << "# of PE: " << PE_counter << "\n";
		if (new_err < fTol) { break; }
	}

	
	//myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-6);
	string aStr = to_string(int(fLamb+1e-12));
	string aStr2 = to_string(int(fRho*1000 + 1e-12));
	string filename1 = "output/Max_prob_vertical_valuefn_rho" + aStr2 + "_lamb" + aStr + ".dat";
	string filename2 = "output/Max_prob_vertical_policy_rho" + aStr2 + "_lamb" + aStr + ".dat";
	io::writeToFile2D<double>(filename1, fMyValfn);
	io::writeToFile2D<bool>(filename2, fMyPolicy);
	double xloc = 0.5;
	double yloc = 0.1;
	int xindex = int(xloc + 1e-12 / fDx);
	int yindex = int(yloc + 1e-12 / fDx);
	return fMyValfn[xindex][yindex];
}

void PDE_Solver::writeDomainToFile(const string aFilename) {

	array_1D params(boost::extents[11]); // Creating boostarray "row vector"

	params[0] = fNx;
	cout << "fNx = " << params[0] << endl;

	params[1] = fNy;
	cout << "fNy = " << params[1] << endl;

	params[2] = fDx;
	cout << "Delta x = " << params[2] << endl;

	params[3] = fDy;
	cout << "Delta y = " << params[3] << endl;

	params[4] = fEpsilon;
	cout << "Eps = " << params[4] << endl;

	params[5] = fRks;
	cout << "Rks = " << params[5] << endl;

	params[6] = fGamma;
	cout << "Gamma = " << params[6] << endl;

	params[7] = fLamb;
	cout << "Lamb = " << params[7] << endl;

	params[8] = fRho;
	cout << "Rho = " << params[8] << endl;

	params[9] = fVictory_thres;
	cout << "Victory thres = " << params[9] << endl;

	params[10] = fDefeat_thres;
	cout << "Defeat thres = " << params[10] << endl;


	// Writing vector of parameters to file
	io::writeToFile1D<double>(aFilename + "_DomainParameters.dat", params);
}
