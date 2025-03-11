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
			double xi = fx0 + double(i) * fDx;
			double yj = fy0 + double(j) * fDy;
			if (j == 0) {//zero population
				fMyValfn[i][j] = 0;
			}
			else if (xi * yj < fDefeat_thres) {
				fMyValfn[i][j] = 0; //Defeat zone
			}
			else if (xi * yj > fVictory_thres) {
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
	for (int i = 1; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			double xi = fx0 + double(i) * fDx;
			double yj = fy0 + double(j) * fDy;

			if (xi * yj >= fDefeat_thres && xi * yj <= fVictory_thres) {
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

					for (int jj = 0; jj < 3; jj++) {
						fullstate_loc_mat[i][j][a][jj] = future_state[jj];
						fullstate_ind_mat[i][j][a][jj] = future_indx[jj];
					}
				}

			}
		}
	}

}

void PDE_Solver::Read_policy() {
	if (fPolicy_name == "constitutive") {
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop

				fMyPolicy[i][j] = 1;
			}
		}
	}
	else if (fPolicy_name == "myopic") {
		// Open the binary file
		std::ifstream infile(fFile_name, std::ios::binary);
		if (!infile.is_open()) {
			std::cerr << "Failed to open file." << std::endl;
		}

		// Read the binary data into a buffer
		std::vector<char> buffer((fNx + 1) * (fNy + 1));
		infile.read(buffer.data(), buffer.size());

		// Close the file
		infile.close();

		// Populate the boost multi_array with the data from the buffer
		for (int i = 0; i < fNx + 1; ++i) {
			for (int j = 0; j < fNy + 1; ++j) {
				fMyPolicy[i][j] = static_cast<bool>(buffer[i * (fNx + 1) + j]);
			}
		}
	}
}

void PDE_Solver::Main_PE_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);
	CPolicyEvaluation myPEclass(fNx, fNy, fx0, fy0, fDx, fDy, fDefeat_thres, fVictory_thres, fDefeat_indx, fVic_indx, fProb_not_arrival, fProb_arrival);

	InitializeMat();

	Read_policy();

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	//std::abort();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";


	myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-9);

	string str1 = to_string(int(fLamb * 1000 + 1e-12));
	string str2 = to_string(int(fRho * 1000 + 1e-12));
	string filename1, filename2;
	if (fPolicy_name == "constitutive") {
		filename1 = "output/Perf_hyper_const_valuefn_rho" + str2 + "_lamb" + str1 + ".dat";
		filename2 = "output/Perf_hyper_const_policy_rho" + str2 + "_lamb" + str1 + ".dat";
	}
	else if (fPolicy_name == "myopic") {
		filename1 = "output/Perf_hyper_myo_valuefn_rho" + str2 + "_lamb" + str1 + ".dat";
		filename2 = "output/Perf_hyper_myo_policy_rho" + str2 + "_lamb" + str1 + ".dat";
	}

	io::writeToFile2D<double>(filename1, fMyValfn);
	io::writeToFile2D<bool>(filename2, fMyPolicy);

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
