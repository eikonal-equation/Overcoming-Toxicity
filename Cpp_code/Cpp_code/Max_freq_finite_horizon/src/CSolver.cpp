#include "CSolver.h"
#include "WriteToFile.h"
#include<omp.h>

inline array<double, 2> PDE_Solver::f_drift(const double aX, const double aY, const int aCtrl)
{
	array<double, 2> drift;
	drift[0] = aX * (1 - aX) * ((1 - aY) * (fRks * (1 - fEpsilon * double(aCtrl)) - 1) + double(aCtrl) * fGamma * aX * aY);
	drift[1] = aY * (1 - aY) * (1 + aX * (fRks * (1 - fEpsilon * double(aCtrl)) - 1)) - double(aCtrl) * fGamma * pow(aY, 2) * aX * (1 - aX);
	return drift;
}

void PDE_Solver::ode_system(const state_type& state, state_type& dstate_dt, double t, int aCtrl) {
	array<double, 2> drift = f_drift(state[0], state[1], aCtrl);
	dstate_dt[0] = drift[0];    // dx/dt = fdot
	dstate_dt[1] = drift[1];   // dy/dt = Ndot
}

void PDE_Solver::InitializeMat(int aTime_slice) {
	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx;
		for (int j = 0; j < fNy + 1; j++) {// y loop

			if (aTime_slice == fM - 1) { // terminal condition
				if (j == 0) {//zero population
					fMyValfn_previous[i][j] = 0;
				}
				else {
					fMyValfn_previous[i][j] = xi;
				}

				if (j == 0 || i == 0) {//zero population or zero killers
					fMyValfn_current[i][j] = 0;
				}
				else if (i == fNx) {
					fMyValfn_current[i][j] = 1; // all killers
				}
				else {
					fMyValfn_current[i][j] = 1e-6;
				}

			}
			else {// lower time slices

				fMyValfn_previous[i][j] = fMyValfn_current[i][j];

				if (j == 0 || i == 0) {//zero population or zero killers
					fMyValfn_current[i][j] = 0;
				}
				else if (i == fNx) {
					fMyValfn_current[i][j] = 1; // all killers
				}
				else {
					fMyValfn_current[i][j] = 1e-6;
				}


			}
			fMyPolicy_current[i][j] = 0;
		}
	}
}



void PDE_Solver::InitializeMat_N(int aTime_slice) {
	double start_time = double(aTime_slice) * fTau;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx;
		for (int j = 0; j < fNy + 1; j++) {// y loop
			double yj = fy0 + double(j) * fDy;

			if (aTime_slice == fM - 1) { // terminal condition
				if (j == 0) {//zero population
					fMyValfn_previous_N[i][j] = 0;
				}
				else {
					fMyValfn_previous_N[i][j] = yj;
				}

				if (j == 0) {//zero population
					fMyValfn_current_N[i][j] = 0;
				}
				else if (i == 0 || i == fNx ) {//zero killer or zero sensitive; logistic growth
					std::vector<double> my_times;
					std::vector<double> my_f;
					std::vector<double> my_N;
					array<double, 2> xy{ xi,yj };
					// Lambda function as observer
					auto observer = [&my_f, &my_N, &my_times](const state_type& x, double t) {
						//my_states.push_back(x);
						my_f.push_back(x[0]);
						my_N.push_back(x[1]);
						my_times.push_back(t);
					};

					auto system_with_ctrl = [this](const state_type& state, state_type& dstate_dt, double t) {
						this->ode_system(state, dstate_dt, t, 0); };
					//integrate the ode
					integrate_adaptive(ctrl_rkck54(), system_with_ctrl, xy, start_time, fTf, fTau, observer);
					fMyValfn_current_N[i][j] = my_N.back();
				}
				else {
					fMyValfn_current_N[i][j] = 1e-6;
				}

			}
			else {// lower time slices

				fMyValfn_previous_N[i][j] = fMyValfn_current_N[i][j];

				if (j == 0) {//zero population or zero killers
					fMyValfn_current_N[i][j] = 0;
				}
				else if (i == 0 || i == fNx) {//zero killer or zero sensitive; logistic growth
					std::vector<double> my_times;
					std::vector<double> my_f;
					std::vector<double> my_N;
					array<double, 2> xy{ xi,yj };
					// Lambda function as observer
					auto observer = [&my_f, &my_N, &my_times](const state_type& x, double t) {
						//my_states.push_back(x);
						my_f.push_back(x[0]);
						my_N.push_back(x[1]);
						my_times.push_back(t);
					};

					auto system_with_ctrl = [this](const state_type& state, state_type& dstate_dt, double t) {
						this->ode_system(state, dstate_dt, t, 0); };
					//integrate the ode
					integrate_adaptive(ctrl_rkck54(), system_with_ctrl, xy, start_time, fTf, fTau, observer);
					fMyValfn_current_N[i][j] = my_N.back();
				}
				else {
					fMyValfn_current_N[i][j] = 1e-6;
				}


			}
		}
	}
}

//4D: i -> x ; j -> y ; a (0,1) -> ctrl ; m (0,1) -> (x,y)
void PDE_Solver::precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat)
{
	for (int i = 1; i < fNx; i++) { // x loop
		double xi = fx0 + double(i) * fDx;
		for (int j = 1; j < fNy + 1; j++) {// y loop
			double yj = fy0 + double(j) * fDy;

			array<double, 2> future_state;
			array<double, 2> drift;
			for (int a = 0; a < 2; a++) {

				future_state[0] = xi;
				future_state[1] = yj;

				drift = f_drift(xi, yj, a);

				for (int ii = 0; ii < 2; ii++) {
					future_state[ii] += fTau * drift[ii];
				}

				int x_indx = find_index(future_state[0], fDx, fx0, fxmax, fNx);
				int y_indx = find_index(future_state[1], fDy, fy0, fymax, fNy);


				array<int, 2> future_indx{ x_indx, y_indx };

				for (int jj = 0; jj < 2; jj++) {
					fullstate_loc_mat[i][j][a][jj] = future_state[jj];
					fullstate_ind_mat[i][j][a][jj] = future_indx[jj];
				}


			}
		}
	}

}


inline array<double, 2> PDE_Solver::TimeMarch_single_point(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	array<double, 2> returnVal;
	int counter = 0;
	for (int a = 0; a < 2; a++) {
		//extract loc, indx, and tau
		int x_indx = afullstate_ind_mat[aI][aJ][a][0];
		int y_indx = afullstate_ind_mat[aI][aJ][a][1];

		double x_loc = afullstate_loc_mat[aI][aJ][a][0];
		double y_loc = afullstate_loc_mat[aI][aJ][a][1];

		array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn_previous, x_indx, y_indx);
		double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

		returnVal[a] = aVal;
	}
	return returnVal;
}


inline double PDE_Solver::TimeMarch_single_point_N(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	bool policy = fMyPolicy_current[aI][aJ];


	int x_indx = afullstate_ind_mat[aI][aJ][int(policy)][0];
	int y_indx = afullstate_ind_mat[aI][aJ][int(policy)][1];

	double x_loc = afullstate_loc_mat[aI][aJ][int(policy)][0];
	double y_loc = afullstate_loc_mat[aI][aJ][int(policy)][1];

	array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn_previous_N, x_indx, y_indx);
	double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

	return aVal;
}

array<double,2> PDE_Solver::Ite_single_point(const int aI, const int aJ, CGrid_2D& aGrid)
{
	double N = fMyN_lim_old[aI][aJ];
	double f = fMyfreq_lim_old[aI][aJ];
	double newN = fRho * N;
	double aVal_f = -1;
	double aVal_N = -1;
	/*cout << "f = " << f << endl;
	cout << "N = " << N << endl;
	cout << "newN = " << newN << endl;*/
	if (aI == 0 || aI == fNx){
		if (aI == 0){
			aVal_f = 0;
		}
		else if (aI == fNx){
			aVal_f = 1;
		}
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);
		array<double, 2> aVal_interval{fMyN_lim_old[aI][N_indx - 1], fMyN_lim_old[aI][N_indx]};
		aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}
	else{
		int f_indx = find_index(f, fDx, fx0, fxmax, fNx);
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);

		array<double, 4> aVal_sq_f = aGrid.Stencil_for_Bilinear_Interp(fMyfreq_lim_old, f_indx, N_indx);
		aVal_f = aGrid.Bilinear_Interp(aVal_sq_f, f_indx, N_indx, f, newN);

		if (f_indx == fNx)
		{// project it back to f = 1 - dx for biological reasons
			array<double, 2> aVal_interval{fMyN_lim_old[fNx - 1][N_indx - 1], fMyN_lim_old[fNx - 1][N_indx]};
			aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		}
		else{
			array<double, 4> aVal_sq_N = aGrid.Stencil_for_Bilinear_Interp(fMyN_lim_old, f_indx, N_indx);
			aVal_N = aGrid.Bilinear_Interp(aVal_sq_N, f_indx, N_indx, f, newN);
		}
		

		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}	
}


array<double,2> PDE_Solver::Ite_single_point_one_time_update(const int aI, const int aJ, CGrid_2D& aGrid)
{
	double N = fMyN_lim_old[aI][aJ];
	double f = fMyfreq_lim_old[aI][aJ];
	double newN = fRho * N;

	double aVal_f = -1;
	double aVal_N = -1;

	if (aI == 0 || aI == fNx){
		if (aI == 0){
			aVal_f = 0;
		}
		else if (aI == fNx){
			aVal_f = 1;
		}
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);
		array<double, 2> aVal_interval{fMyN_map[aI][N_indx - 1], fMyN_map[aI][N_indx]};
		aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}
	else{
		int f_indx = find_index(f, fDx, fx0, fxmax, fNx);
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);

		array<double, 4> aVal_sq_f = aGrid.Stencil_for_Bilinear_Interp(fMyfreq_map, f_indx, N_indx);
		double aVal_f = aGrid.Bilinear_Interp(aVal_sq_f, f_indx, N_indx, f, newN);

		if (f_indx == fNx)
		{// project it back to f = 1 - dx for biological reasons
			array<double, 2> aVal_interval{fMyN_map[fNx - 1][N_indx - 1], fMyN_map[fNx - 1][N_indx]};
			aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		}
		else{
			array<double, 4> aVal_sq_N = aGrid.Stencil_for_Bilinear_Interp(fMyN_map, f_indx, N_indx);
			double aVal_N = aGrid.Bilinear_Interp(aVal_sq_N, f_indx, N_indx, f, newN);

		}
		
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}	
}


double PDE_Solver::Ite_single_step(CGrid_2D& aGrid)
{
	double new_err = 0;
	array<double,2> Interped_val;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			
			
			Interped_val = Ite_single_point(i, j,aGrid);

			fMyfreq_lim_new[i][j] = Interped_val[0];
			fMyN_lim_new[i][j] = Interped_val[1];

			double this_change = abs(fMyfreq_lim_new[i][j] - fMyfreq_lim_old[i][j]);

			new_err = max(new_err, this_change);
		}
	}
	return new_err;
}


double PDE_Solver::Ite_single_step_one_time_update(CGrid_2D& aGrid)
{
	double new_err = 0;
	array<double,2> Interped_val;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			
			
			Interped_val = Ite_single_point_one_time_update(i, j,aGrid);

			fMyfreq_lim_new[i][j] = Interped_val[0];
			fMyN_lim_new[i][j] = Interped_val[1];

			double this_change = abs(fMyfreq_lim_new[i][j] - fMyfreq_lim_old[i][j]);

			new_err = max(new_err, this_change);
		}
	}
	return new_err;
}



void PDE_Solver::Iterative_map(CGrid_2D& aGrid) {
	double new_err = -1;
	long MaxIter = 10000;

	//Main GS Value Iterations on the full state model
	for (long iter = 1; iter < MaxIter; ++iter)
	{
		new_err = Ite_single_step(aGrid);

		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyfreq_lim_new[i][j];
				fMyN_lim_old[i][j] = fMyN_lim_new[i][j];
			}
		}

		cout << "Iteration #: " << iter << " with err = " << new_err << "\n";\

		if (new_err < fTol) { break; }
	}
	
}



void PDE_Solver::Iterative_map_one_time_update(CGrid_2D& aGrid) {
	double new_err = -1;
	long MaxIter = 10000;

	string filename1 = "Max_freq_dilution_f.dat";
	string filename2 = "Max_freq_dilution_N.dat";
	io::writeToFile2D<double>(filename1, fMyfreq_lim_old);
	io::writeToFile2D<double>(filename2, fMyN_lim_old);

	//Main GS Value Iterations on the full state model
	for (long iter = 1; iter < MaxIter; ++iter)
	{
		new_err = Ite_single_step_one_time_update(aGrid);

		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyfreq_lim_new[i][j];
				fMyN_lim_old[i][j] = fMyN_lim_new[i][j];
			}
		}
		if (iter % 10 == 0){
			io::AppendToFile2D<double>(filename1, fMyfreq_lim_old);
			io::AppendToFile2D<double>(filename2, fMyN_lim_old);
		}
		
		cout << "Iteration #: " << iter << " with err = " << new_err << "\n";\
		fNum_dilution = iter;

		if (new_err < fTol) { break; }
	}
	if (fNum_dilution % 10  != 0){
		io::AppendToFile2D<double>(filename1, fMyfreq_lim_old);
		io::AppendToFile2D<double>(filename2, fMyN_lim_old);
	}
	
}


void PDE_Solver::Main_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);

	string filename1 = "Max_freq_finite_horizon_valuefn.dat";
	string filename2 = "Max_freq_finite_horizon_policy.dat";

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	//std::abort();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	//Main Time Marching Solver
	for (int t = fM - 1; t >= 0; t--)
	{

		InitializeMat(t);

		if (t == fM - 1) {
			io::writeToFile2D<double>(filename1, fMyValfn_previous);
			io::writeToFile2D<bool>(filename2, fMyPolicy_current);
		}

		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop

				//array<double, 2> w_possible;
				array<double, 2> InterpVals;
				InterpVals
					= TimeMarch_single_point(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);

				//for (int ii = 0; ii < 2; ii++) {
				//	w_possible[ii] =  InterpVals[ii];
				//}


				if (InterpVals[0] >= InterpVals[1]) {
					fMyValfn_current[i][j] = InterpVals[0];
					fMyPolicy_current[i][j] = 0;
				}
				else {
					fMyValfn_current[i][j] = InterpVals[1];
					fMyPolicy_current[i][j] = 1;
				}
				/*cout << " i = " << i << endl;
				cout << "j = " << j << endl;*/
			}

		}
		io::AppendToFile2D<double>(filename1, fMyValfn_current);
		io::AppendToFile2D<bool>(filename2, fMyPolicy_current);
		cout << "Time slice: " << t << endl;
	}
}



void PDE_Solver::Inf_limit_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);

	string filename1 = "Max_freq_finite_horizon_inf_limit_f.dat";
	string filename2 = "Max_freq_finite_horizon_inf_limit_N.dat";

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	//std::abort();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	//Main Time Marching Solver
	for (int t = fM - 1; t >= 0; t--)
	{

		InitializeMat(t);

		if (t == fM - 1) {
			io::writeToFile2D<double>(filename1, fMyValfn_previous);
			io::writeToFile2D<double>(filename2, fMyValfn_previous_N);
		}

		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop

				array<double, 2> InterpVals;
				InterpVals
					= TimeMarch_single_point(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);


				if (InterpVals[0] >= InterpVals[1]) {
					fMyValfn_current[i][j] = InterpVals[0];
					fMyPolicy_current[i][j] = 0;
				}
				else {
					fMyValfn_current[i][j] = InterpVals[1];
					fMyPolicy_current[i][j] = 1;
				}

			}

		}

		InitializeMat_N(t);
		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				double InterpVal
					= TimeMarch_single_point_N(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);
				fMyValfn_current_N[i][j] = InterpVal;
			}
		}


		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyValfn_current[i][j];
				fMyfreq_lim_new[i][j] = fMyfreq_lim_old[i][j];
				fMyN_lim_old[i][j] = fMyValfn_current_N[i][j];
				fMyN_lim_new[i][j] = fMyN_lim_old[i][j];

				fMyfreq_map[i][j] = fMyValfn_current[i][j];
				fMyN_map[i][j] = fMyValfn_current_N[i][j];
			}
		}


		Iterative_map(myGrid);

		io::AppendToFile2D<double>(filename1, fMyfreq_lim_new);
		io::AppendToFile2D<double>(filename2, fMyN_lim_new);

		cout << "Time slice: " << t << endl;
	}

	//Iterative_map_one_time_update(myGrid);
	// Iterative_map(myGrid);
	// string aStr = to_string(int((1.0/fTf*1000 + 1e-12)));
	// string aStr2 = to_string(int(fRho*1000 + 1e-12));
	// string filename1 = "output/limit_f_rho" + aStr2 + "_lamb" + aStr + ".dat";
	// string filename2 = "output/limit_N_rho" + aStr2 + "_lamb" + aStr + ".dat";
	// io::writeToFile2D<double>(filename1, fMyfreq_lim_new);
	// io::writeToFile2D<double>(filename2, fMyN_lim_new);
	// cout << aStr2 << endl;
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

	params[7] = fTf;
	cout << "Tf = " << params[7] << endl;

	params[8] = fM;
	cout << "fM = " << params[8] << endl;

	params[9] = fTau;
	cout << "Tau = " << params[9] << endl;

	params[10] = fNum_dilution;
	cout << "Num of dilution = " << params[10] << endl;


	// Writing vector of parameters to file
	io::writeToFile1D<double>(aFilename + "_DomainParameters.dat", params);
}
