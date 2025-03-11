/*==============================================================================
 * File: main.cpp
 *
 * Author: MingYi Wang
 *
 * Description: This file initializes all the global variables and executes the
 * the corresponding example from the command line.
 *
 *============================================================================*/

 //-----------------------Project specific header files---------------------------
#include "CSolver.h"

//----------------------------Libraries------------------------------------------
#include <chrono>

using namespace std;

int main()
{
	setbuf(stdout, NULL);
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 16;
	const double gMax_population = 1;
	const double gMax_frequency = 1;
	const double gEpsilon = 0.2;
	const double gRks = 0.85;
	const double gGamma = 1.0;
	//const double gLamb = 1.0;
	//const double gRho = 0.65;
	const double gTol = 1e-6;
	const double gDefeat_thres = 0.005;
	const double gVictory_thres = 0.95;
	const string gPolicy_name = "myopic";
	array<double,10> lamb_list{0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0};
	//array<double,18> rho_list{0.5,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.60};
	array<double,2> rho_list{0.65, 0.7};
	
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();

	for (int i = 0; i < 10; i++){
		double gLamb = lamb_list[i];
		const string str1 = to_string(int(gLamb * 1000 + 1e-12));
		const string gFile_name = "../Max_freq_1stdilution/output/Sto_myopic_policy_lamb" + str1 + ".dat";
		for (int j = 0; j < 2; j++){
			double gRho = rho_list[j];
			cout << "Running Probabilistic Performance Metric Solver: " << endl;
			//Build the main class of our PDE solver
			PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol,
				gDefeat_thres, gVictory_thres, gEpsilon, gRks, gGamma, gLamb, gRho, gPolicy_name, gFile_name);

			//Calling the main PDE solver
			myExample.Main_PE_Solver();

			// Writing Grid Parameters to file
			std::string file_label = std::to_string(gRefinement_factor);
			myExample.writeDomainToFile("output/Perf_hyper");

			cout << "Successfully completed!" << endl;
		}
	}

	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
