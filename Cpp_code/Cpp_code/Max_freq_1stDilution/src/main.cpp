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
	const double gTol = 1e-6;
	array<double,10> lamb_list{0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0};
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();

	for (int i = 0; i < 10; i++){
		double gLamb = lamb_list[i];
		cout << "Lamb = " << gLamb << endl;
		cout << "Running Max Freq 1st Dilution Solver: " << endl;
		//Build the main class of our PDE solver
		PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gLamb);

		//Calling the main PDE solver
		myExample.Main_VPI_Solver();


		// Writing Grid Parameters to file
		std::string file_label = std::to_string(gRefinement_factor);
		myExample.writeDomainToFile("output/Sto_myopic");
		cout << "Successfully completed!" << endl;
	}
	
	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	//=============================================END of Our PDE Solver=====================================================
	
	return 0;
}
