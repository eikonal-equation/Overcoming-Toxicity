# Overcoming Toxicity

This repository contains the software source code and experimental data analysis used in the research paper "Overcomingtoxicity: how non-antagonistic
 microbes manage to thrive in boom-and-bust environments." The code and data here are intended for reproducible experiments, analysis, and visualization described in our manuscript.

# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

## Abstract of our manuscript
Antagonistic interactions are critical determinants of microbial community stability and composition, offering host benefits such as pathogen protection and providing avenues for antimicrobial control. While the ability to eliminate competitors confers an advantage to antagonistic microbes, it often incurs a fitness cost. Consequently, many microbes only produce toxins or engage in antagonistic behavior in response to specific cues like quorum sensing molecules or environmental stress. In laboratory settings, antagonistic microbes typically dominate over sensitive ones, raising the question of why both antagonistic and non-antagonistic microbes are found in natural environments and host microbiomes. Here, using both theoretical models and experiments with killer strains of _Saccharomyces cerevisiae_, we show that "boom-and-bust" dynamics -- periods of rapid growth punctuated by episodic mortality events -- caused by temporal environmental fluctuations can favor non-antagonistic microbes that do not incur the growth rate cost of toxin production. Additionally, using control theory, we derive bounds on the competitive performance and identify optimal regulatory toxin-production strategies in various boom-and-bust environments where population dilutions occur either deterministically or stochastically over time. Our mathematical investigation reveals that optimal toxin regulation is much more beneficial to killers in stochastic, rather than deterministic, boom-and-bust environments. Overall, our findings offer a new perspective on how both antagonistic and non-antagonistic microbes can thrive under varying environmental conditions.


# Manuscript
The "Overcomingtoxicity: how non-antagonistic microbes manage to thrive in boom-and-bust environments" manuscript can be found [here](#)  

# Contributions & Acknowledgements
  * AG suggested the topic and the modeling approach, conducted the wetlab experiments, and contextualized.
  * AV guided the optimization framework and contributed several algorithmic and modeling ideas. 
  * MW contributed additional algorithmic ideas, developed all numerical codes, conducted in silic oexperiments, and wrote most of the SI Appendix.
  * All three authors have analyzed the experimental results and equally contributed to writing the main manuscript.
  * We thank Andrew W. Murray, Stephen P. Ellner, Tobias Dorr, the editor and two anonymous reviewers for insightful comments on the manuscript, and Elisa Garabello for valuable discussions. AG thanks Marie F. Gorwa-Grauslund for gifting CyOFP1opt. AG acknowledges support by the National Institute of General Medical Sciences of the National Institutes of Health, United States of America under award number 1R35GM147493.
MW and AV acknowledge support by the National Science Foundation (awards DMS-1645643 and DMS-2111522) as well as the Air Force Office of Scientific Research (award FA9550-22-1-0528).

# Instructions

## Requirements
* The codebase primarily uses C++, with additional components in Jupyter Notebooks (Python), Mathematica, and MATLAB.
* Ensure you have the following installed:
    * C++ compiler.
    * MATLAB and Mathematica for their respective scripts.
    * Python (>=3.9) with the necessary packages listed in requirements.txt or in each notebook.
    * Jupyter Notebook or JupyterLab.
      
* The C++ code requires the users to install the "[Boost](https://www.boost.org/)" library (external) and "[Eigen](https://eigen.tuxfamily.org/)" library (external). 
    * We install it directly under the directory `Software_Code_for_Numerical_Experiments/Solver_Cpp_code/lib`. Please feel free to install it anywhere you want on your end but you would have to modify its path in the ***Makefile*** described below.

## Running the C++ Code: ##
The following instructions explain how to run the Solver for several PDEs constructed in the manuscript using the ***Makefile***. 

To change compilers, edit `CC=` by putting your desired compiler after the `=` sign. The default compiler is set to `g++`. 

To update the path/version of the "Boost" library you've installed, edit `BOOSTPATH` by putting your own path/version after the `=` sign. The current path/version is `./boost_1_79_0`.

To compile the code, type `make main` in the terminal in this folder. 
* Please note that you need to enter the `Risk-Neutral` folder to run the PDE solver for the risk-neutral value function/policy, and the `Risk-Aware` folder to run the PDE solver for the risk-aware value function/policy, respectively.
* In either case, after compilation, type `./main` to run the respective PDE solver.
* 
To delete the executable and the data files, type `make clean`.

The datafiles generated by running these solvers are stored in `Software_Code_for_Numerical_Experiments/Solver_Cpp_code/ouput` as default.

## Running MATLAB/Mathematica Code
- For MATLAB:  
  Open the `demo.m` scripts in MATLAB and run them as described in code comments. (Feel free to change the global variable values if you have modified the code yourself.)
- For Mathematica:  
  Open the `.nb` files in Mathematica.

## Running the Jupyter Notebooks
1. Clone this repository:
    ```bash
    git clone https://github.com/eikonal-equation/Overcoming-Toxicity.git
    cd Overcoming-Toxicity
    ```
2. Install required Python packages:
    ```bash
    pip install -r requirements.txt
    ```
3. Launch Jupyter Notebook or JupyterLab:
    ```bash
    jupyter notebook
    ```
4. Open and execute the desired notebook(s) in your browser.
