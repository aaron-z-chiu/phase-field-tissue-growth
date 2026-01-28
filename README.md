# phase-field-tissue-growth
This repository contains the source code for the paper. It is used to generate the following results:

* **Convergence Test:** Section 4.1
* **2D Visualizations:**
    * Section 4.2.1 (Figure 2, left)
    * Section 4.2.2 (Figure 3, left)

## Authors
Zecheng Qiu, Yutong Wu, Junxiang Yang* \
\*Corresponding Author

## File Descriptions
* `PhaseFieldModel.h` & `PhaseFieldModel.cpp`: Core C++ source and header files implementing the semi-implicit finite difference solver.
* `main.cpp`: The main driver program.
* `lambda.m`: MATLAB script to visualize the effect of $\lambda_0$.
* `alpha.m`:   MATLAB script to visualize the effect of $\alpha$.
* `makefile`:  For macOS, Linux, WSL (Windows) systems, automates the compilation and cleaning process.

## Compilation and Execution
### MacOS / Linux / WSL (Windows)

#### Compiling and running the calculation codes
Open your terminal, and run the following commands:
1. `cd /path/to/your/codes`
2. `make all`
3. `./phase_field_model`

#### Plot Graphs & Check Results
* To visualize the growth rate effect ($\lambda_0 = 30, 100$), run `lambda.m` in MATLAB.
* To visualize the curvature effect ($\alpha = 30, 100$), run `alpha.m` in MATLAB.
* For convergence test results, see `temporal_convergence_results.txt` in the root folder.

#### Deleting all the results
Run `make delete` in your terminal.

**Warning**: Executing this command will **permanently remove** all generated data files (e.g., `.dat` files containing calculated results). This operation is irreversible. Please ensure you have backed up any critical results before proceeding.

