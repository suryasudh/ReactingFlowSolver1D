# ReactingFlowSolver1D
This solver is developed as a part of the course work <br>
of PH354: Computational Physics <br>
taken in the January-April Semester 2025. <br><br>

Author: Surya Datta Sudhakar <br>
SR No: 24309 <br>
Dept: Computational and Data Sciences <br><br>

The following code package needs <br>
cantera c++ library. <br>
<a href=https://cantera.org/stable/install/index.html#installing-the-cantera-c-interface-fortran-90-module>Cantera Installation</href> <br> <br>

This package was tested on 3 different cases. <br>
Contains 3 cases: <br>
1. Standing Wave (Non-Reacting) (enter case 1)
2. Homogeneous Auto Ignition (Reacting - Uniform) (enter case 3)
3. Gaussian Distribution based initial temperature condition (Reacting - Non Uniform) (enter case 2)

<br>
The input conditions need to be passed from `configuration.json` in the code folder. <br>
From the same folder, use `make` commands (examples provided below)

Change of plans
---------------
Now running only on Ubuntu (linux based systems) <br>
Mainly because of the limitation of c++ cantera library <br>
which is more suited for linux environments.

Ubuntu based runs
-----------------
Command to run the solver with memory analysis <br>
<code> /usr/bin/time -v ./main01.cxx 2>&1 | tee -a "/media/ssudhakar/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver without memory analysis <br>
<code> ./main01.cxx | tee -a "/media/ssudhakar/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver with memory analysis <br>
<code> valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main01.cxx 2>&1 | tee -a "/media/ssudhakar/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver with memory analysis <br>
<code> valgrind --tool=massif --massif-out-file=/media/ssudhakar/data_solver/valgrindLog/massif_output_$(date +\%Y_\%m_\%d_\%H_\%M_\%S).%p ./main01.cxx 2>&1 | tee -a "/media/ssudhakar/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

In case your PC is not set to perform profiling earlier, this might help.
<code> sudo sysctl -w kernel.perf_event_paranoid=1 </code>
<code> perf stat ./main01.cxx 2>&1 | tee -a "/media/ssudhakar/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

--------

Windows based runs
------------------
(deprecated)
First make were done on windows <br>
Keep the code files and CMakeLists.txt file in the code folder. <br>
Create a build directory in the code folder. <br>
Move into the folder. <br>
Run the following command (on Windows) <br>
<code> >>> cmake .. -G "MinGW Makefiles" </code>

Now, once the Makefile is created, <br>
(This is also the command to run after any changes to the files) <br>
run the next command (on Windows) <br>
<code> >>> mingw32-make </code>

To clean the code files generated, <br>
run the following command (on Windows) <br>
<code> >>> mingw32-make clean </code>

Once the make command is run without any errors <br>
The executables would be created. <br>
To run the program, execute the following command (on Windows) <br>
<code> >>> .\main001.exe </code>




