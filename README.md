# dummyTheMummy
dummy The Mummy is a simple and innocent mummy who became dummy


Change of plans
---------------
Now running only on Ubuntu (linux based systems) <br>
Mainly because of the limitation of c++ cantera library <br>
which is more suited for linux environments.

Ubuntu based runs
-----------------
Command to run the solver with memory analysis <br>
<code> /usr/bin/time -v ./main01.cxx 2>&1 | tee -a "/media/ssudhakar/DATA_10TB/data_solver/logs/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver without memory analysis <br>
<code> ./main01.cxx | tee -a "../logfiles/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver with memory analysis <br>
<code> valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main01.cxx 2>&1 | tee -a "../logfiles/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

Command to run the solver with memory analysis <br>
<code> valgrind --tool=massif --massif-out-file=/media/ssudhakar/DATA_10TB/data_solver/valgrindLog/massif_output_$(date +\%Y_\%m_\%d_\%H_\%M_\%S).%p ./main01.cxx 2>&1 | tee -a "../logfiles/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

<code> sudo sysctl -w kernel.perf_event_paranoid=1 </code>
<code> perf stat ./main01.cxx 2>&1 | tee -a "../logfiles/log_output_$(date '+%Y_%m_%d_%H_%M_%S').txt" </code>

--------

Windows based runs
------------------
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




