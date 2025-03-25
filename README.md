# dummyTheMummy
dummy The Mummy is a simple and innocent mummy who became dummy


Change of plans
---------------
Now running only on Ubuntu (linux based systems) <br>
Mainly because of the limitation of c++ cantera library <br>
which is more suited for linux environments.

--------

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


