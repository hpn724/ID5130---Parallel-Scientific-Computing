Kindly use the make file to run the program. 

Run the program from the project directory itself and none of its subdirectories for a smooth execution. If the working directory is not set as project directory, the code may not find the data to import and plot.

make all will run the serial code and parallel code together.

make run_GCN_serial will run the serial code alone

make run_GCN_omp will run the parallel code for the specified number of threads in the src/omp/main.cpp function.
