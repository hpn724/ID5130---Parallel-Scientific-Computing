CXX		  := g++
CXX_FLAGS := -std=c++17 -fopenmp
PYTHON := python3

BIN		:= bin
SRC		:= src
INCLUDE	:= include
LIB		:= lib

LIBRARIES	:=
EXECUTABLE	:= main


all: data_import serial_GCN serial_GCN omp_GCN omp_plot

run_GCN_serial: data_import serial_GCN serial_plot

run_GCN_omp : data_import omp_GCN omp_plot

data_import : 
	$(PYTHON) src/data_import.py

serial_GCN : src/Serial/main.cpp src/Serial/GradDescentOptim.cpp src/Serial/GCNLayer.cpp src/Serial/SoftmaxLayer.cpp src/Serial/GCNNetwork.cpp

	$(CXX) $(CXX_FLAGS) -o src/Serial/main src/Serial/main.cpp src/Serial/operations.cpp src/Serial/GradDescentOptim.cpp src/Serial/GCNLayer.cpp src/Serial/SoftmaxLayer.cpp src/Serial/GCNNetwork.cpp 
	src/Serial/main

serial_plot :
	$(PYTHON) src/Serial/serial_plot.py

omp_GCN : src/omp/main.cpp src/omp/GradDescentOptim_omp.cpp src/omp/GCNLayer_omp.cpp src/omp/SoftmaxLayer_omp.cpp src/omp/GCNNetwork_omp.cpp

	$(CXX) $(CXX_FLAGS) -o src/omp/main src/omp/main.cpp src/omp/operations_omp.cpp src/omp/GradDescentOptim_omp.cpp src/omp/GCNLayer_omp.cpp src/omp/SoftmaxLayer_omp.cpp src/omp/GCNNetwork_omp.cpp
	src/omp/main

omp_plot :
	$(PYTHON) src/omp/omp_plot.py
