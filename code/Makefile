# Declare phony targets (not files)
.PHONY: all clean build

# Default target
all: clean build

# Clean build artifacts
clean:
	rm -f main01.cxx

# Build the executable
build: main01.cxx

# Explicit build rule with dependencies
main01.cxx: main01.cpp
	g++ utils_solver.cpp log_heap.cpp derivatives.cpp solver_funcs.cpp main01.cpp -o main01.cxx -O3 $$(pkg-config --cflags --libs cantera)
