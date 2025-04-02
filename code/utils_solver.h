// utils_solver.h
#ifndef UTILS_SOLVER_H
#define UTILS_SOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <iomanip>
#include <unordered_map>

// Structure to hold the parsed JSON data.
// You can modify this structure based on the keys you expect in your JSON file.
struct JsonData {
    std::unordered_map<std::string, std::string> stringValues;
    std::unordered_map<std::string, int> integerValues;
    std::unordered_map<std::string, double> floatValues;
};

// Function to read and parse a JSON file using standard library functions.
// Takes the filename as input and populates the JsonData structure.
bool readAndParseJson(const std::string& filename, JsonData& data);


void file_writer(std::vector<double> time, std::vector<std::vector<double>> vec1, std::string fileName);

void file_writer2(std::vector<double> x, std::vector<std::vector<double>> vec1, std::string fileName);

void file_writer3(std::string filename_op, std::vector<double> t, std::vector<std::array<double, 2>> x);


#endif // UTILS_SOLVER_H