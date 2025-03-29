// utils_solver.h
#ifndef UTILS_SOLVER_H
#define UTILS_SOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
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

#endif // UTILS_SOLVER_H