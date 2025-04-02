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

template <int Precision>
struct PrecisionToType;

template <>
struct PrecisionToType<32> {
    using Type = float;
};

template <>
struct PrecisionToType<64> {
    using Type = double;
};

template <>
struct PrecisionToType<128> {
    using Type = long double;
};



/**
 * @brief Creates an array initialized with a specified value (type converted to T)
 * @tparam T Type of array elements
 * @tparam U Type of initialization value (automatically deduced)
 * @param size Number of elements in the array
 * @param value Initialization value (will be static_cast to T)
 * @return T* Newly allocated array
 */
template <typename T, typename U>
T* create_initialized_array(int size, U value) {
    if(size <= 0) {
        throw std::invalid_argument("Array size must be positive");
    }
    
    T* arr = new T[size];
    T cast_value = static_cast<T>(value);
    
    for(int i = 0; i < size; ++i) {
        arr[i] = cast_value;
    }
    
    return arr;
}


// 1. For raw pointers + size
template <typename T>
void print_array(const T* arr, size_t size) {
    std::cout << "[ ";
    bool first = true;
    for(size_t i = 0; i < size; ++i) {
        std::cout << (first ? "" : ", ") << arr[i];
        first = false;
    }
    std::cout << " ]\n";
}

template <typename Container>
void print_array(const Container& cont) {
    std::cout << "[ ";
    bool first = true;
    for(const auto& elem : cont) {
        std::cout << (first ? "" : ", ") << elem;
        first = false;
    }
    std::cout << " ]\n";
}


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