// utils_solver.h
#ifndef UTILS_SOLVER_H
#define UTILS_SOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <numeric> 
#include <cmath>
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

template <typename T>
std::vector<T> linspace_utils(T start, T end, std::size_t num) {
    std::vector<T> result(num);
    T step = (end - start) / (num - 1);

    for (std::size_t i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    
    return result;
}

template <typename T>
std::vector<T> vector_typecast(double* arr, int N) {
    std::vector<T> vec(N);
    for (int i=0; i<N; i++){
        vec[i] = static_cast<T>(arr[i]);
    }

    return vec;
}


// Add two vectors point wise
template <typename T>
std::vector<T> add_vectors(std::vector<T> arr1, std::vector<T> arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] + arr2[i];
    }

    return res;
}

// Subtract two vectors point wise (arr1 - arr2)
template <typename T>
std::vector<T> subtract_vectors(std::vector<T> arr1, std::vector<T> arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] - arr2[i];
    }

    return res;
}

// Multiply two vectors point wise
template <typename T>
std::vector<T> multiply_vectors(std::vector<T> arr1, std::vector<T> arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] * arr2[i];
    }

    return res;
}

// Multiply scalar to vector point wise
template <typename T>
std::vector<T> multiply_vector_scalar(std::vector<T> arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] * value1;
    }

    return res;
}

// Raise each element of the vector to the passed value
template <typename T>
std::vector<T> vector_float_power(std::vector<T> arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = std::pow(arr1[i], value1);
    }

    return res;
}


// Divide two vectors point wise
template <typename T>
std::vector<T> divide_vectors(std::vector<T> arr1, std::vector<T> arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] / arr2[i];
    }

    return res;
}

// Divide vector point wise by a scalar
template <typename T>
std::vector<T> divide_vector_scalar(std::vector<T> arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] / value1;
    }

    return res;
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

// For vectors, arrays from std library
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