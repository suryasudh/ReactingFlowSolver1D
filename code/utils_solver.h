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

// Structure to hold the parsed JSON data.
// You can modify this structure based on the keys you expect in your JSON file.
struct JsonData {
    std::unordered_map<std::string, std::string> stringValues;
    std::unordered_map<std::string, int> integerValues;
    std::unordered_map<std::string, double> floatValues;
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
std::vector<T> vector_vector_typecast(std::vector<double> arr1){
    std::vector<T> arr2(arr1.size());
    for (int i=0; i<arr1.size(); i++){
        arr2[i] = static_cast<T>(arr1[i]);
    }

    return arr2;
}


template <typename T>
std::vector<double> vector_vector_typecast_inv(std::vector<T> arr1){
    std::vector<double> arr2(arr1.size());
    for (int i=0; i<arr1.size(); i++){
        arr2[i] = static_cast<double>(arr1[i]);
    }

    return arr2;
}


template <typename T>
std::vector<T> vector_typecast(double* arr, int N) {
    std::vector<T> vec(N);
    for (int i=0; i<N; i++){
        vec[i] = static_cast<T>(arr[i]);
    }

    return vec;
}


template <typename T>
std::tuple<double*, int> vector_typecast_inv(std::vector<T> vec) {
    int N = vec.size();
    double* arr = new double[N];
    for (int i=0; i<N; i++){
        arr[i] = static_cast<double>(vec[i]);
    }

    return std::make_tuple(arr, N);
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


// Divide a scalar by vector point wise
template <typename T>
std::vector<T> divide_scalar_vector(T value1, std::vector<T> arr1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        if (arr1[i] != 0){
            res[i] =  value1 / arr1[i];
        } else {
            res[i] = static_cast<T>(999999999);
        }
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

template <typename T>
void print_vector(std::vector<T> vec) {
    std::cout << "[ ";
    for (int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << "] " << std::endl;
}

// Function to read and parse a JSON file using standard library functions.
// Takes the filename as input and populates the JsonData structure.
bool readAndParseJson(const std::string& filename, JsonData& data);

void file_writer(std::vector<double> time, std::vector<std::vector<double>> vec1, std::string fileName);

void file_writer2(std::vector<double> x, std::vector<std::vector<double>> vec1, std::string fileName);

void file_writer3(std::string filename_op, std::vector<double> t, std::vector<std::array<double, 2>> x);



// // --- Helper Functions ---

// // Template function for power, equivalent to np.float_power
// template<typename T>
// T float_power(T base, T exp) {
//     return std::pow(base, exp);
// }

// // Template function for element-wise maximum of an array and a scalar, equivalent to np.maximum(arr, value)
// template<typename T, size_t N>
// std::array<T, N> array_maximum(const std::array<T, N>& arr, T value) {
//     std::array<T, N> result;
//     for (size_t i = 0; i < N; ++i) {
//         result[i] = std::max(arr[i], value);
//     }
//     return result;
// }

// // Template function to sum elements of a std::array, equivalent to np.sum(arr)
// template<typename T, size_t N>
// T array_sum(const std::array<T, N>& arr) {
//     // Start summation with 0.0 of type T
//     return std::accumulate(arr.begin(), arr.end(), static_cast<T>(0.0));
// }



#endif // UTILS_SOLVER_H