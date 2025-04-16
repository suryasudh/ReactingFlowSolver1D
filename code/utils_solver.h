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
    T step = (end - start) / (num);

    for (std::size_t i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    
    return result;
}

template <typename T>
std::vector<T> vector_vector_typecast(const std::vector<double>& arr1){
    std::vector<T> arr2(arr1.size());
    for (int i=0; i<arr1.size(); i++){
        arr2[i] = static_cast<T>(arr1[i]);
    }

    return arr2;
}


template <typename T>
std::vector<double> vector_vector_typecast_inv(const std::vector<T>& arr1){
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
std::tuple<double*, int> vector_typecast_inv(const std::vector<T>& vec) {
    int N = vec.size();
    double* arr = new double[N];
    for (int i=0; i<N; i++){
        arr[i] = static_cast<double>(vec[i]);
    }

    return std::make_tuple(arr, N);
}



// Add two vectors point wise
template <typename T>
std::vector<T> add_vectors(const std::vector<T>& arr1, const std::vector<T>& arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] + arr2[i];
    }

    return res;
}

template <typename T>
std::vector<std::vector<T>> add_vectors(const std::vector<std::vector<T>>& arr1, const std::vector<std::vector<T>>& arr2){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j = 0; j < arr1[0].size(); j++){
            res[i][j] = arr1[i][j] + arr2[i][j];
        }
    }

    return res;
}



// Subtract two vectors point wise (arr1 - arr2)
template <typename T>
std::vector<T> subtract_vectors(const std::vector<T>& arr1, const std::vector<T>& arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] - arr2[i];
    }

    return res;
}


// Subtract two vectors point wise (arr1 - arr2)
// @overload 
template <typename T>
std::vector<std::vector<T>> subtract_vectors(const std::vector<std::vector<T>>& arr1, const std::vector<std::vector<T>>& arr2){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j = 0; j < arr1[0].size(); j++){
            res[i][j] = arr1[i][j] - arr2[i][j];
        }
    }

    return res;
}


// Multiply two vectors point wise
template <typename T>
std::vector<T> multiply_vectors(const std::vector<T>& arr1, const std::vector<T>& arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] * arr2[i];
    }

    return res;
}

// Multiply scalar to vector point wise
template <typename T>
std::vector<T> multiply_vector_scalar(const std::vector<T>& arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] * value1;
    }

    return res;
}


// Multiply two vectors point wise
template <typename T>
std::vector<std::vector<T>> multiply_vectors_2nd_order(const std::vector<std::vector<T>>& arr1, const std::vector<std::vector<T>>& arr2){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j = 0; j < arr1[0].size(); j++){
            res[i][j] = arr1[i][j] * arr2[i][j];
        }
    }

    return res;
}

// Multiply vector-scalar point wise
template <typename T>
std::vector<std::vector<T>> multiply_vector_2nd_order_scalar(const std::vector<std::vector<T>>& arr1, T value1){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j = 0; j < arr1[0].size(); j++){
            res[i][j] = arr1[i][j] * value1;
        }
    }

    return res;
}


// Reduce vector of 2nd order to vector of 1st order by adding across the 1st axis
template <typename T>
std::vector<T> vector_reduction_by_sum1(const std::vector<std::vector<T>>& arr){
    int grid_size = arr.size();
    int n_species = arr[0].size();
    
    std::vector<T> res(grid_size);
    
    for (int j = 0; j < grid_size; j++){
        for (int i = 0; i < n_species; i++){
            res[j] += arr[j][i];
        }
    }

    return res;
}

// Raise each element of the vector to the passed value
template <typename T>
std::vector<T> vector_float_power(const std::vector<T>& arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = std::pow(arr1[i], value1);
    }

    return res;
}


// Divide two vectors point wise
template <typename T>
std::vector<T> divide_vectors(const std::vector<T>& arr1, const std::vector<T>& arr2){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] / arr2[i];
    }

    return res;
}

// Divide vector point wise by a scalar
template <typename T>
std::vector<T> divide_vector_scalar(const std::vector<T>& arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        res[i] = arr1[i] / value1;
    }

    return res;
}


// Divide a scalar by vector point wise
template <typename T>
std::vector<T> divide_scalar_vector(T value1, const std::vector<T>& arr1){
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


// Divide two vectors point wise
template <typename T>
std::vector<std::vector<T>> divide_vectors_2nd_order(const std::vector<std::vector<T>>& arr1, const std::vector<std::vector<T>>& arr2){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j = 0; j < arr1[0].size(); j++){
            res[i][j] = arr1[i][j] / arr2[i][j];
        }
    }

    return res;
}


// Divide two vectors point wise
template <typename T>
std::vector<std::vector<T>> divide_vector_2nd_order_1st_order(const std::vector<std::vector<T>>& arr1, const std::vector<T>& arr2){
    int grid_size = arr1.size();
    int n_species = arr1[0].size();
    
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < grid_size; i++){
        for (int j = 0; j < n_species; j++){
            res[i][j] = arr1[i][j] / arr2[i];
        }
    }

    return res;
}

// returns the selected column number
template <typename T>
std::vector<T> column_returner(const std::vector<std::vector<T>>& arr, int j){
    std::vector<T> res(arr.size());

    for (int i = 0; i < arr.size(); i++){
        res[i] = arr[i][j];
        }

    return res;
}


template <typename T>
std::vector<std::vector<T>> tiling_vector(const std::vector<T>& arr, int n_species){
    int grid_size = arr.size();
    std::vector<std::vector<T>> res(grid_size, std::vector<T>(n_species));
    
    for (int j = 0; j < grid_size; j++){
        for (int i = 0; i < n_species; i++){
            res[j][i] = arr[j];
        }
    }

    return res;
}


template <typename T>
std::vector<T> vector_value_minimum_limiter(const std::vector<T>& arr1, T value1){
    std::vector<T> res(arr1.size());

    for (int i = 0; i < arr1.size(); i++){
        if (arr1[i] < value1){
            res[i] = value1;
        } else {
            res[i] = arr1[i];
        }
    }

    return res;
}



// @overload
template <typename T>
std::vector<std::vector<T>> vector_value_minimum_limiter(const std::vector<std::vector<T>>& arr1, T value1){
    std::vector<std::vector<T>> res(arr1.size(), std::vector<T>(arr1[0].size()));

    for (int i = 0; i < arr1.size(); i++){
        for (int j=0; j < arr1[0].size(); j++){
            if (arr1[i][j] < value1){
                res[i][j] = value1;
            } else {
                res[i][j] = arr1[i][j];
            }
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
void print_vector(const std::vector<T>& vec) {
    std::cout << "[ ";
    for (int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << "] " << std::endl;
}

template <typename T>
void print_vector(const std::vector<std::vector<T>>& vec) {
    std::cout << "[ ";
    for (int i = 0; i < vec.size(); i++){
        for (int j=0; j<vec[0].size(); j++){
            std::cout << vec[i][j] << " ";
        }
        std::cout << "" << std::endl;      
    }
    std::cout << "] " << std::endl;
}


template <typename T>
void file_writer_vector(T time, const std::vector<T>& vec1, std::string fileName, int printflag1){
    std::ofstream output_file(fileName, std::ios::app);
    
    if (output_file.is_open()) {
        output_file << std::fixed << std::setprecision(15) << time << ",";
        for (int i = 0; i < vec1.size(); i++){
            output_file <<        std::fixed << std::setprecision(15) << vec1[i] << "," ;
        }
        output_file << std::endl;
        output_file.close();

        if (printflag1 == 1){
            std::cout << std::unitbuf << "File written successfully at timestep: " << time << "\n";
        }
    } else {
        std::cerr << "Error opening file for writing.\n";
    }
}




// Function to read and parse a JSON file using standard library functions.
// Takes the filename as input and populates the JsonData structure.
bool readAndParseJson(const std::string& filename, JsonData& data);

std::string make_filename(std::string folder, std::string type, std::string run_num);

void clear_file(const std::string& fileName);


#endif // UTILS_SOLVER_H