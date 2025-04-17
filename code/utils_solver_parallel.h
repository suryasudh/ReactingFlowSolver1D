#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <thread>
#include <functional>
#include <vector>
#include <list>
#include <array>
#include <numeric> 
#include <cmath>
#include <iomanip>
#include <unordered_map>

// One problem with multithreading is that the overhead is more than the problem
template <typename T>
std::vector<T> add_vectors(const std::vector<T>& arr1, const std::vector<T>& arr2) {
    std::vector<T> res(arr1.size());

    const int n_threads = 16; //std::thread::hardware_concurrency(); // or set manually
    const std::size_t N = arr1.size();
    const std::size_t chunk_size = (N + n_threads - 1) / n_threads;

    auto worker = [&](std::size_t start, std::size_t end) {
        for (std::size_t i = start; i < end && i < N; ++i) {
            res[i] = arr1[i] + arr2[i];
        }
    };

    std::vector<std::thread> threads;
    for (int t = 0; t < n_threads; ++t) {
        std::size_t start = t * chunk_size;
        std::size_t end = start + chunk_size;
        threads.emplace_back(worker, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return res;
}