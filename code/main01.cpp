// main.cpp
#include "json_parser.h"
#include <iostream>

int main() {
    JsonData config;
    std::string filename = "../configuration.json"; // Replace with your JSON filename

    if (readAndParseJson(filename, config)) {
        std::cout << "configuration.json file parsed successfully!" << std::endl;

        // Access the parsed values:
        if (config.stringValues.count("data_folder01")) {
            std::cout << "Path to cantera files: " << config.stringValues["data_folder01"] << std::endl;
        }

        if (config.integerValues.count("run_num")) {
            std::cout << "Run number chosen: " << config.integerValues["run_num"] << std::endl;
        }

        if (config.integerValues.count("float_precision")) {
            std::cout << "Precision chosen for the run: " << config.integerValues["float_precision"] << std::endl;
        }

        if (config.stringValues.count("anotherString")) {
            std::cout << "Another String Value: " << config.stringValues["anotherString"] << std::endl;
        }

        if (config.floatValues.count("floatValue")) {
            std::cout << "Float Value: " << config.floatValues["floatValue"] << std::endl;
        }

        // Access other values similarly based on your JSON structure.
    } else {
        std::cerr << "Failed to parse JSON file." << std::endl;
    }




    return 0;
}