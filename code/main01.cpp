// main.cpp
#include "json_parser.h"
#include <iostream>

int main() {
    JsonData config;
    std::string filename = "D:/actual files in this 2/to_sort/computational_physics/frojekt/dummyTheMummy/code/configuration.json"; // Replace with your JSON filename

    if (readAndParseJson(filename, config)) {
        std::cout << "JSON file parsed successfully!" << std::endl;

        // Access the parsed values:
        if (config.stringValues.count("stringValue")) {
            std::cout << "String Value: " << config.stringValues["stringValue"] << std::endl;
        }

        if (config.stringValues.count("anotherString")) {
            std::cout << "Another String Value: " << config.stringValues["anotherString"] << std::endl;
        }

        if (config.floatValues.count("floatValue")) {
            std::cout << "Float Value: " << config.floatValues["floatValue"] << std::endl;
        }

        if (config.integerValues.count("intValue")) {
            std::cout << "Integer Value: " << config.integerValues["intValue"] << std::endl;
        }

        if (config.floatValues.count("anotherFloat")) {
            std::cout << "Float Value: " << config.floatValues["anotherFloat"] << std::endl;
        }

        // Access other values similarly based on your JSON structure.
    } else {
        std::cerr << "Failed to parse JSON file." << std::endl;
    }

    return 0;
}