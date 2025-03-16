// json_parser.cpp
#include "json_parser.h"
#include <cctype> // For isspace and isdigit

// Helper function to trim leading and trailing whitespace from a string.
std::string trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t\n\r");
    if (std::string::npos == first) {
        return s;
    }
    size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, (last - first + 1));
}

// Helper function to check if a string can be converted to an integer.
bool isInteger(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) {
        return false;
    }
    for (size_t i = 1; i < s.length(); ++i) {
        if (!isdigit(s[i])) {
            return false;
        }
    }
    return true;
}

// Helper function to check if a string can be converted to a floating-point number.
bool isFloat(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+') && (s[0] != '.'))) {
        return false;
    }
    bool decimalPointFound = false;
    for (size_t i = 0; i < s.length(); ++i) {
        if (isdigit(s[i])) continue;
        if (s[i] == '.' && !decimalPointFound) {
            decimalPointFound = true;
            continue;
        }
        return false;
    }
    return true;
}

bool readAndParseJson(const std::string& filename, JsonData& data) {
    // Open the JSON file for reading.
    std::ifstream inputFile(filename);

    // Check if the file was opened successfully.
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // Read the entire content of the file into a string stream.
    std::stringstream buffer;
    buffer << inputFile.rdbuf();
    std::string jsonContent = buffer.str();

    // Close the input file.
    inputFile.close();

    // Remove leading and trailing whitespace from the JSON content.
    jsonContent = trim(jsonContent);

    // Check if the content starts and ends with curly braces, indicating a JSON object.
    if (jsonContent.empty() || jsonContent.front() != '{' || jsonContent.back() != '}') {
        std::cerr << "Error: Invalid JSON format in file " << filename << std::endl;
        return false;
    }

    // Remove the outer curly braces.
    jsonContent = jsonContent.substr(1, jsonContent.length() - 2);

    // Use a string stream to easily process the content.
    std::stringstream contentStream(jsonContent);
    std::string token;

    // Loop through the key-value pairs in the JSON object.
    while (std::getline(contentStream, token, ',')) {
        // Trim whitespace from the token.
        token = trim(token);

        // Find the position of the colon, which separates the key and the value.
        size_t colonPos = token.find(':');
        if (colonPos == std::string::npos) {
            std::cerr << "Error: Invalid JSON format - missing colon in '" << token << "'" << std::endl;
            return false;
        }

        // Extract the key and the value.
        std::string key = trim(token.substr(0, colonPos));
        std::string value = trim(token.substr(colonPos + 1));

        // Remove quotes from the key if it's a string.
        if (key.length() >= 2 && key.front() == '"' && key.back() == '"') {
            key = key.substr(1, key.length() - 2);
        } else {
            std::cerr << "Error: Invalid JSON format - key not enclosed in quotes in '" << token << "'" << std::endl;
            return false;
        }

        // Remove quotes from the value if it's a string.
        if (value.length() >= 2 && value.front() == '"' && value.back() == '"') {
            data.stringValues[key] = value.substr(1, value.length() - 2);
        } else if (isInteger(value)) {
            try {
                data.integerValues[key] = std::stoi(value);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid integer value for key '" << key << "': " << e.what() << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Error: Integer value out of range for key '" << key << "': " << e.what() << std::endl;
                return false;
            }
        } else if (isFloat(value)) {
            try {
                data.floatValues[key] = std::stod(value);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid floating-point value for key '" << key << "': " << e.what() << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Error: Floating-point value out of range for key '" << key << "': " << e.what() << std::endl;
                return false;
            }
        } else {
            std::cerr << "Error: Invalid JSON value type for key '" << key << "': " << value << std::endl;
            return false;
        }
    }

    return true;
}