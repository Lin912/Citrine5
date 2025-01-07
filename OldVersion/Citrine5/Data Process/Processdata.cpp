#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>

const int KEEP_ROWS = 2100;

bool isNumber(const std::string& str) {
    try {
        std::stod(str);
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    }
    return true;
}

std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    
    if (!file.is_open()) {
        std::cerr << "Can not Open the file: " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;
        
        while (std::getline(ss, value, ',')) {
            if (isNumber(value)) {
                row.push_back(std::stod(value));
            } else {
                row.push_back(0.0);
            }
        }
        
        data.push_back(row);
    }
    
    file.close();
    return data;
}

void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Can not Open the file: " << filename << std::endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }
    
    file.close();
}

int main() {
    std::string inputFile = "output.csv";
    std::string outputFile = "output1.csv";

    std::vector<std::vector<double>> data = readCSV(inputFile);

    size_t actualRows = data.size();
    size_t actualCols = actualRows > 0 ? data[0].size() : 0;

    if (actualRows < KEEP_ROWS || actualCols == 0) {
        std::cerr << "The CSV file does not have enough rows or is not in the correct format! Actual number of rows: " << actualRows << ", Actual number of rows: " << actualCols << std::endl;
        return 1;
    }

    for (size_t i = KEEP_ROWS; i < actualRows; ++i) {
        for (size_t j = 0; j < actualCols; ++j) {
            data[i][j] = 0.0;
        }
    }

    writeCSV(outputFile, data);

    std::cout << "PorcessData Done! " << outputFile << std::endl;

    return 0;
}
