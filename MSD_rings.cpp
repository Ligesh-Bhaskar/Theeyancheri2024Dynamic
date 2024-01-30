#include <iostream>
#include <fstream>
#include <vector>

int main() {
    const int NSTEP = 2000001, NPART = 80, Ntau = 1000000;
    std::vector<std::vector<double>> X(NSTEP, std::vector<double>(NPART)),
                                      Y(NSTEP, std::vector<double>(NPART)),
                                      Z(NSTEP, std::vector<double>(NPART)),
                                      ID(NSTEP, std::vector<double>(NPART)),
                                      TYP(NSTEP, std::vector<double>(NPART));

    std::ifstream inputFile("name.dat");
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    std::string input, output;
    inputFile >> input >> output;
    inputFile.close();

    std::ifstream dataFile(input);
    if (!dataFile.is_open()) {
        std::cerr << "Error opening data file." << std::endl;
        return 1;
    }

    std::ofstream outputFile(output);

    for (int j1 = 0; j1 < NSTEP; ++j1) {
        for (int j2 = 0; j2 < NPART; ++j2) {
            dataFile >> X[j1][j2] >> Y[j1][j2] >> Z[j1][j2];
        }
    }

    for (int l = 0; l < Ntau; ++l) {
        double tot = 0.0;

        for (int j2 = 0; j2 < NPART; ++j2) {
            double dist = 0.0;

            for (int k = 0; k < NSTEP - l; ++k) {
                double dx = X[l + k][j2] - X[k][j2];
                double dy = Y[l + k][j2] - Y[k][j2];
                // double dz = Z[l + k][j2] - Z[k][j2];

                dist += dx * dx + dy * dy; // + dz * dz;
            }

            tot += dist / static_cast<double>(NSTEP - l);
        }

        outputFile << static_cast<double>(l) << " " << tot / static_cast<double>(NPART) << std::endl;
        // std::cout << log10(static_cast<double>(l)) << " " << log10(tot / static_cast<double>(NPART)) << std::endl;
    }

    dataFile.close();
    outputFile.close();

    return 0;
}

