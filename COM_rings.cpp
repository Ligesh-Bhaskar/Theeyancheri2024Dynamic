#include <iostream>
#include <fstream>
#include <vector>

int main() {
    const int Nstep = 400000, Np = 4000, Ntot = Nstep * Np, Nr = 50;
    std::vector<float> P(Np), Q(Np), M(Np);
    std::vector<float> t1(Ntot), t2(Ntot);

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

    for (int i = 0; i < Ntot; ++i) {
        dataFile >> t1[i] >> t2[i];
    }

    for (int k = 0; k < Ntot; k += Np) {
        for (int j = 0; j < Np; j += Nr) {
            for (int i = 0; i < Nr; ++i) {
                P[j] += t1[i + j + k - 2];
                Q[j] += t2[i + j + k - 2];
                // M[j] += Rz[i + j + k - 2];
            }

            outputFile << P[j] / Nr << " " << Q[j] / Nr << std::endl; // << " " << M[j] / Nr << std::endl;
        }
    }

    dataFile.close();
    outputFile.close();

    return 0;
}

