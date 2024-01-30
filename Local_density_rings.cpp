#include <iostream>
#include <fstream>

const int Nstep = 4000;
const int ilines = 10;

int main() {
    std::string input, output;
    double x[Nstep], y[Nstep];

    std::ifstream inputFile("name.dat");
    inputFile >> input >> output;
    inputFile.close();

    std::ifstream dataFile(input);
    std::ofstream outputFile(output);

    for (int k = 0; k < Nstep; ++k) {
        dataFile >> x[k] >> y[k];
    }

    for (int i = -35; i <= 35 - ilines; i += ilines) {
        for (int j = -35; j <= 35 - ilines; j += ilines) {
            int m = i + ilines;
            int n = j + ilines;
            double nlines = 0.0;

            for (int l = 0; l < Nstep; ++l) {
                if ((x[l] >= i && x[l] <= m) && (y[l] >= j && y[l] <= n)) {
                    nlines += 1.0;
                }
            }

            outputFile << (nlines * 0.785) / (ilines * ilines) << std::endl;
        }
    }

    dataFile.close();
    outputFile.close();

    return 0;
}
