#include <iostream>
#include <fstream>
#include <vector>

int main() {
    const int N = 50;
    const int n_frames = 2000001;
    const double mass = 1.0;
    std::vector<std::vector<double>> coord(2, std::vector<double>(N));
    std::vector<std::vector<double>> coord_M(2, std::vector<double>(N));
    std::vector<double> CM(2), gyr_ten(4);

    std::ifstream inputFile("Pos_AF100_ChainL70.dat");
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    std::ofstream outputFile("Gyration_Tensor_Active_ChainL70_F100_L50_Topol_80Chains_WCA_K0_Long.dat");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

    for (int j = 0; j < n_frames; ++j) {
        for (int i = 0; i < N; ++i) {
            inputFile >> coord[0][i] >> coord[1][i];
            coord_M[0][i] = coord[0][i] * mass;
            coord_M[1][i] = coord[1][i] * mass;
        }

        // Center of mass calculation
        CM[0] = std::accumulate(coord_M[0].begin(), coord_M[0].end(), 0.0) / N;
        CM[1] = std::accumulate(coord_M[1].begin(), coord_M[1].end(), 0.0) / N;

        // Gyration tensor calculation
        gyr_ten[0] = gyr_ten[1] = gyr_ten[2] = gyr_ten[3] = 0.0;
        for (int i = 0; i < N; ++i) {
            gyr_ten[0] += (coord[0][i] - CM[0]) * (coord[0][i] - CM[0]);
            gyr_ten[1] += (coord[1][i] - CM[1]) * (coord[0][i] - CM[0]);
            gyr_ten[2] += (coord[0][i] - CM[0]) * (coord[1][i] - CM[1]);
            gyr_ten[3] += (coord[1][i] - CM[1]) * (coord[1][i] - CM[1]);
        }

        for (int k = 0; k < 4; ++k) {
            gyr_ten[k] /= N;
            outputFile << std::setw(15) << std::setprecision(5) << std::scientific << gyr_ten[k];
        }
        outputFile << '\n';
    }

    inputFile.close();
    outputFile.close();

    return 0;
}
