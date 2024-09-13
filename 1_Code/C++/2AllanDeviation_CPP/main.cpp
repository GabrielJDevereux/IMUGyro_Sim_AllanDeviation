#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unistd.h>

// Load gyro values from csv file
class LoadGyroValues {
public:
    std::vector<double> gyro_x, gyro_y, gyro_z;

    void load_from_csv(const std::string& csv_file) {
        std::ifstream file(csv_file);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << csv_file << std::endl;
            return;
        }

        std::string line;
        std::getline(file, line); // skip header

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string item;

            std::getline(ss, item, ',');
            gyro_x.push_back(std::stod(item));

            std::getline(ss, item, ',');
            gyro_y.push_back(std::stod(item));

            std::getline(ss, item, ',');
            gyro_z.push_back(std::stod(item));

            if (gyro_x.size() == 100000) {
                break; // limit samples
            }
        }
        file.close();
    }
};

// Class to calculate Allan Deviation
class AllanDeviation {
public:
    int F_s = 100; // Hz
    double T_s;
    std::vector<double> theta_x;
    std::vector<double> allanvar, tau, deviation;

    AllanDeviation() {
        T_s = 1.0 / F_s;
    }

    // Perform Euler integration on angular velocities.
    void euler_integration(const std::vector<double>& x) {
        theta_x.clear();
        double sum = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            sum += x[i] * T_s;
            theta_x.push_back(sum);
        }
    }

    // Calculate Allan Variance and then root for deviation.
    void allan_deviation() {
        size_t L = theta_x.size();
        size_t max_m = 32770;  // Adjust max_m
        allanvar.clear();
        tau.clear();

        for (size_t m = 1; m < max_m; ++m) {
            long double sum_sq_diff = 0;

            for (size_t k = 0; k < L - 2 * m; ++k) {
                double diff = theta_x[k + 2 * m] - 2 * theta_x[k + m] + theta_x[k];
                sum_sq_diff += diff * diff;
            }

            // Calculate tau directly
            double tau_m = m * T_s;

            long double variance = sum_sq_diff / (2.0 * tau_m * tau_m * (static_cast<double>(L) - 2.0 * static_cast<double>(m)));

            allanvar.push_back(variance);
            tau.push_back(tau_m);  // Store tau

            // Debug output for the specific value of m
            if (m == 32768) {
                std::cout << "At m = " << m << ", sum_sq_diff = " << sum_sq_diff << std::endl;
                std::cout << "m = " << m << ", tau = " << tau.back() << ", variance = " << variance << std::endl;
            }
        }

        // Allan deviation
        deviation.clear();
        for (double var : allanvar) {
            deviation.push_back(sqrt(var));
        }

        // Print results for debugging
        for (size_t i = 0; i < deviation.size(); ++i) {
            if (i % 100 == 0 || std::abs(tau[i] - 327.68) < 0.01) {  // Print every 100th value and when tau is around 327.68
                std::cout << "Final Allan deviation at tau = " << tau[i] << " is " << deviation[i] << std::endl;
            }
        }
    }


    // Write the results to a CSV file
    void write_to_csv(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        file << "Averaging Time (tau),Allan Deviation\n";
        for (size_t i = 0; i < tau.size(); ++i) {
            file << tau[i] << "," << deviation[i] << "\n";
        }

        file.close();
        std::cout << "Results written to " << filename << std::endl;
    }
};

int main() {
    // Directory check
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != __nullptr) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        std::cerr << "Failed to get current directory" << std::endl;
    }
    // Load gyro values from csv file
    std::string csv_file = "gyro_data.csv";
    LoadGyroValues gen_gyro;
    gen_gyro.load_from_csv(csv_file);

    // Check if data is loaded correctly
    if (gen_gyro.gyro_x.empty()) {
        std::cerr << "Failed to load data from CSV file." << std::endl;
        return 1;
    }

    // Instant Allan Deviation class and compute deviation
    AllanDeviation allan_dev_calc;
    allan_dev_calc.euler_integration(gen_gyro.gyro_x);  // Using Gyro x axis data
    allan_dev_calc.allan_deviation();

    // Write the Allan deviation results to a CSV file
    std::string output_file = "allan_deviation_output.csv";
    allan_dev_calc.write_to_csv(output_file);

    return 0;
}






