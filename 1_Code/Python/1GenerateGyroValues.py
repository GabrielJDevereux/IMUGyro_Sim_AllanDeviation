import csv
import numpy as np

class GyroscopeDeviationSimulator:
    def __init__(self, samples=100000, F_s=100):
        # Sim parameters
        self.samples = samples  
        self.F_s = F_s  # in Hz
        self.T_s = 1 / self.F_s  # Sampling period
        
        # Noise generation deviation standards
        self.arw_std = 0.001  
        self.rrw_std = 0.0001  
        self.bias_instability_std = 0.00001  
        
        #Bias and rate random walk variables
        self.bias_x = self.bias_y = self.bias_z = 0.0
        self.rate_rw_x = self.rate_rw_y = self.rate_rw_z = 0.0
        
        # Arrays for gyro data XYZ
        self.gyro_x = np.zeros(self.samples)
        self.gyro_y = np.zeros(self.samples)
        self.gyro_z = np.zeros(self.samples)

    # Generate noise based on (ARW) / white noise
    def generate_angle_random_walk(self):
        arw_noise_x = np.random.normal(0, self.arw_std, self.samples)
        arw_noise_y = np.random.normal(0, self.arw_std, self.samples)
        arw_noise_z = np.random.normal(0, self.arw_std, self.samples)
        return arw_noise_x, arw_noise_y, arw_noise_z

    def apply_rate_random_walk(self, rate_rw, rrw_std, t):
        # Apply Rate Random Walk (RRW) as a random drift in the rate
        rrw_std = rrw_std * (1 + 0.000005 * t)  # Reduced increase over time
        rate_rw += np.random.normal(0, rrw_std) * np.sqrt(self.T_s)
        return rate_rw

    def apply_bias_instability(self, bias, bias_instability_std, t):
        # Apply Bias Instability 
        bias_instability_std = bias_instability_std * (1 + 0.0000025 * t)  # Reduced increase over time
        bias += np.random.normal(0, bias_instability_std) * np.sqrt(self.T_s)
        return bias

    # Sim gyro values with ARW RRW and BI
    def simulate_gyro_values(self):
        arw_noise_x, arw_noise_y, arw_noise_z = self.generate_angle_random_walk()
        
        # Loop through to apply RRW and BI
        for t in range(self.samples):
            self.rate_rw_x = self.apply_rate_random_walk(self.rate_rw_x, self.rrw_std, t)
            self.rate_rw_y = self.apply_rate_random_walk(self.rate_rw_y, self.rrw_std, t)
            self.rate_rw_z = self.apply_rate_random_walk(self.rate_rw_z, self.rrw_std, t)

            self.bias_x = self.apply_bias_instability(self.bias_x, self.bias_instability_std, t)
            self.bias_y = self.apply_bias_instability(self.bias_y, self.bias_instability_std, t)
            self.bias_z = self.apply_bias_instability(self.bias_z, self.bias_instability_std, t)

            self.gyro_x[t] = self.rate_rw_x + self.bias_x + arw_noise_x[t]
            self.gyro_y[t] = self.rate_rw_y + self.bias_y + arw_noise_y[t]
            self.gyro_z[t] = self.rate_rw_z + self.bias_z + arw_noise_z[t]

    def compute_deviation(self):
        # Troubleshooting / deviation for each axis
        deviation_x = np.std(self.gyro_x)
        deviation_y = np.std(self.gyro_y)
        deviation_z = np.std(self.gyro_z)

        print(f"Deviation Gyro X-axis: {deviation_x}")
        print(f"Deviation Gyro Y-axis: {deviation_y}")
        print(f"Deviation Gyro Z-axis: {deviation_z}")

    def write_to_csv(self, filename="gyro_data.csv"):
        # Write xyz values to a csv file
        with open(filename, mode='w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Gyro_X', 'Gyro_Y', 'Gyro_Z']) 
            for i in range(self.samples):
                csvwriter.writerow([self.gyro_x[i], self.gyro_y[i], self.gyro_z[i]])
        print(f"Data written to {filename}")

def main():
    gyro_sim = GyroscopeDeviationSimulator()
    gyro_sim.simulate_gyro_values()
    gyro_sim.compute_deviation()  # Compute and print the deviation
    gyro_sim.write_to_csv("gyro_data.csv")

if __name__ == "__main__":
    main()





   