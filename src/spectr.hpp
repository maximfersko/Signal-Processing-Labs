#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

class SpectrumCalculator {
public:
    SpectrumCalculator() {}

    void GetInputParameters() {
        std::cout << "Enter value for A: ";
        std::cin >> A_;
        std::cout << "Enter value for alpha: ";
        std::cin >> alpha_;
        std::cout << "Enter value for T: ";
        std::cin >> T_;
        std::cout << "Enter number of spectrum points: ";
        std::cin >> num_points_;
    }

    std::complex<double> ContinuousSpectrum(double omega) {
        return A_ / (alpha_ + std::complex<double>(0, omega));
    }

    std::complex<double> DiscreteSpectrum(int k) {
        int N = static_cast<int>(0.4 / T_); // tu = 0.4s
        std::complex<double> Xk = 0;

        for (int n = 0; n < N; ++n) {
            double t = n * T_;
            std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * k * n / N));
            Xk += A_ * std::exp(-alpha_ * t) * exp_term;
        }

        Xk *= T_;
        return Xk;
    }

    void CalculateDiscreteSpectrum() {
        double mse = 0.0;

        for (int k = 0; k < num_points_; ++k) {
            std::complex<double> Xk = DiscreteSpectrum(k);
            std::complex<double> Xak = ContinuousSpectrum(2 * M_PI * k / (num_points_ * T_));
            double error = std::abs(Xk - Xak);
            mse += error * error;
            std::cout << "Discrete Spectrum X(" << k << "): real = " << Xk.real() << ", imaginary = " << Xk.imag() << std::endl;
        }
        std::cout << std::endl;

        mse = std::sqrt(mse / num_points_);
        std::cout << "RMSE for discrete spectrum: " << mse << std::endl;
        std::cout << std::endl;
    }

    void CalculateContinuousSpectrum() {
        double mse = 0.0;
        double omega_step = 2 * M_PI / (num_points_ * T_);

        for (int i = 0; i < num_points_; ++i) {
            double omega = i * omega_step;
            std::complex<double> Xk = DiscreteSpectrum(0);
            std::complex<double> Xak = ContinuousSpectrum(omega);
            double error = std::abs(Xk - Xak);
            mse += error * error;
            std::cout << "Continuous Spectrum X(" << omega << "): real = " << Xak.real() << ", imaginary = " << Xak.imag() << std::endl;
        }
        std::cout << std::endl;

        mse = std::sqrt(mse / num_points_);
        std::cout << "RMSE for continuous spectrum: " << mse << std::endl;
    }

    void Calculate() {
        GetInputParameters();
        CalculateDiscreteSpectrum();
        CalculateContinuousSpectrum();
    }

private:
    double A_, alpha_, T_;
    int num_points_;
};