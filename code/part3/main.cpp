/**
 * Solve the heat equation using the LDL^T decomposition.
 *
 * To compile the code, you should have cmake installed and a C++ compiler that supports C++20.
 * Run the following commands to compile the code:
 * cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
 * cmake --build build
 */

#include <cassert>
#include <chrono>
#include <cmath>
#include <format>
#include <iostream>
#include <vector>

#include "util.h"

// The value of PI. On some platforms, the macro M_PI may not be defined.
constexpr double PI = 3.14159265358979323846;


using namespace std;
using namespace chrono;

/**
 * A class to view the memory-optimized matrix as a normal matrix.
 * It is used to map the index between the memory-optimized matrix and the normal matrix.
 */
class MatrixView {
    vector<vector<double>> *original;
    size_t dim;
    size_t stripSize;

public:
    explicit MatrixView(vector<vector<double>> &original) : original(&original) {
        dim = original.size();
        assert(original[0].size() % 2 == 1);
        stripSize = (original[0].size() - 1) / 2;
    }

    /**
     * Get the element at the position (i, j) in the memory-optimized matrix
     */
    double &operator()(const int i, const int j) const { return (*original)[i][j - i + stripSize]; }

    size_t getDim() const { return dim; }

    size_t getStripSize() const { return stripSize; }
};


/**
 * LDL^T decomposition.
 *
 * @param A the coefficient matrix. This matrix will be modified in place to store the L matrix.
 * @param v the vector to store the diagonal elements of the matrix D.
 */
void LDL_decomposition(vector<vector<double>> &A, vector<double> &v) {
    const MatrixView view(A);
    const int d = static_cast<int>(view.getDim());
    const int n = static_cast<int>(view.getStripSize()) + 1;
    for (int j = 0; j < d; ++j) {
        const int temp2 = max(j - n + 1, 0);
        double temp = 0;
        for (int i = temp2; i <= j - 1; ++i) {
            v[i] = view(j, i) * view(i, i);
            temp += view(j, i) * v[i];
        }
        view(j, j) -= temp;
        const int temp1 = min(d, j + n);
        for (int i = j + 1; i < temp1; ++i) {
            temp = 0;
            const int start = max(i - n + 1, temp2);
            const int end = min(j - 1, i + n - 1);
            for (int k = start; k <= end; ++k) {
                temp += view(i, k) * v[k];
            }
            view(i, j) -= temp;
            view(i, j) /= view(j, j);
        }
    }

    for (int i = 0; i < d; ++i) {
        v[i] = view(i, i);
        view(i, i) = 1;
    }
    for (int i = 0; i < d; ++i) {
        const int temp = min(d, i + n);
        for (int j = i + 1; j < temp; ++j) {
            view(i, j) = view(j, i);
        }
    }
}

/**
 * Forward elimination after the LDL^T decomposition
 */
void forward_elimination(vector<vector<double>> &A, vector<double> &b) {
    const MatrixView view(A);
    const int d = static_cast<int>(view.getDim());
    const int n = static_cast<int>(view.getStripSize()) + 1;
    for (int i = 0; i < d - 1; ++i) {
        b[i] /= view(i, i);
        const int temp = min(d, i + n);
        for (int j = i + 1; j < temp; ++j) {
            b[j] -= view(j, i) * b[i];
        }
    }
    b[d - 1] /= view(d - 1, d - 1);
}

/**
 * Backward elimination after the LDL^T decomposition
 */
void backward_elimination(vector<vector<double>> &A, vector<double> &b) {
    const auto view = MatrixView(A);
    const int d = static_cast<int>(view.getDim());
    const int n = static_cast<int>(view.getStripSize()) + 1;
    for (int i = d - 1; i > 0; --i) {
        b[i] /= view(i, i);
        const int temp = max(-1, i - n);
        for (int j = i - 1; j > temp; --j) {
            b[j] -= view(j, i) * b[i];
        }
    }
    b[0] /= view(0, 0);
}

/**
 * Solve the linear system after the LDL^T decomposition
 *
 * @param A The coefficient matrix modified by the LDL^T decomposition
 * @param v The diagonal elements of the matrix D
 * @param b The right-hand side vector to be solved
 * @return The solution of the linear system
 */
vector<double> LDL_solution(vector<vector<double>> &A, const vector<double> &v, vector<double> &b) {
    const size_t d = b.size();
    forward_elimination(A, b);
    for (int i = 0; i < d; ++i) {
        b[i] /= v[i];
    }
    backward_elimination(A, b);
    return b;
}

/**
 * Generate the coefficient matrix A
 */
vector<vector<double>> generate_A(const int n, const Method method, const int t_n) {
    const double dt = 1.0 / t_n;
    const double dx = 1.0 / n;
    const int m = static_cast<int>(method);
    const double delta = dt / (dx * dx) * m / 2;
    const double a = 4.0 / 9 + delta * 8 / 3.0;
    const double b = 1.0 / 9 - delta * 1 / 3.0;
    const double c = 1.0 / 36 - delta * 1 / 3.0;
    auto A = vector((n - 1) * (n - 1), vector(2 * n + 1, 0.0));
    const MatrixView view(A);
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            const int index = i * (n - 1) + j;
            view(index, index) = a;
            if (i != 0) {
                view(index, index - n + 1) = b;
            }
            if (i != n - 2) {
                view(index, index + n - 1) = b;
            }
            if (j != 0) {
                view(index, index - 1) = b;
            }
            if (j != n - 2) {
                view(index, index + 1) = b;
            }
            if (i != 0 && j != 0) {
                view(index, index - n) = c;
            }
            if (i != 0 && j != n - 2) {
                view(index, index - n + 2) = c;
            }
            if (i != n - 2 && j != 0) {
                view(index, index + n - 2) = c;
            }
            if (i != n - 2 && j != n - 2) {
                view(index, index + n) = c;
            }
        }
    }
    return A;
}

/**
 * Generate the initial values according to the initial condition
 */
vector<double> generate_initial(const int n) {
    auto u = vector((n - 1) * (n - 1), 0.0);
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1; j++) {
            u[i * (n - 1) + j] = sin(PI * (i + 1) / n) * sin(PI * (j + 1) / n);
        }
    }
    return u;
}

/**
 * Generate the right-hand side vector b for each time step
 */
vector<double> generate_b(const int n, const Method method, const vector<double> &u,
                          const int t_n) {
    auto res = vector((n - 1) * (n - 1), 0.0);
    const double dx = 1.0 / n;
    const int m = 2 - static_cast<int>(method);
    const double dt = 1.0 / t_n;
    const double delta = dt / (dx * dx) * m / 2;
    const double a = 4.0 / 9 - delta * 8 / 3.0;
    const double b = 1.0 / 9 + delta * 1 / 3.0;
    const double c = 1.0 / 36 + delta * 1 / 3.0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            const int index = i * (n - 1) + j;
            res[index] = a * u[index];
            if (i != 0) {
                res[index] += b * u[index - n + 1];
            }
            if (i != n - 2) {
                res[index] += b * u[index + n - 1];
            }
            if (j != 0) {
                res[index] += b * u[index - 1];
            }
            if (j != n - 2) {
                res[index] += b * u[index + 1];
            }
            if (i != 0 && j != 0) {
                res[index] += c * u[index - n];
            }
            if (i != 0 && j != n - 2) {
                res[index] += c * u[index - n + 2];
            }
            if (i != n - 2 && j != 0) {
                res[index] += c * u[index + n - 2];
            }
            if (i != n - 2 && j != n - 2) {
                res[index] += c * u[index + n];
            }
        }
    }
    return res;
}

/**
 * Solve the heat equation
 */
vector<double> solve(const int n, const Method method, const int t_n) {
    auto A = generate_A(n, method, t_n);
    auto u = generate_initial(n);
    auto v = vector((n - 1) * (n - 1), 0.0);
    LDL_decomposition(A, v);
    cout << endl;
    cout << "\033[?25l";
    for (int t = 0; t < t_n; ++t) {
        if (t % 100 == 0) {
            progress_bar(t, t_n);
        }
        auto b = generate_b(n, method, u, t_n);
        u = LDL_solution(A, v, b);
    }
    cout << "\033[?25h";
    return u;
}


int main(int argc, char *argv[]) {
    int n;
    Method method;
    parse_argument(argc, argv, n, method);
    cout << "n = " << n << endl;
    const auto method_str = method == Method::EXPLICIT         ? "explicit"
                            : method == Method::CRANK_NICOLSON ? "Crank-Nicolson"
                                                               : "implicit";
    cout << format("Method: {}\n", method_str);
    const auto start = high_resolution_clock::now();
    const auto u = solve(n, method, 12 * n * n);
    const auto end = high_resolution_clock::now();
    cout << endl;
    const auto duration =
            static_cast<double>(duration_cast<milliseconds>(end - start).count()) / 1000.0;
    cout << format("Time: {:.4f}s\n", duration);
    cout << format("Error on the center = {:.5e}\n", abs(u[n * (n / 2 - 1)] - exp(-2 * PI * PI)));
    write_to_file(u, n, method_str);
    return 0;
}
