#include <cassert>
#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

constexpr double PI = 3.14159265358979323846;


using namespace std;
using namespace chrono;

/**
 * A class to view the memory-optimized matrix as a normal matrix
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

    double &operator()(const int i, const int j) const { return (*original)[i][j - i + stripSize]; }

    size_t getDim() const { return dim; }

    size_t getStripSize() const { return stripSize; }
};


/**
 * LDL^T decomposition
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

vector<double> LDL_solution(vector<vector<double>> &A, const vector<double> &v, vector<double> &b) {
    const size_t d = b.size();
    forward_elimination(A, b);
    for (int i = 0; i < d; ++i) {
        b[i] /= v[i];
    }
    backward_elimination(A, b);
    return b;
}

enum class Method : int { EXPLICIT = 0, IMPLICIT = 2, CRANK_NICOLSON = 1 };

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

vector<double> generate_initial(const int n) {
    auto u = vector((n - 1) * (n - 1), 0.0);
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1; j++) {
            u[i * (n - 1) + j] = sin(PI * (i + 1) / n) * sin(PI * (j + 1) / n);
        }
    }
    return u;
}

// Optimize
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

void progress_bar(const int current, const int total) {
    static time_point<steady_clock> lastCallTime;
    static int lastItemCount = 0;
    constexpr int barWidth = 70;
    cout << "[";
    const int pos = static_cast<int>(barWidth * current / static_cast<double>(total));
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            cout << "=";
        }
        else if (i == pos) {
            cout << ">";
        }
        else {
            cout << " ";
        }
    }
    cout << "] " << fixed << setprecision(2) << static_cast<int>(100.0 * current / total) << "%";
    if (lastItemCount != 0) {
        const auto now = high_resolution_clock::now();
        const auto duration = duration_cast<milliseconds>(now - lastCallTime).count();
        const int itemCount = current - lastItemCount;
        cout << "  " << itemCount / (static_cast<double>(duration) / 1000.0) << " it/s";
        lastCallTime = now;
        lastItemCount = current;
    }
    else {
        cout << "  0.00 it/s";
        lastCallTime = high_resolution_clock::now();
        lastItemCount = current;
    }
    cout << '\r';
    cout.flush();
}

vector<double> solve(const int n, const Method method, const int t_n) {
    auto A = generate_A(n, method, t_n);
    auto u = generate_initial(n);
    auto v = vector((n - 1) * (n - 1), 0.0);
    LDL_decomposition(A, v);
    cout << endl;
    for (int t = 0; t < t_n; ++t) {
        if (t % 100 == 0) {
            progress_bar(t, t_n);
        }
        auto b = generate_b(n, method, u, 1);
        u = LDL_solution(A, v, b);
    }
    return u;
}

void write_to_file(const vector<double> &u, const int n) {
    ofstream file(format("output_{}.txt", n));
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            file << scientific << setprecision(5) << u[i * (n - 1) + j] << " ";
        }
        file << endl;
    }
    file.close();
}

int main() {
    constexpr int n = 128;
    cout << "n = " << n << endl;
    const auto start = high_resolution_clock::now();
    const auto u = solve(n, Method::IMPLICIT, 10 * n * n);
    const auto end = high_resolution_clock::now();
    cout << endl;
    cout << "Time: " << fixed << setprecision(4)
         << static_cast<double>(duration_cast<milliseconds>(end - start).count()) / 1000.0 << "s"
         << endl;
    cout << scientific << setprecision(5) << u[n * (n / 2 - 1)] << endl;
    cout << scientific << setprecision(5)
         << "error = " << abs(u[n * (n / 2 - 1)] - exp(-2 * PI * PI)) << endl;
    write_to_file(u, n);
    return 0;
}
