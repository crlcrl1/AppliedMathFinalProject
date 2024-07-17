#include "util.h"

#include <chrono>
#include <format>
#include <fstream>
#include <iostream>


using namespace std;
using namespace chrono;

void write_to_file(const vector<double> &u, const int n, const string &method_str) {
    ofstream file(format("output_{}_{}.txt", n, method_str));
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            file << format("{:.5e} ", u[i * (n - 1) + j]);
        }
        file << endl;
    }
    file.close();
}

void print_usage() {
    cout << "Usage: ./main -n <n> -m <method>" << endl;
    cout << "n: the number of grid points" << endl;
    cout << "method: 0 for explicit, 1 for Crank-Nicolson, 2 for implicit" << endl;
    exit(1);
}

void parse_argument(const int argc, char *argv[], int &n, Method &method) {
    if (argc != 5) {
        print_usage();
    }
    if (strcmp(argv[1], "-n") != 0 || strcmp(argv[3], "-m") != 0) {
        print_usage();
    }
    char *end;
    n = strtol(argv[2], &end, 10);
    if (*end != '\0' || n <= 0) {
        print_usage();
    }
    const int m = strtol(argv[4], &end, 10);
    if (*end != '\0' || m < 0 || m > 2) {
        print_usage();
    }
    method = static_cast<Method>(m);
}

void progress_bar(const int current, const int total) {
    if (current != 0) {
        cout << "\033[K\033[1G";
    }
    static time_point<steady_clock> lastCallTime;
    static int lastItemCount = 0;
    constexpr int barWidth = 70;
    cout << "[";
    const int pos = static_cast<int>(barWidth * current / static_cast<double>(total));
    const auto bar = string(pos, '=') + '>' + string(barWidth - pos - 1, ' ');
    cout << bar;
    cout << format("] {}%", static_cast<int>(100.0 * current / total));
    if (lastItemCount != 0) {
        const auto now = high_resolution_clock::now();
        const auto duration = duration_cast<milliseconds>(now - lastCallTime).count();
        const int itemCount = current - lastItemCount;
        cout << format("  {:.2f} it/s", itemCount / (static_cast<double>(duration) / 1000.0));
        lastCallTime = now;
        lastItemCount = current;
    }
    else {
        cout << "  0.00 it/s";
        lastCallTime = high_resolution_clock::now();
        lastItemCount = current;
    }
    cout.flush();
}
