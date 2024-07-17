/**
 * Some utility functions
 */

#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <vector>

enum class Method : int { EXPLICIT = 0, IMPLICIT = 2, CRANK_NICOLSON = 1 };

/**
 * Write the solution to a file
 */
void write_to_file(const std::vector<double> &u, int n, const std::string &method_str);

/**
 * Parse the command line arguments to get the number of grid points and the method
 */
void parse_argument(int argc, char *argv[], int &n, Method &method);

/**
 * Draw a progress bar
 * @param current the current progress
 * @param total the total progress
 */
void progress_bar(int current, int total);

#endif // UTIL_H
