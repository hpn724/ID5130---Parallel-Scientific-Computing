#ifndef VECTOR_OPERATIONS_PARALLEL_H
#define VECTOR_OPERATIONS_PARALLEL_H

#include <vector>

using namespace std;

vector<vector<double>> glorot_init(int nin, int nout);

template<typename T>
void show_vector(vector<T> A);

void show_matrix(vector<vector<double>> A);

template<typename T>
vector<T> flatten(const vector<vector<T>>& nested_vector);


vector<vector<double>> reshape(const vector<double>& input_vector, int rows, int cols);

double norm(vector<double> A);

vector<double> vector_diff(vector<double> A, vector<double> B);

double norm_diff(vector<double> A, vector<double> B);

vector<vector<double>> transpose(const vector<vector<double>>& matrix);

double dotProduct(vector<double> A, vector<double> B);

template<typename T>
double mean(vector<T> A);

vector<vector<double>> identity_matrix(int n);

vector<vector<double>> square_root_inverse_diag_matrix(vector<vector<double>> diag);

vector<vector<double>> matrix_sum(const vector<vector<double>>& A, const vector<vector<double>>& B);

vector<vector<double>> matrix_product(const vector<vector<double>> A, const vector<vector<double>> B);

#endif // VECTOR_OPERATIONS_H
