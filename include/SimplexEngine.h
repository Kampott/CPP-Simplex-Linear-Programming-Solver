#pragma once
#include <Matrix.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#define N 10000000



class SimplexEngine {
public:
    std::vector<double> function, function_copy, restrictions, restrictions_copy;
    std::vector<int> operators, operators_copy, positive_coefs, positive_coefs_copy, negative_coefs, negative_coefs_copy, basis_ex, artificial_basis;
    std::vector<int> leq, leq_copy;
    std::vector<double> leq_column, leq_column_copy;
    std::vector<std::pair<int, int>> expanded_vars;
    std::string taskType, taskType_copy;
    double function_constant, function_constant_copy;
    int doubletrigger = 0;
    int r, r2, c, c2;
    Matrix mate, mate_copy;
    SimplexEngine(const std::string& filename);
    long long fact(int n)
    {
        if (n == 0)
            return 1;
        long long res = 1;
        for (int i = 2; i <= n; i++) {
            res = res * i;
        }
        return res;
    }

    void printData(const Matrix& mate, const std::vector<int>& operators, const std::vector<double>& function, const std::vector<double>& restrictions, const std::vector<int>& positive_coefs, const double& function_constant) const;

    void Canonize(std::vector<double>& canon_function, std::vector<double>& canon_restrictions, Matrix& canon_matrix, std::vector<int>& canon_operators, std::vector<int>& canon_positive_coefs, std::vector<double>& canon_negative_coefs, std::vector<std::pair<int, int>>& canon_expanded_vars, std::vector<int>& canon_basis_ex, std::vector<double>& canon_leq_column, double& canon_function_constant, std::vector<int>& canon_leq, const std::string& canon_TastkType, const bool& isDouble, const bool& v);
    void artificialBasis(Matrix& artificial_mate, std::vector<int>& artificial_leq, std::vector<int>& artificial_basis_ex, std::vector<double>& artificial_function, std::vector<int>& artificial_artificial_basis, std::vector<int>& artificial_positive_coefs, const std::vector<int>& artificial_operators, const std::vector<double>& artificial_restrictions, const double& artificial_function_constant, const double& v);
};