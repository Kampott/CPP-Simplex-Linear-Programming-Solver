#pragma once
#include <vector>
#include "SimplexEngine.h"
class SimplexSolver {
private:
    std::vector<double> restrictions, leq_column, negative_coefs, function;
    std::vector<int> basis_ex, artificial_basis, leq, operators, positive_coefs;
    std::vector<std::pair<int, int>> expanded_vars;
    SimplexEngine* method;
    std::string taskType;
    double function_constant;
    Matrix mate;
    bool isDouble, v;
public:
    SimplexSolver(const std::vector<double>& function_orig, const std::vector<double>& restrictions_orig, const std::vector<int>& operators_orig, const std::string& taskType_orig, const Matrix& mate_orig, SimplexEngine* method, const std::vector<int>& positive_coefs_orig, const double& function_constant_orig, const bool& isDouble, const bool& verbose);
    void Solve();
    void SolveDouble();
};