#include <SimplexEngine.h>
SimplexEngine::SimplexEngine(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("File not found.");
    }
    file >> mate.rows >> mate.cols;
    mate = Matrix(mate.rows, mate.cols);

    r = mate.rows;
    r2 = mate.rows;
    c = mate.cols;
    c2 = mate.cols;

    for (int i = 0; i < mate.rows; ++i) {
        for (int j = 0; j < mate.cols; ++j) {
            file >> mate.data[i][j];
        }
        std::string op;
        if (file >> op) {
            if (op == "!") { operators.push_back(0); }
            else {
                if (op == "=") operators.push_back(1);
                else if (op == "<") operators.push_back(2);
                else if (op == ">") operators.push_back(3);
                else if (op == "<=") operators.push_back(4);
                else if (op == ">=") operators.push_back(5);
                double restriction;
                if (file >> restriction) {
                    restrictions.push_back(restriction);
                }
            }
        }
        else {
            operators.push_back(0);
            restrictions.push_back(0);
        }
    }
    mate_copy = mate;
    operators_copy = operators;
    restrictions_copy = restrictions;
    int index = 1;
    int iszero = 0;
    while (index <= mate.cols) {
        file >> iszero;
        if (iszero == 1) {
            positive_coefs.push_back(index - 1);
        }
        index++;
    }
    positive_coefs_copy = positive_coefs;
    file >> taskType;
    taskType_copy = taskType;

    int i = 1;
    function.resize(1);
    while (file >> function[i - 1] && function.size() < mate.rows) {
        i++;
        function.resize(i);
    }
    function_copy = function;
    file >> function_constant;
    function_constant_copy = function_constant;
}

void SimplexEngine::printData(const Matrix& mate, const std::vector<int>& operators, const std::vector<double>& function, const std::vector<double>& restrictions, const std::vector<int>& positive_coefs, const double& function_constant) const {
    std::cout << "Restrictions matrix:\n";


    int width = 5; 

    for (int i = 0; i < mate.rows; i++) {
        for (int j = 0; j < mate.cols; j++) {
            std::cout << std::setw(width) << std::fixed << std::setprecision(2) << mate.data[i][j] << " ";
        }


        if (operators[i] == 1) std::cout << " = ";
        else if (operators[i] == 2) std::cout << " < ";
        else if (operators[i] == 3) std::cout << " > ";
        else if (operators[i] == 4) std::cout << " <= ";
        else if (operators[i] == 5) std::cout << " >= ";


        if (operators[i] != 0) std::cout << std::setw(width) << restrictions[i];

        std::cout << std::endl;
    }

    std::cout << "Target function:\n";
    for (double f : function) std::cout << std::setw(width) << std::fixed << std::setprecision(2) << f << " ";
    std::cout << " + constant: " << std::setw(width) << function_constant << std::endl;

    std::cout << "Positive coefficients:\n";
    for (int pos : positive_coefs) std::cout << std::setw(5) << pos + 1 << " ";
    std::cout << std::endl;

    std::cout << "Task type: " << taskType << std::endl;
    std::cout << "\n" << std::endl;
    

}

void SimplexEngine::Canonize(std::vector<double>& canon_function, std::vector<double>& canon_restrictions, Matrix& canon_matrix, std::vector<int>& canon_operators, std::vector<int>& canon_positive_coefs, std::vector<double>& canon_negative_coefs, std::vector<std::pair<int, int>>& canon_expanded_vars, std::vector<int>& canon_basis_ex, std::vector<double>& canon_leq_column, double& canon_function_constant, std::vector<int>& canon_leq, const std::string& canon_TaskType, const bool& isDouble, const bool& v) {
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;


    if (!isDouble) {
        if(v) std::cout << "Simplex Problem Canonization" << std::endl;
        if(v) std::cout << "Number of coefficients with no restrictions: " << r - canon_positive_coefs.size() << std::endl;
        for (int i = 0; i < canon_restrictions.size(); i++) {
            if (canon_restrictions[i] < 0) {
                canon_matrix.mulRow(i, -1);
                canon_restrictions[i] *= -1;
                if (canon_operators[i] == 4) {
                    canon_operators[i] = 5;
                }
                else if (canon_operators[i] == 5) {
                    canon_operators[i] = 4;
                }
                else {
                    continue;
                }
            }
        }
    }
    if (canon_matrix.rows - canon_positive_coefs.size() != 0) {
        if(v) std::cout << "Adding the necessary number of rows and copying them:\n";
        std::vector<double> canon_skipped_coefs;
        int checker = 0;
        if(v) std::cout << "Number of negative coefficients: ";
        canon_negative_coefs.resize(0);
        for (int pos : canon_positive_coefs) if(v) std::cout << pos + 1 << " ";
        if(v) std::cout << std::endl;
        for (int i = 0; i < canon_matrix.cols; i++) {
            checker = 0;
            for (int j = 0; j < canon_positive_coefs.size(); j++) {
                if (canon_positive_coefs[j] + 1 == i + 1) {
                    checker = 1;
                }
            }
            if (!checker) {
                canon_negative_coefs.push_back(i + 1);
            }
        }
        if(v) std::cout << "Number of coefficients with no restrictions: ";
        for (auto j : canon_negative_coefs) if(v) std::cout << j << " ";
        if(v) std::cout << std::endl;
        for (auto j : canon_negative_coefs) canon_expanded_vars.push_back({ j - 1, 0 });
        for (int i = 0; i < canon_negative_coefs.size(); i++) {
            canon_matrix.addColumn(1);

            canon_matrix.copyColumn(canon_negative_coefs[i] - 1, canon_matrix.cols - 1, -1);
            canon_function.push_back((-1) * canon_function[canon_negative_coefs[i] - 1]);
            canon_positive_coefs.push_back(canon_negative_coefs[i] - 1);
            canon_positive_coefs.push_back(canon_matrix.cols - 1);
            canon_expanded_vars[i].second = canon_matrix.cols - 1;
            std::sort(canon_positive_coefs.begin(), canon_positive_coefs.end());

        }
        if(v) printData(canon_matrix, canon_operators, canon_function, canon_restrictions, canon_positive_coefs, canon_function_constant);
    }
    if(v) std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::vector<int> canon_ineq_coefs;
    if (!isDouble) {
        for (int i = 0; i < canon_operators.size(); i++) {
            if (canon_operators[i] == 4) {
                canon_leq.push_back(1);
            }
            else {
                canon_leq.push_back(0);
            }
            if (canon_operators[i] == 4 || canon_operators[i] == 5) {
                canon_ineq_coefs.push_back(i);
            }
        }
    }
    else {
        for (int i = 0; i < canon_operators.size(); i++) {
            if (canon_operators[i] == 4 || canon_operators[i] == 5) {
                canon_ineq_coefs.push_back(i);
            }
        }
    }

    for (int i = 0; i < canon_ineq_coefs.size(); i++) {
        canon_matrix.addColumn(1);
        if (canon_operators[canon_ineq_coefs[i]] == 4) {
            canon_matrix.data[canon_ineq_coefs[i]][canon_matrix.cols - 1] = 1;
            canon_basis_ex.push_back(canon_matrix.cols - 1);
            canon_leq_column.push_back(canon_matrix.cols);
        }
        else if (canon_operators[canon_ineq_coefs[i]] == 5) {
            canon_matrix.data[canon_ineq_coefs[i]][canon_matrix.cols - 1] = -1;
        }
        canon_operators[canon_ineq_coefs[i]] = 1;
        canon_positive_coefs.push_back(canon_matrix.cols - 1);
        canon_function.push_back(0);
        if(v) printData(canon_matrix, canon_operators, canon_function, canon_restrictions, canon_positive_coefs, canon_function_constant);
    }

    if(isDouble) {
        if (!canon_TaskType.compare("max")) {
            for (int i = 0; i < canon_function.size(); i++) {
                canon_function[i] *= -1;
            }
        }
        for (int i = 0; i < canon_restrictions.size(); i++) {
            if (canon_restrictions[i] < 0) {
                canon_matrix.mulRow(i, -1);
                canon_restrictions[i] *= -1;
                if (canon_operators[i] == 4) canon_operators[i] = 5;
                else if (canon_operators[i] == 5) canon_operators[i] = 4;
                else continue;
            }
        }
    }
}

void SimplexEngine::artificialBasis(Matrix& artificial_mate, std::vector<int>& artificial_leq, std::vector<int>& artificial_basis_ex, std::vector<double>& artificial_function, std::vector<int>& artificial_artificial_basis, std::vector<int>& artificial_positive_coefs, const std::vector<int>& artificial_operators, const std::vector<double>& artificial_restrictions, const double& artificial_function_constant, const double& v) {
    int current_cols = artificial_mate.cols - 1;
    for (int i = 0; i < artificial_leq.size(); i++) {
        if (artificial_leq[i] == 0) {
            artificial_mate.addColumn(1);
            artificial_mate.data[i][artificial_mate.cols - 1] = 1;
            artificial_basis_ex.push_back(artificial_mate.cols - 1);
            artificial_artificial_basis.push_back(artificial_mate.cols - 1);
            artificial_function.push_back(N);
            artificial_positive_coefs.push_back(artificial_mate.cols - 1);
        }
    }
    if(v) std::cout << "==================== Construction of aritificial basis ====================" << std::endl;
    if(v) printData(artificial_mate, artificial_operators, artificial_function, artificial_restrictions, artificial_positive_coefs, artificial_function_constant);
    if(v) std::cout << "List of basis vectors: ";
    for (int i : artificial_basis_ex) {
        if(v) std::cout << i + 1 << " ";
    }
    if(v) std::cout << std::endl;
}