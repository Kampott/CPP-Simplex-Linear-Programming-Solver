#include "SimplexSolver.h"




SimplexSolver::SimplexSolver(const std::vector<double>& function_orig, const std::vector<double>& restrictions_orig, const std::vector<int>& operators_orig, const std::string& taskType_orig, const Matrix& mate_orig, SimplexEngine* engine, const std::vector<int>& positive_coefs_orig, const double& function_constant_orig, const bool& isDouble_orig, const bool& verbose_orig) : function(function_orig), restrictions(restrictions_orig), operators(operators_orig), taskType(taskType_orig), mate(mate_orig), method(engine), positive_coefs(positive_coefs_orig), function_constant(function_constant_orig), isDouble(isDouble_orig), v(verbose_orig) {};

void SimplexSolver::SolveDouble() {

    if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
    if(v) std::cout << "=================================== Double Task Preparation ===================================" << std::endl;


    if(v) std::cout << "Changin the inequality sign, by multiplying the row by -1" << std::endl;
    if (!taskType.compare("max")) {
        for (int i = 0; i < operators.size(); i++) {
            if (operators[i] == 5) {
                mate.mulRow(i, -1);
                restrictions[i] *= -1;
                operators[i] = 4;
            }
        }
        taskType = "min";
        if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
    }
    else if (!taskType.compare("min")) {
        for (int i = 0; i < operators.size(); i++) {
            if (operators[i] == 4) {
                mate.mulRow(i, -1);
                restrictions[i] *= -1;
                operators[i] = 5;
            }
        }
        taskType = "max";
        if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
    }


    if(v) std::cout << "\n\nTransponing the restrictions matrix" << std::endl;
    mate.transpose();
    std::vector<double> temp;
    for (int i = 0; i < restrictions.size(); i++) {
        temp.push_back(restrictions[i]);
    }
    restrictions.resize(0);
    for (int i = 0; i < function.size(); i++) {
        restrictions.push_back(function[i]);
    }
    function.resize(0);
    for (int i = 0; i < temp.size(); i++) {
        function.push_back(temp[i]);
    }
    std::vector<int> tempcoefs = positive_coefs;
    positive_coefs.resize(0);
    for (int i = 0; i < operators.size(); i++) {
        if (operators[i] != 1) {
            positive_coefs.push_back(i);
        }
    }

    std::vector<int> ineq_coefs;
    int c = operators.size();
    operators.resize(0);
    ineq_coefs.resize(0);
    leq.resize(0);
    for (int i = 0; i < c; i++) {
        if (i < tempcoefs.size()) {
            operators.push_back(4);
            leq.push_back(1);
        }
        else {
            operators.push_back(1);
            leq.push_back(0);
        }
        if (operators[i] == 4 || operators[i] == 5) {
            ineq_coefs.push_back(i);
        }
    }
    for (int i = 0; i < operators.size(); i++) {
        if (restrictions[i] < 0) {
            restrictions[i] *= -1;
            mate.mulRow(i, -1);
            if (operators[i] == 4) operators[i] = 5;
            else if (operators[i] == 5) operators[i] = 4;
            else continue;
        }
    }
    leq.resize(0);
    ineq_coefs.resize(0);
    for (int i = 0; i < operators.size(); i++) {
        if (operators[i] == 4 || operators[i] == 5) {
            ineq_coefs.push_back(i);
            if (operators[i] == 4) {
                leq.push_back(1);
            }
            else {
                leq.push_back(0);
            }
        }
        else {
            leq.push_back(0);
        }
    }
    if(v) std::cout << std::setw(5) << "Columns where the inequality is of <= type";
    for (auto i : leq) {
        if(v) std::cout << std::setw(5) << i << " ";
    }
    if(v) std::cout << std::endl;
    if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
}
void SimplexSolver::Solve() {
    if (isDouble) {
        SolveDouble();
    }
    method->Canonize(function, restrictions, mate, operators, positive_coefs, negative_coefs, expanded_vars, basis_ex, leq_column, function_constant, leq, taskType, isDouble, v);
    method->artificialBasis(mate, leq, basis_ex, function, artificial_basis, positive_coefs, operators, restrictions, function_constant, v);
    std::vector<int> one_not_null_coeff;
    std::vector<std::vector<double>> new_coefs;
    long stop_counter = ((method->fact(mate.cols)) / (method->fact(mate.rows) * method->fact(mate.cols - mate.rows)));
    if(v) std::cout << "====================================================================\nExpressing the artificial variables:" << std::endl;
    int checkmark = 0;
    for (int i = 0; i < function.size(); i++) {
        function[i] *= -1;
    }
    function.push_back(0);
    if(v) std::cout << std::fixed << std::setprecision(2);
    int iterator = 0;
    while (!checkmark) {
        checkmark = 1;
        one_not_null_coeff.resize(0);
        double main_element = 0;
        int count = 0;
        if (iterator == 0) {
            std::vector<double> current_coeffs;
            for (int i = 0; i < artificial_basis.size(); i++) {
                if(v) std::cout << "x_" << artificial_basis[i] + 1 << "= ";
                current_coeffs.resize(0);
                int c = 0;
                for (int j = 0; j < mate.rows; j++) {
                    if (mate.data[j][artificial_basis[i]] == 1) {
                        c = j;
                        break;
                    }
                }

                current_coeffs.push_back(restrictions[c]);
                for (int j = 0; j < mate.cols; j++) {
                    if (std::find(artificial_basis.begin(), artificial_basis.end(), j) != artificial_basis.end()) {
                        continue;
                    }
                    else {
                        if (j != 0) if(v) std::cout << " + ";
                        if(v) std::cout << -mate.data[c][j] << "x_" << j + 1;
                        current_coeffs.push_back(-mate.data[c][j]);
                    }
                }
                if(v) std::cout << std::endl;
                new_coefs.push_back(current_coeffs);
            }
            for (int i = 0; i < new_coefs.size(); i++) {
                function[function.size() - 1] += N * new_coefs[i][0];
                for (int j = 1; j < new_coefs[i].size(); j++) {
                    function[j - 1] -= new_coefs[i][j] * N;
                }
            }
            for (auto i : artificial_basis) {
                function[i] = 0;
            }
            if(v) std::cout << "=========================================================" << " iteration " << iterator+1 << " ============================================================================" << std::endl;
            if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
        }
        else {
            for (int j = 0; j < mate.cols; j++) {
                count = 0;
                for (int i = 0; i < mate.rows; i++) {
                    if (mate.data[i][j] != 0) count += 1;
                }
                if (count <= 1) one_not_null_coeff.push_back(j);
            }
            if(v) std::cout << std::endl;
            double max_coef = -100000000000000;
            int max_coef_count = 0;
            for (int i = 0; i < function.size(); i++) {
                if (function[i] > max_coef && (i != function.size() - 1)) {
                    max_coef = function[i];
                    max_coef_count = i;
                }
            }
            if(v) std::cout << "The number of the leading row: " << max_coef_count + 1 << "\n";
            double min_restriction = 10000000;
            int min_restriction_count = 0;
            int negative_count = 0;
            for (int i = 0; i < restrictions.size(); i++) {
                if (mate.data[i][max_coef_count] > 0) {
                    double curr = (restrictions[i] / mate.data[i][max_coef_count]);
                    if (curr >= 0) {
                        if (curr < min_restriction) {
                            min_restriction = curr;
                            min_restriction_count = i;
                        }
                    }
                    else if (curr < 0) {
                        negative_count++;
                    }
                }
                else {
                    negative_count++;
                }
            }
            if (negative_count == restrictions.size()) {
                throw std::runtime_error("The target function is unlimited.");
            }
            if(v) std::cout << "The number of the column with the least coefficient of division by the vector of free members: " << min_restriction_count + 1 << std::endl;
            main_element = mate.data[min_restriction_count][max_coef_count];
            if(v) std::cout << "Leading element: " << std::setprecision(2) << main_element << "[" << min_restriction_count + 1 << "]" << "[" << max_coef_count + 1 << "]" << std::endl;
            if(v) std::cout << "replacing the basis vector form " << min_restriction_count + 1 << " to " << max_coef_count + 1 << std::endl;
            if(v) std::cout << basis_ex.size() << " " << min_restriction_count << std::endl;
            basis_ex[min_restriction_count] = max_coef_count;
            if(v) std::cout << "Basis vectors: ";
            for (auto i : basis_ex) {
                if(v) std::cout << i + 1 << " ";
            }
            if(v) std::cout << std::endl;
            if(v) std::cout << "=========================================================" << " Iteration " << iterator+1 << " ============================================================================" << std::endl;
            if (main_element != 0) {
                mate.mulRow(min_restriction_count, (1 / main_element));
                restrictions[min_restriction_count] *= (1 / main_element);
            }
            else {
                if(v) std::cout << "Main element is equal to 0!" << std::endl;
            }
            double current_main = mate.data[min_restriction_count][max_coef_count];

            double coef = 0;
            for (int i = min_restriction_count + 1; i < mate.rows; i++) {
                coef = mate.data[i][max_coef_count] / current_main;
                for (int j = 0; j < mate.cols; j++) {
                    mate.data[i][j] -= coef * mate.data[min_restriction_count][j];
                }
                restrictions[i] -= coef * restrictions[min_restriction_count];
            }
            for (int i = 0; i < min_restriction_count; i++) {
                coef = mate.data[i][max_coef_count] / current_main;
                for (int j = 0; j < mate.cols; j++) {
                    mate.data[i][j] -= coef * mate.data[min_restriction_count][j];
                }
                restrictions[i] -= coef * restrictions[min_restriction_count];
            }
            coef = function[max_coef_count] / current_main;
            for (int j = 0; j < mate.cols; j++) {
                function[j] -= coef * mate.data[min_restriction_count][j];
            }
            function[mate.cols] -= coef * restrictions[min_restriction_count];
            if(v) method->printMatrix(mate, restrictions, function, operators, positive_coefs);
        }
        iterator++;
        for (int i = 0; i < function.size(); i++) {
            if (i != function.size() - 1) {
                if (function[i] > 0) {
                    checkmark = 0;
                }
            }
        }
        if (iterator > stop_counter) {
            checkmark = 1;
            if(v) std::cout << "Exceeded maximum number of iterations!" << iterator << ">" << stop_counter << "\n";
        }

    }
    function[function.size() - 1] -= function_constant;
    if(v) std::cout << "Optimal result: ";
    if(v) std::cout << std::setprecision(2);
    if(v) std::cout << function[function.size() - 1];
    if(v) std::cout << std::endl;
    std::vector<double> final_vector;
    for (int i = 0; i < basis_ex.size(); i++) {
        final_vector.push_back(restrictions[i] / mate.data[i][basis_ex[i]]);
    }
    if(v) std::cout << std::endl;
    std::vector<double> number_of_x;
    for (int i = 0; i < mate.cols; i++) {
        number_of_x.push_back(0);
    }
    for (int i = 0; i < basis_ex.size(); i++) {
        number_of_x[basis_ex[i]] = final_vector[i];
    }
    for (auto i : number_of_x) {
        if(v) std::cout << i << " ";
    }
    if(v) std::cout << std::endl;
    for (int i = 0; i < expanded_vars.size(); i++) {
        if(v) std::cout << expanded_vars[i].first + 1 << " " << expanded_vars[i].second + 1 << " | ";
    }
    if(v) std::cout << std::endl;
    std::vector<double> very_final_vector;
    very_final_vector.resize(0);
    int trigger = 0;
    double counter = 0;
    for (int i = 0; i < restrictions.size(); i++) {
        trigger = 0;
        counter = 0;
        for (int j = 0; j < expanded_vars.size(); j++) {
            if (i == expanded_vars[j].first) {
                counter = number_of_x[expanded_vars[j].first] - number_of_x[expanded_vars[j].second];
                trigger = 1;

            }
        }
        if (trigger == 1) {
            very_final_vector.push_back(counter);
        }
        else {
            very_final_vector.push_back(number_of_x[i]);
        }
    }
    std::cout << "Vector of optimal results: (";
    std::cout << very_final_vector[0];
    for (int i = 1; i < very_final_vector.size(); i++) {
        std::cout << "," << very_final_vector[i];
    }
    std::cout << ")";
    std::cout << std::endl;
    if (!isDouble) {
        std::cout << "Finished solving canonical LP task\n" << "==========================================================================\n\n\n" << std::endl;
    }
    else {
        std::cout << "Finished solving double LP task\n" << "=========================================================================\n\n\n" << std::endl;
    }
}