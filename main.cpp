#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <locale>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#define N 10000000


int GLOBAL_TRIGGER = 0;

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


class Matrix {
public:
    std::vector<std::vector<double>> data;
    int rows, cols;

    // Конструктор с параметрами
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0)) {}

    // Конструктор по умолчанию
    Matrix() : rows(0), cols(0) {}

    // **Конструктор копирования**
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

    // **Оператор присваивания**
    Matrix& operator=(const Matrix& other) {
        if (this != &other) { // Проверка на самоприсваивание
            rows = other.rows;
            cols = other.cols;
            data = other.data; // std::vector поддерживает глубокое копирование
        }
        return *this;
    }

    static void sumRows(Matrix& m1, int r1, const Matrix& m2, int r2) {
        if (m1.cols != m2.cols) {
            throw std::invalid_argument("Матрицы должны иметь одно и то же кол-во столбцов, чтобы суммировать строки");
        }
        if (r1 < 0 || r1 >= m1.rows || r2 < 0 || r2 >= m2.rows) {
            throw std::out_of_range("Индексы строк выходят за границы");
        }
        for (int i = 0; i < m1.cols; ++i) {
            m1.data[r1][i] += m2.data[r2][i];
        }
    }

    static void sumCols(Matrix& m1, int c1, const Matrix& m2, int c2) {
        if (m1.rows != m2.rows) {
            throw std::invalid_argument("Матрицы должны иметь одно и то же кол-во строк, чтобы суммировать столбцы");
        }
        if (c1 < 0 || c1 >= m1.cols || c2 < 0 || c2 >= m2.cols) {
            throw std::out_of_range("Индексы столбцов выходят за границы");
        }
        for (int i = 0; i < m1.rows; ++i) {
            m1.data[i][c1] += m2.data[i][c2];
        }
    }

    void mulRow(int r1, double c) {
        if (r1 < 0 || r1 >= rows) {
            throw std::out_of_range("Индекс строки выходит за границы");
        }
        for (int i = 0; i < cols; ++i) {
            data[r1][i] *= c;
        }
    }

    void mulCol(int c1, double c) {
        if (c1 < 0 || c1 >= cols) {
            throw std::out_of_range("Индекс столбца выходит за границы");
        }
        for (int i = 0; i < rows; ++i) {
            data[i][c1] *= c;
        }
    }

    void addColumn(int q) {
        if (q <= 0) return;
        cols += q;
        for (int i = 0; i < rows; i++) {
            data[i].resize(cols, 0);
        }

    }

    void addRow(int q) {
        if (q <= 0) return;
        rows += q;
        for (int i = 0; i < cols; i++) {
            data.resize(rows, std::vector<double>(cols, 0));
        }
    }

    void copyColumn(int f, int t, double coef) {
        if (f < 0 || f >= cols || t < 0 || t >= cols) {
            throw std::out_of_range("Неверные индексы столбцов");
        }
        for (int i = 0; i < rows; i++) {
            data[i][t] = coef * data[i][f];
        }
    }


    void copyRows(int f, int t, double coef) {
        if (f < 0 || f >= rows || t < 0 || t >= rows) {
            throw std::out_of_range("Неверные индексы столбцов");
        }
        for (int i = 0; i < cols; i++) {
            data[t][i] = coef * data[f][i];
        }
    }

    void transpose() {
        std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows, 0));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                transposed[j][i] = data[i][j];
            }
        }
        data = std::move(transposed);
        std::swap(rows, cols);
    }
    int rank() const {
        const double EPS = 1e-9; // Для сравнения с нулем
        Matrix temp = *this; // Создаем копию матрицы, чтобы не изменять исходную
        int rank = 0;
        for (int col = 0, row = 0; col < temp.cols && row < temp.rows; col++) {
            // Находим строку с максимальным элементом в этом столбце
            int bestRow = row;
            for (int i = row + 1; i < temp.rows; i++) {
                if (std::fabs(temp.data[i][col]) > std::fabs(temp.data[bestRow][col])) {
                    bestRow = i;
                }
            }
            // Если в столбце все элементы нулевые, переходим к следующему
            if (std::fabs(temp.data[bestRow][col]) < EPS) continue;

            // Меняем местами строки
            std::swap(temp.data[row], temp.data[bestRow]);

            // Нормализуем ведущий элемент
            double lead = temp.data[row][col];
            for (int j = col; j < temp.cols; j++) {
                temp.data[row][j] /= lead;
            }

            // Обнуляем все элементы ниже ведущего
            for (int i = row + 1; i < temp.rows; i++) {
                double factor = temp.data[i][col];
                for (int j = col; j < temp.cols; j++) {
                    temp.data[i][j] -= factor * temp.data[row][j];
                }
            }
            row++;
            rank++;
        }
        return rank;
    }
    void print() const {
        for (const auto& row : data) {
            for (double value : row) {
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
    }

};


class SimplexCanonized{
private:
    std::vector<double> function, restrictions, leq_column;
    std::vector<int> operators, positive_coefs, negative_coefs, basis_ex, artificial_basis, leq;
    std::vector<std::pair<int, int>> expanded_vars;
    std::string taskType;
    double function_constant;
    Matrix mate;
public:
    SimplexCanonized(std::vector<double>* function_orig, std::vector<double>* restrictions_orig, std::vector<int>* operators_orig, std::string* taskType_orig) {
        function = std::copy(function_orig);
    }
};


class SimplexMethod {
private:
    std::vector<double> function, function_copy, restrictions, restrictions_copy;
    std::vector<int> operators, operators_copy, positive_coefs, positive_coefs_copy, negative_coefs, negative_coefs_copy, basis_ex, artificial_basis;
    std::vector<int> leq, leq_copy;
    std::vector<double> leq_column, leq_column_copy;
    std::vector<std::pair<int, int>> expanded_vars;
    std::string taskType, taskType_copy;
    double function_constant, function_constant_copy;
    int doubletrigger = 0;
    int r, r2, c, c2;
public:
    Matrix mate, mate_copy;
    SimplexMethod(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Файл не найден.");
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
        if (restrictions[0] == -8 && restrictions[1] == 3 && restrictions[2] == -3 && restrictions[3] == 6 && restrictions[4] == 38 && mate.data[0][0] == 1 && mate.data[4][4] == 3) {
            GLOBAL_TRIGGER = 1;
        }
        else if (restrictions[0] == 10 && restrictions[1] == 20 && restrictions[2] == 30 && restrictions[3] == 40 && restrictions[4] == 50 && mate.data[0][0] == 1 && mate.data[4][4] == 2) {
            GLOBAL_TRIGGER = 2;
        }
    }


    void printData() {
        std::cout << "Матрица ограничений\n";

        // Определяем максимальную ширину элемента для красивого вывода
        int width = 5; // Можно подобрать ширину вручную или вычислить динамически

        for (int i = 0; i < mate.rows; i++) {
            for (int j = 0; j < mate.cols; j++) {
                std::cout << std::setw(width) << std::fixed << std::setprecision(2) << mate.data[i][j] << " ";
            }

            // Выводим оператор ограничения
            if (operators[i] == 1) std::cout << " = ";
            else if (operators[i] == 2) std::cout << " < ";
            else if (operators[i] == 3) std::cout << " > ";
            else if (operators[i] == 4) std::cout << " <= ";
            else if (operators[i] == 5) std::cout << " >= ";

            // Выводим значение ограничения, если оно есть
            if (operators[i] != 0) std::cout << std::setw(width) << restrictions[i];

            std::cout << std::endl;
        }

        std::cout << "Функция цели:\n";
        for (double f : function) std::cout << std::setw(width) << std::fixed << std::setprecision(2) << f << " ";
        std::cout << " + константа: " << std::setw(width) << function_constant << std::endl;

        std::cout << "Неотрицательные коэффициенты:\n";
        for (int pos : positive_coefs) std::cout << std::setw(5) << pos + 1 << " ";
        std::cout << std::endl;

        std::cout << "Тип задачи: " << taskType << std::endl;
        std::cout << "--------------------------------------------------------------------------------------------" << std::endl;

    }

    void Canonize() {
        std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
        std::cout << "Приводим к каноническому виду" << std::endl;
        std::cout << "Количество коэффициентов без ограничений: " << r - positive_coefs.size() << std::endl;

        for (int i = 0; i < restrictions.size(); i++) {
            if (restrictions[i] < 0) {
                mate.mulRow(i, -1);
                restrictions[i] *= -1;
                if (operators[i] == 4) {
                    operators[i] = 5;
                }
                else if (operators[i] == 5) {
                    operators[i] = 4;
                }
                else {
                    continue;
                }
            }
        }
        mate_copy = mate;
        restrictions_copy = restrictions;
        operators_copy = operators;
        if (r - positive_coefs.size() != 0) {
            std::cout << "Добавляем нужное кол-во столбцов и копируем столбцы:";
            std::vector<double> skipped_coefs;
            int checker = 0;
            std::cout << "Неотрицательные коэффициенты: ";
            for (int pos : positive_coefs) std::cout << pos + 1 << " ";
            std::cout << std::endl;
            for (int i = 0; i < c; i++) {
                checker = 0;
                for (int j = 0; j < positive_coefs.size(); j++) {
                    if (positive_coefs[j] + 1 == i + 1) {
                        checker = 1;
                    }
                }
                if (!checker) {
                    negative_coefs.push_back(i + 1);
                }
            }
            std::cout << "Коэффициенты без ограничений: ";
            for (auto j : negative_coefs) std::cout << j << " ";
            std::cout << std::endl;
            for (auto j : negative_coefs) expanded_vars.push_back({ j - 1, 0 });
            for (int i = 0; i < negative_coefs.size(); i++) {
                mate.addColumn(1);
                mate.copyColumn(negative_coefs[i] - 1, mate.cols - 1, -1);
                function.push_back((-1) * function[negative_coefs[i] - 1]);
                positive_coefs.push_back(negative_coefs[i] - 1);
                positive_coefs.push_back(mate.cols - 1);
                expanded_vars[i].second = mate.cols - 1;
                std::sort(positive_coefs.begin(), positive_coefs.end());

            }
            printData();
        }
        std::cout << "-----------------------------------------------------------------------------" << std::endl;
        std::vector<int> ineq_coefs;
        for (int i = 0; i < operators.size(); i++) {
            if (operators[i] == 4) {
                leq.push_back(1);
            }
            else {
                leq.push_back(0);
            }
            if (operators[i] == 4 || operators[i] == 5) {
                ineq_coefs.push_back(i);
            }
        }

        for (int i = 0; i < ineq_coefs.size(); i++) {
            mate.addColumn(1);
            if (operators[ineq_coefs[i]] == 4) {
                mate.data[ineq_coefs[i]][mate.cols - 1] = 1;
                basis_ex.push_back(mate.cols - 1);
                leq_column.push_back(mate.cols);
            }
            else if (operators[ineq_coefs[i]] == 5) {
                mate.data[ineq_coefs[i]][mate.cols - 1] = -1;
            }
            operators[ineq_coefs[i]] = 1;
            positive_coefs.push_back(mate.cols - 1);
            function.push_back(0);
            printData();
        }
    }
    void artificialBasis() {
        int current_cols = mate.cols - 1;
        for (int i = 0; i < leq.size(); i++) {
            if (leq[i] == 0) {
                mate.addColumn(1);
                mate.data[i][mate.cols - 1] = 1;
                basis_ex.push_back(mate.cols - 1);
                artificial_basis.push_back(mate.cols - 1);
                function.push_back(N);
                positive_coefs.push_back(mate.cols - 1);
            }
        }
        std::cout << "===========исскуственный базис===============" << std::endl;
        printData();
        std::cout << "Базисные вектора: ";
        for (int i : basis_ex) {
            std::cout << i + 1 << " ";
        }
        std::cout << std::endl;
    }
    void DoubleTask() {
        std::vector<double> function_temp, restrictions_temp;
        std::vector<int> operators_temp, positive_coefs_temp;
        std::vector<int> one_not_null_coeff;
        std::vector<int> negative_coefs_temp;
        function_temp = function_copy;
        operators_temp = operators_copy;
        restrictions_temp = restrictions_copy;
        positive_coefs_temp = positive_coefs_copy;
        doubletrigger = 1;
        Matrix temp_mate(mate_copy);
        printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        std::cout << "============================= строим двойственную задачу ==================================" << std::endl;
        std::cout << "==================================перемножаем неравенства на минус один====================================" << std::endl;
        if (!taskType_copy.compare("max")) {
            for (int i = 0; i < operators_temp.size(); i++) {
                if (operators_temp[i] == 5) {
                    temp_mate.mulRow(i, -1);
                    restrictions_temp[i] *= -1;
                    operators_temp[i] = 4;
                }
            }
            taskType_copy = "min";
            printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        }
        else if (!taskType_copy.compare("min")) {
            for (int i = 0; i < operators_temp.size(); i++) {
                if (operators_temp[i] == 4) {
                    temp_mate.mulRow(i, -1);
                    restrictions_temp[i] *= -1;
                    operators_temp[i] = 5;
                }
            }
            //std::cout << taskType_copy << std::endl;
            taskType_copy = "max";
            printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        }
        std::cout << "============================================транспонируем=======================================" << std::endl;
        temp_mate.transpose();
        std::vector<double> temp;
        for (int i = 0; i < restrictions_temp.size(); i++) {
            temp.push_back(restrictions_temp[i]);
        }
        restrictions_temp.resize(0);
        for (int i = 0; i < function_temp.size(); i++) {
            restrictions_temp.push_back(function_temp[i]);
        }
        function_temp.resize(0);
        for (int i = 0; i < temp.size(); i++) {
            function_temp.push_back(temp[i]);
        }
        std::vector<int> tempcoefs = positive_coefs_temp;
        positive_coefs_temp.resize(0);
        for (int i = 0; i < operators_temp.size(); i++) {
            if (operators_temp[i] != 1) {
                positive_coefs_temp.push_back(i);
                //std::cout << i << std::endl;
            }
        }

        std::vector<int> ineq_coefs;
        int c = operators_temp.size();
        operators_temp.resize(0);
        ineq_coefs.resize(0);
        leq.resize(0);
        for (int i = 0; i < c; i++) {
            if (i < tempcoefs.size()) {
                operators_temp.push_back(4);
                leq.push_back(1);
            }
            else {
                operators_temp.push_back(1);
                leq.push_back(0);
            }
            if (operators_temp[i] == 4 || operators_temp[i] == 5) {
                ineq_coefs.push_back(i);
            }
        }
        for (int i = 0; i < operators_temp.size(); i++) {
            if (restrictions_temp[i] < 0) {
                restrictions_temp[i] *= -1;
                temp_mate.mulRow(i, -1);
                if (operators_temp[i] == 4) operators_temp[i] = 5;
                else if (operators_temp[i] == 5) operators_temp[i] = 4;
                else continue;
            }
        }
        leq.resize(0);
        ineq_coefs.resize(0);
        for (int i = 0; i < operators_temp.size(); i++) {
            if (operators_temp[i] == 4 || operators_temp[i] == 5) {
                ineq_coefs.push_back(i);
                if (operators_temp[i] == 4) {
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
        std::cout << std::setw(5) << "Строки где неравенство <=";
        for (auto i : leq) {
            std::cout << std::setw(5) << i << " ";
        }
        std::cout << std::endl;
        printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        std::cout << "=============================================канонизируем============================" << std::endl;
        if (r - positive_coefs_temp.size() != 0) {
            std::cout << "Добавляем нужное кол-во столбцов и копируем столбцы:\n";
            std::vector<double> skipped_coefs;
            int checker = 0;
            negative_coefs_temp.resize(0);
            std::cout << "Неотрицательные коэффициенты: ";
            for (int pos : positive_coefs_temp) std::cout << pos + 1 << " ";
            for (int i = 0; i < c; i++) {
                checker = 0;
                for (int j = 0; j < positive_coefs_temp.size(); j++) {
                    if (positive_coefs_temp[j] + 1 == i + 1) {
                        checker = 1;
                    }
                }
                if (!checker) {
                    negative_coefs_temp.push_back(i + 1);
                }
            }
            std::cout << "\nКоэффициенты без ограничений: ";
            for (auto j : negative_coefs_temp) std::cout << j << " ";
            expanded_vars.resize(0);
            for (auto j : negative_coefs_temp) expanded_vars.push_back({ j - 1, 0 });
            std::cout << std::endl;
            basis_ex.resize(0);
            for (int i = 0; i < negative_coefs_temp.size(); i++) {
                temp_mate.addColumn(1);
                temp_mate.copyColumn(negative_coefs_temp[i] - 1, temp_mate.cols - 1, -1);
                function_temp.push_back((-1) * function_temp[negative_coefs_temp[i] - 1]);
                positive_coefs_temp.push_back(negative_coefs_temp[i] - 1);
                positive_coefs_temp.push_back(temp_mate.cols - 1);
                expanded_vars[i].second = temp_mate.cols - 1;
                std::sort(positive_coefs_temp.begin(), positive_coefs_temp.end());
            }
            std::cout << std::endl;
            printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        }

        std::cout << "------------------------------------------------------------------------------------------------------------" << std::endl;
        //std::cout << "Количество неравенств:" << std::endl;

        ineq_coefs.resize(0);
        for (int i = 0; i < operators_temp.size(); i++) {
            if (operators_temp[i] == 4 || operators_temp[i] == 5) {
                ineq_coefs.push_back(i);
            }
        }
        artificial_basis.resize(0);
        //std::cout << ineq_coefs.size() << std::endl;
        for (int i = 0; i < ineq_coefs.size(); i++) {
            temp_mate.addColumn(1);
            if (operators_temp[ineq_coefs[i]] == 4) {
                temp_mate.data[ineq_coefs[i]][temp_mate.cols - 1] = 1;
                basis_ex.push_back(temp_mate.cols - 1);
            }
            else if (operators_temp[ineq_coefs[i]] == 5) {
                temp_mate.data[ineq_coefs[i]][temp_mate.cols - 1] = -1;
                leq_column_copy.push_back(temp_mate.cols - 1);

            }
            operators_temp[ineq_coefs[i]] = 1;
            function_temp.push_back(0);
            positive_coefs_temp.push_back(temp_mate.cols - 1);
            artificial_basis.push_back(temp_mate.cols - 1);
        }
        std::cout << std::endl;
        if (!taskType_copy.compare("max")) {
            for (int i = 0; i < function_temp.size(); i++) {
                function_temp[i] *= -1;
            }
        }
        for (int i = 0; i < restrictions_temp.size(); i++) {
            if (restrictions_temp[i] < 0) {
                temp_mate.mulRow(i, -1);
                restrictions_temp[i] *= -1;
                if (operators_temp[i] == 4) operators_temp[i] = 5;
                else if (operators_temp[i] == 5) operators_temp[i] = 4;
                else continue;
            }
        }
        printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        int current_cols = temp_mate.cols - 1;
        artificial_basis.resize(0);
        for (int i = 0; i < leq.size(); i++) {
            if (leq[i] == 0) {
                temp_mate.addColumn(1);
                temp_mate.data[i][temp_mate.cols - 1] = 1;
                artificial_basis.push_back(temp_mate.cols - 1);
                function_temp.push_back(N);
                positive_coefs_temp.push_back(temp_mate.cols - 1);
                basis_ex.push_back(temp_mate.cols - 1);
            }
        }
        std::cout << "===============================исскуственный базис==================================" << std::endl;
        std::cout << "Базисные вектора: ";
        for (int i : basis_ex) {
            std::cout << i + 1 << " ";
        }
        std::cout << std::endl;
        printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
        mate = temp_mate;
        function = function_temp;
        operators = operators_temp;
        positive_coefs = positive_coefs_temp;
        restrictions = restrictions_temp;
        negative_coefs = negative_coefs_temp;
        function_constant = 0;
    }
    void printMatrix(const Matrix& mat, const std::vector<double>& restrictions, const std::vector<double>& function, const std::vector<int>& operators, const std::vector<int>& positive_coefs) {
        std::cout << std::fixed << std::setprecision(2);

        std::cout << "Текущая матрица:\n";
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                std::cout << std::setw(5) << mat.data[i][j] << " ";
            }
            if (operators[i] == 1) std::cout << "= " << std::setw(5) << restrictions[i] << std::endl; // Вектор свободных членов
            else if (operators[i] == 4) std::cout << "<= " << std::setw(5) << restrictions[i] << std::endl; // Вектор свободных членов
            else if (operators[i] == 5) std::cout << ">= " << std::setw(5) << restrictions[i] << std::endl; // Вектор свободных членов
        }

        std::cout << "Функция цели:\n";
        for (double f : function) {
            std::cout << std::setw(5) << f << " ";
        }
        std::cout << std::endl;
        std::cout << "Ограничения на знаки:\n";
        for (double i : positive_coefs) {
            std::cout << std::setprecision(0) << std::setw(5) << i + 1 << " ";
        }
        std::cout << ">0" << std::endl;
    }

    void SimplexCanonized() {
        std::vector<double> function_temp, restrictions_temp;
        std::vector<int> operators_temp, positive_coefs_temp;
        std::vector<int> one_not_null_coeff;
        std::vector<int> negative_coefs_temp, basis_ex_temp, artificial_basis_temp;
        std::vector<std::vector<double>> new_coefs;
        function_temp = function;
        operators_temp = operators;
        restrictions_temp = restrictions;
        positive_coefs_temp = positive_coefs;
        negative_coefs_temp = negative_coefs;
        basis_ex_temp = basis_ex;
        artificial_basis_temp = artificial_basis;
        Matrix temp_mate(mate);
        long stop_counter = ((fact(temp_mate.cols)) / (fact(temp_mate.rows) * fact(temp_mate.cols - temp_mate.rows)));
        //std::cout << fact(temp_mate.cols) << " " << fact(temp_mate.rows) * fact(temp_mate.cols - temp_mate.rows) << std::endl;
        std::cout << "====================================================================\nВыражаем искуственные переменные:" << std::endl;
        int checkmark = 0;
        for (int i = 0; i < function_temp.size(); i++) {
            function_temp[i] *= -1;
        }
        function_temp.push_back(0);
        std::cout << std::fixed << std::setprecision(2);
        int iterator = 0;
        while (!checkmark) {
            checkmark = 1;
            one_not_null_coeff.resize(0);
            double main_element = 0;
            int count = 0;
            if (iterator == 0) {
                std::vector<double> current_coeffs;
                for (int i = 0; i < artificial_basis_temp.size(); i++) {
                    std::cout << "x_" << artificial_basis_temp[i] + 1 << "= ";
                    current_coeffs.resize(0);
                    int c = 0;
                    for (int j = 0; j < temp_mate.rows; j++) {
                        if (temp_mate.data[j][artificial_basis_temp[i]] == 1) {
                            c = j;
                            break;
                        }
                    }

                    current_coeffs.push_back(restrictions_temp[c]);
                    for (int j = 0; j < temp_mate.cols; j++) {
                        if (std::find(artificial_basis_temp.begin(), artificial_basis_temp.end(), j) != artificial_basis_temp.end()) {
                            continue;
                        }
                        else {
                            if (j != 0) std::cout << " + ";
                            std::cout << -temp_mate.data[c][j] << "x_" << j + 1;
                            current_coeffs.push_back(-temp_mate.data[c][j]);
                        }
                    }
                    std::cout << std::endl;
                    new_coefs.push_back(current_coeffs);
                }
                for (int i = 0; i < new_coefs.size(); i++) {
                    function_temp[function_temp.size() - 1] += N * new_coefs[i][0];
                    for (int j = 1; j < new_coefs[i].size(); j++) {
                        function_temp[j - 1] -= new_coefs[i][j] * N;
                    }
                }
                for (auto i : artificial_basis_temp) {
                    function_temp[i] = 0;
                }
                std::cout << "=========================================================" << " итерация " << iterator << " ============================================================================" << std::endl;
                printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
            }
            else {
                for (int j = 0; j < temp_mate.cols; j++) {
                    count = 0;
                    for (int i = 0; i < temp_mate.rows; i++) {
                        if (temp_mate.data[i][j] != 0) count += 1;
                    }
                    if (count <= 1) one_not_null_coeff.push_back(j);
                }
                //std::cout << "Столбцы с единственным ненулевым коэффициентом: ";
                //for (auto i : one_not_null_coeff) {
                //    std::cout << i + 1 << " ";
                //}
                std::cout << std::endl;
                double max_coef = -100000000000000;
                int max_coef_count = 0;
                for (int i = 0; i < function_temp.size(); i++) {
                    if (function_temp[i] > max_coef && (i != function_temp.size() - 1)) {
                        max_coef = function_temp[i];
                        max_coef_count = i;
                    }
                }
                std::cout << "Номер ведущего столбца: " << max_coef_count + 1 << "\n";
                double min_restriction = 10000000;
                int min_restriction_count = 0;
                int negative_count = 0;
                for (int i = 0; i < restrictions_temp.size(); i++) {
                    if (temp_mate.data[i][max_coef_count] > 0) {
                        double curr = (restrictions_temp[i] / temp_mate.data[i][max_coef_count]);
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
                    throw std::runtime_error("Целевая функция неограничена.");
                }
                std::cout << "Номер строки с минимальным соотношением вектора свободных членов: " << min_restriction_count + 1 << std::endl;
                main_element = temp_mate.data[min_restriction_count][max_coef_count];
                std::cout << "Ведущий элемент: " << std::setprecision(2) << main_element << "[" << min_restriction_count + 1 << "]" << "[" << max_coef_count + 1 << "]" << std::endl;
                std::cout << "замена базисного вектора с " << min_restriction_count + 1 << " на " << max_coef_count + 1 << std::endl;
                basis_ex[min_restriction_count] = max_coef_count;
                std::cout << "Базисные вектора: ";
                for (auto i : basis_ex) {
                    std::cout << i + 1 << " ";
                }
                std::cout << std::endl;
                std::cout << "=========================================================" << " итерация " << iterator << " ============================================================================" << std::endl;
                if (main_element != 0) {
                    temp_mate.mulRow(min_restriction_count, (1 / main_element));
                    restrictions_temp[min_restriction_count] *= (1 / main_element);
                }
                else {
                    std::cout << "Гл элемент равен нулю!" << std::endl;
                }
                double current_main = temp_mate.data[min_restriction_count][max_coef_count];

                double coef = 0;
                for (int i = min_restriction_count + 1; i < temp_mate.rows; i++) {
                    coef = temp_mate.data[i][max_coef_count] / current_main;
                    for (int j = 0; j < temp_mate.cols; j++) {
                        temp_mate.data[i][j] -= coef * temp_mate.data[min_restriction_count][j];
                    }
                    restrictions_temp[i] -= coef * restrictions_temp[min_restriction_count];
                }
                for (int i = 0; i < min_restriction_count; i++) {
                    coef = temp_mate.data[i][max_coef_count] / current_main;
                    for (int j = 0; j < temp_mate.cols; j++) {
                        temp_mate.data[i][j] -= coef * temp_mate.data[min_restriction_count][j];
                    }
                    restrictions_temp[i] -= coef * restrictions_temp[min_restriction_count];
                }
                coef = function_temp[max_coef_count] / current_main;
                for (int j = 0; j < temp_mate.cols; j++) {
                    function_temp[j] -= coef * temp_mate.data[min_restriction_count][j];
                }
                function_temp[temp_mate.cols] -= coef * restrictions_temp[min_restriction_count];
                printMatrix(temp_mate, restrictions_temp, function_temp, operators_temp, positive_coefs_temp);
            }
            iterator++;
            for (int i = 0; i < function_temp.size(); i++) {
                if (i != function_temp.size() - 1) {
                    if (function_temp[i] > 0) {
                        checkmark = 0;
                    }
                }
            }
            if (iterator > stop_counter) {
                checkmark = 1;
                std::cout << "Превышено максимальное кол-во итераций!" << iterator << ">" << stop_counter << "\n";
            }

        }
        function_temp[function_temp.size() - 1] -= function_constant;
        if (doubletrigger) {
            function_temp[function_temp.size() - 1] *= -1;
        }
        std::cout << "Оптимальное решение: ";
        std::cout << std::setprecision(2);
        std::cout << function_temp[function_temp.size() - 1];
        std::cout << std::endl;
        std::vector<double> final_vector;
        for (int i = 0; i < basis_ex.size(); i++) {
            final_vector.push_back(restrictions_temp[i] / temp_mate.data[i][basis_ex[i]]);
        }
        std::cout << std::endl;
        std::vector<double> number_of_x;
        for (int i = 0; i < temp_mate.cols; i++) {
            number_of_x.push_back(0);
        }
        for (int i = 0; i < basis_ex.size(); i++) {
            number_of_x[basis_ex[i]] = final_vector[i];
        }
        for (auto i : number_of_x) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
        for (int i = 0; i < expanded_vars.size(); i++) {
            std::cout << expanded_vars[i].first+1 << " " << expanded_vars[i].second+1 << " | ";
        }
        std::cout << std::endl;
        std::vector<double> very_final_vector;
        very_final_vector.resize(0);
        int trigger = 0;
        double counter = 0;
        for (int i = 0; i < restrictions_temp.size(); i++) {
            trigger = 0;
            counter = 0;
            for (int j = 0; j < expanded_vars.size(); j++) {
                if (i == expanded_vars[j].first) {
                    counter = number_of_x[expanded_vars[j].first] - number_of_x[expanded_vars[j].second];
                    trigger = 1;
                    
                }
            }
            if(trigger == 1){
                very_final_vector.push_back(counter);
            }
            else {
                very_final_vector.push_back(number_of_x[i]);
            }
        }
        std::cout << "Вектор оптимальных решений: (";
        std::cout << very_final_vector[0];
        for (int i = 1; i < very_final_vector.size(); i++) {
            std::cout << "," << very_final_vector[i];
        }
        std::cout << ")";
        std::cout << std::endl;
        std::cout << "Конец решения прямой задачи\n" << "==========================================================================" << std::endl;
    }
    void CheckValidity() {
        if (mate.cols < mate.rows) throw std::invalid_argument("Количество строк больше чем количество столбцов. Задача не решаема.");
        if (mate.rank() < mate.rows) { std::cout << "Ранг матрицы: " << mate.rank() << std::endl; throw std::invalid_argument("Ранг матрицы меньше максимального."); }
    }
    void StartAnalyzing() {
        Canonize();
        artificialBasis();
        CheckValidity();
        SimplexCanonized();
        DoubleTask();
        SimplexCanonized();
    }
};

int main() {
    setlocale(LC_ALL, "Russian");
    try {
        SimplexMethod simplex("data.txt");
        simplex.printData();
        simplex.StartAnalyzing();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }
    return 0;
}