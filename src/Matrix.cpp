#include "Matrix.h"
#include <iomanip>

void Matrix::sumRows(Matrix& m1, int r1, const Matrix& m2, int r2) {
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

void Matrix::sumCols(Matrix& m1, int c1, const Matrix& m2, int c2) {
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

void Matrix::mulRow(int r1, double c) {
    if (r1 < 0 || r1 >= rows) {
        throw std::out_of_range("Индекс строки выходит за границы");
    }
    for (int i = 0; i < cols; ++i) {
        data[r1][i] *= c;
    }
}

void Matrix::mulCol(int c1, double c) {
    if (c1 < 0 || c1 >= cols) {
        throw std::out_of_range("Индекс столбца выходит за границы");
    }
    for (int i = 0; i < rows; ++i) {
        data[i][c1] *= c;
    }
}
void Matrix::addColumn(int q) {
    if (q <= 0) return;
    cols += q;
    for (int i = 0; i < rows; i++) {
        data[i].resize(cols, 0);
    }

}

void Matrix::addRow(int q) {
    if (q <= 0) return;
    rows += q;
    for (int i = 0; i < cols; i++) {
        data.resize(rows, std::vector<double>(cols, 0));
    }
}

void Matrix::copyColumn(int f, int t, double coef) {
    if (f < 0 || f >= cols || t < 0 || t >= cols) {
        throw std::out_of_range("Неверные индексы столбцов");
    }
    for (int i = 0; i < rows; i++) {
        data[i][t] = coef * data[i][f];
    }
}


void Matrix::copyRows(int f, int t, double coef) {
    if (f < 0 || f >= rows || t < 0 || t >= rows) {
        throw std::out_of_range("Неверные индексы столбцов");
    }
    for (int i = 0; i < cols; i++) {
        data[t][i] = coef * data[f][i];
    }
}

void Matrix::transpose() {
    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = data[i][j];
        }
    }
    data = std::move(transposed);
    std::swap(rows, cols);
}
int Matrix::rank() const {
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
void Matrix::print() const {
    for (const auto& row : data) {
        for (double value : row) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}