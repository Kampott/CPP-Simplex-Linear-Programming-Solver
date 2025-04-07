#pragma once
#include <vector>
#include <stdexcept>
#include <iostream>
class Matrix {
public:
    std::vector<std::vector<double>> data;
    int rows, cols;
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0)) {}
    Matrix() : rows(0), cols(0) {}
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

    Matrix& operator=(const Matrix& other) {
        if (this != &other) { // ѕроверка на самоприсваивание
            rows = other.rows;
            cols = other.cols;
            data = other.data; // std::vector поддерживает глубокое копирование
        }
        return *this;
    }

    void sumRows(Matrix& m1, int r1, const Matrix& m2, int r2);

    static void sumCols(Matrix& m1, int c1, const Matrix& m2, int c2);

    void mulRow(int r1, double c);
 

    void mulCol(int c1, double c);

    void addColumn(int q);

    void addRow(int q);

    void copyColumn(int f, int t, double coef);


    void copyRows(int f, int t, double coef);

    void transpose();
    int rank() const;
    void print() const;

};
