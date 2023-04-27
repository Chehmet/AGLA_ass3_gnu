#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cmath>

#define GNUPLOT_NAME "gnuplot -persist"

using namespace std;

// class Matrix is for storing matrix and making actions to it
class Matrix {
public:
    int rows, columns;
    vector<vector<double>> data;

    // input (>>) overloading
    friend istream &operator>>(istream &in, Matrix &inputMatrix) {
        vector<vector<double>> values;
        for (int i = 0; i < inputMatrix.rows; i++) {
            vector<double> temp;
            for (int j = 0; j < inputMatrix.columns; j++) {
                double element;
                in >> element;
                temp.push_back(element);
            }
            values.push_back(temp);
        }
        inputMatrix.data = values;
        return in;
    }

    // output (<<) overloading
    friend ostream &operator<<(ostream &out, const Matrix &outputMatrix) {
        for (int i = 0; i < outputMatrix.rows; i++) {
            for (int j = 0; j < outputMatrix.columns; j++) {
                out << fixed << setprecision(4) << outputMatrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }

    Matrix() {

    }

    Matrix(int rows, int cols) {
        this->rows = rows;
        this->columns = cols;
        data.resize(rows, vector<double>(cols));
    }

    explicit Matrix(vector<vector<double>> myMatrix) {
        this->data = myMatrix;
        this->rows = myMatrix.size();
        this->columns = myMatrix[0].size();
    }

    // addition (+) overloading
    Matrix operator+(Matrix &other) {
        Matrix result(rows, columns);
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[0].size(); ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // subtraction (-) overloading
    Matrix operator-(Matrix &other) {
        Matrix result(data.size(), data[0].size());
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[0].size(); ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    // multiplication (*) overloading
    Matrix operator*(Matrix &other) {
        Matrix result(data.size(), other.data[0].size());
        for (int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < other.data[0].size(); ++j) {
                double temp = 0;
                for (int k = 0; k < other.data.size(); ++k) {
                    temp += data[i][k] * 1.0 * other.data[k][j];
                }
                result.data[i][j] = temp;
//                if (j + 1 == data.size())
//                    cout << fixed << setprecision(4) << temp;
//                else
//                    cout << fixed << setprecision(4) << temp << " ";
//            }
//            cout << "\n";
            }
        }
        return result;
    }

    // transpose function
    static Matrix transpose(vector<vector<double>> A, int a, int b) {
        vector<double> t(a, 0);
        vector<vector<double>> T(b, t);

        for (int i = 0; i < a; ++i) {
            vector<double> temp;
            for (int j = 0; j < b; ++j) {
                T[j][i] = A[i][j];
            }
        }

        return (Matrix) T;
    }

    void changeIJ(int i, int j, double newValue) {
        data[i][j] = newValue;
    }

    double getIJ(int i, int j) {
        return data[i][j];
    }

    int correctPivot(int i) {
        double mx = -10000000;
        int index = -1;
        for (int j = i; j < data.size(); ++j) {
            if (fabs(data[j][i]) > mx) {
                mx = abs(data[j][i]);
                index = j;
            }
        }
        return index;
    }

};

class ColumnVector : public Matrix {
public:
    explicit ColumnVector(int n) {
        this->rows = n;
        this->columns = 1;
        data.resize(n, vector<double>(1));
    }

    ColumnVector(int rows, int columns) : Matrix(rows, columns) {
        this->rows = rows;
        this->columns = 1;
        data.resize(rows, vector<double>(1));
    }

    void changeIJ(int i, double newValue) {
        data[i][0] = newValue;
    }

    void swapColumnVector(int i, int j) {
        swap(data[i][0], data[j][0]);
    }

    double getIJ(int i) {
        return data[i][0];
    }

};

// class for storing permutation matrix
class PermutationMatrix : public Matrix {
public:
    // constructor for permutation matrix
    explicit PermutationMatrix(int n) {
        this->rows = n;
        this->columns = n;
        data.resize(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    data[i][j] = 1;
                }
            }
        }
    }

    // function for swapping 2 rows
    void permute(Matrix &A, int firstElementToChange, int secondElementToChange) {
        swap(data[firstElementToChange], data[secondElementToChange]);
    }
};

// class for storing elimination matrix
class EliminationMatrix : public Matrix {
public:
    // constructor for elimination matrix
    explicit EliminationMatrix(int n) {
        this->rows = n;
        this->columns = n;
        data.resize(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    data[i][j] = 1;
                }
            }
        }
    }

    // function for changing elimination matrix
    void eliminate(Matrix &A, int firstElementToChange, int secondElementToChange) {
        double factor = (A.getIJ(firstElementToChange, secondElementToChange) * 1.0) /
                        A.getIJ(secondElementToChange, secondElementToChange);
        factor *= -1;
        data[firstElementToChange][secondElementToChange] = factor;
    }
};

class IdentityMatrix : public Matrix {
public:
    explicit IdentityMatrix(int n) {
        this->rows = n;
        this->columns = n;
        data.resize(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    data[i][j] = 1;
                }
            }
        }
    }

};

// class for Gauss Elimination
class GaussElimination {

public:

    Matrix gaussElimination(int n, Matrix A, Matrix B) {
        int steps = 1;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) {
                    int k = A.correctPivot(i);
                    if (k != i) {
                        PermutationMatrix permutationMatrix(n);
                        permutationMatrix.permute(A, i, k);
                        A = permutationMatrix * A;
                        B = permutationMatrix * B;
                    }
                }
                if (i > j) {
                    if (A.getIJ(i, j) != 0.00) {
                        EliminationMatrix eliminationMatrix(n);
                        eliminationMatrix.eliminate(A, i, j);
                        A = eliminationMatrix * A;
                        B = eliminationMatrix * B;

                    }
                }
            }
        }

        for (int j = n - 1; j >= 0; j--) {
            for (int i = n - 1; i >= 0; --i) {
                if (i < j) {
                    if (A.getIJ(i, j) == 0) continue;
                    EliminationMatrix eliminationMatrix(n);
                    eliminationMatrix.eliminate(A, i, j);
                    A = eliminationMatrix * A;
                    B = eliminationMatrix * B;

                }
            }
        }

        double norm = 1;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    norm = A.getIJ(i, j);
                    A.changeIJ(i, j, A.getIJ(i, j) / norm);
                } else {
                    A.changeIJ(i, j, A.getIJ(i, j) / norm);
                }
            }
            for (int j = 0; j < n; ++j) {
                B.changeIJ(i, j, B.getIJ(i, j) / norm);
            }
        }
        return B;
    }
};


int main() {
    int matrixSize, degreeOfEquation; // size of matrix
    vector<vector<double>> matrixA;
    cin >> matrixSize;
    vector<double>xValues(matrixSize), yValues(matrixSize);
//    double xValues[matrixSize];
//    double yValues[matrixSize];

    for (int i = 0; i < matrixSize; ++i) {
        double x, y;
        cin >> x >> y;
        xValues[i] = x;
        yValues[i] = y;
    }
    cin >> degreeOfEquation;
    ColumnVector columnVector(matrixSize);
    Matrix A(matrixSize, degreeOfEquation + 1);
    for (int i = 0; i < matrixSize; ++i) {
        matrixA.emplace_back(degreeOfEquation + 1, 0);
        for (int j = degreeOfEquation; j >= 0; --j) {
            matrixA[i][j] = pow(xValues[i], j);
            A.changeIJ(i, j, pow(xValues[i], j));
        }
        columnVector.changeIJ(i, yValues[i]);
    }

    cout << "A:" << endl;
    cout << A;
    cout << "A_T*A:" << endl;
    Matrix ATransposed = Matrix::transpose(matrixA, matrixSize, degreeOfEquation + 1);
    Matrix A_T_multipliedBy_A = ATransposed * A;
    cout << A_T_multipliedBy_A;
    GaussElimination gaussElimination;
    IdentityMatrix identityMatrix(degreeOfEquation + 1);
    Matrix A_T_multipliedBy_A_inverse = gaussElimination.gaussElimination(degreeOfEquation + 1, A_T_multipliedBy_A,
                                                                          identityMatrix);
    cout << "(A_T*A)^-1:" << endl;
    cout << A_T_multipliedBy_A_inverse;

    cout << "A_T*b:" << endl;
    Matrix A_T_mult_b = ATransposed * columnVector;
    cout << A_T_mult_b;
    cout << "x~:" << endl;
    Matrix result = A_T_multipliedBy_A_inverse * A_T_mult_b;
    cout << result;

    FILE *pipe = popen(GNUPLOT_NAME, "w");


    if (pipe != nullptr) {
        const double nPoints = 200;
        const double step = 0.1;

        int  xStart = 0, xEnd = 20, yStart = 0, yEnd = 20;

        fprintf(pipe, "set datafile separator '\t'\n");
        ::fprintf(pipe, "set yrange [%f:%f]\n", yStart - 7.0, yEnd +
                                                              7.0);
        ::fprintf(pipe, "set xrange [%f:%f]\n", xStart - 7.0, xEnd +
                                                              7.0);

        fprintf(pipe, "%s\n", "f(x) = 2.000 * x**0 + " "0.7500 * x**1 + " "0.1250 * x**2");
        ::fprintf(pipe, "%s\n",
                  "plot '-' using 1:2 title 'points' with points, "
                  "f(x) title 'least square approximation' with lines");

        for (int i = 0; i < matrixSize; ++i) {
            double x = xValues[i];
            double y = yValues[i];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "e\n");
        fflush(pipe);

        pclose(pipe);
    } else{
        cout << "Couldn't open pipe" << "\n";
    }
    return 0;
};




