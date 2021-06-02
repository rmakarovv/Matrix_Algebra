#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

/**
 * The class for representing matrices
 * Includes:
 * Overloaded operators '+', '-', '*', transposing
 * Overloaded input and output.
 */
class Matrix {
public:
    int n; // number of rows
    int m; // number of columns
    double M[100][100]; // A way to store a matrix

    Matrix(int nToSet, int mToSet) {
        n = nToSet;
        m = mToSet;
    }

    double determinant() {
        double result = 1;
        int step = 1;
        string previousStep = "None";
        for (int i = 0; i < n; i++) {
            int pivotRow = i; double pivot = M[i][i];
            for (int j = i + 1; j < n; ++j) {
                if (abs(M[j][i]) > abs(pivot)) {
                    pivot = M[j][i];
                    pivotRow = j;
                }
            }
            if (pivot == 0.0) return 0;
            else result *= pivot;
            if (pivotRow != i) {
                result *= -1;
                previousStep = "permutation";
                swap(M[i], M[pivotRow]);
                pivotRow = i;
            }

            for (int row = pivotRow + 1; row < n; row++) {
                if (M[row][pivotRow] != 0){
                    previousStep = "elimination";
                    step++;
                    double diff = M[row][pivotRow] / pivot;
                    for (int del = 0; del < n; del++) {
                        M[row][del] -= diff * M[pivotRow][del];
                    }
                }
            }
        }
        return result;
    }
    /*
     * Overloading '>>' operator.
     */
    friend istream & operator>>(istream&  is, Matrix& Matrix) {
        for (int i = 0; i < Matrix.n; i++) {
            for (int j = 0; j < Matrix.m; j++)
                is >> Matrix.M[i][j];
        }
        return is;
    }

    /*
     * Overloading '<<' operator.
     */
    friend ostream& operator<< (std::ostream& out, const Matrix& A) {
        for (int i = 0; i < A.n; i++) {
            for (int j = 0; j < A.m; j++) {
                out << A.M[i][j];
                if (j != A.m - 1)
                    out << " ";
            }
            out << "\n";
        }
        return out;
    }

    /*
     * Overloading '=' operator.
     */
    Matrix& operator=(const Matrix& Matrix) {
        n = Matrix.n;
        m = Matrix.m;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                M[i][j] = Matrix.M[i][j];
        }
        return *this;
    }

    /*
     * Overloading '+' operator.
     */
    Matrix operator+(const Matrix& B) {
        if (m != B.m || n != B.n){
            cout << "Error: the dimensional problem occurred\n";
        } else {
            Matrix D(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++)
                    D.M[i][j] = M[i][j] + B.M[i][j];
            }
            return D;
        }
        return Matrix(0, 0);
    }

    /*
     * Overloading '-' operator.
     */
    Matrix operator-(const Matrix& A) {
        if (m != A.m || n != A.n){
            cout << "Error: the dimensional problem occurred\n";
        } else {
            Matrix E(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++)
                    E.M[i][j] = M[i][j] - A.M[i][j];
            }
            return E;
        }
        return Matrix(0, 0);
    }

    /*
     * Overloading '*' operator.
     */
    Matrix operator*(const Matrix& A) {
        if (this->m != A.n) {
            cout << "Error: the dimensional problem occurred\n";
        } else {
            Matrix F(this->n, A.m);
            for (int i = 0; i < this->n; i++) {
                for (int j = 0; j < A.m; j++) {
                    F.M[i][j] = 0;
                    for (int k = 0; k < this->m; k++)
                        F.M[i][j] += M[i][k] * A.M[k][j];
                }
            }
            return F;
        }
        return Matrix(0, 0);
    }

    /*
     * A function for getting a transpose matrix.
     */
    Matrix Transpose() {
        Matrix G(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                G.M[i][j] = M[j][i];
            }
        }
        return G;
    }

    /*
     * Making an identity matrix from an initial.
     */
    void makeMatrixIdentity(){
        for (int i = 0; i < this->n; i++) {
            this->M[i][i] = 1;
        }
    }

    /*
     * Making an elimination matrix from an initial.
     */
    void makeMatrixElimination(int elRow, int elCol, const Matrix &toEliminate) {
        makeMatrixIdentity();
        elRow--; elCol--;
        this->M[elRow][elCol] =  - toEliminate.M[elRow][elCol] / toEliminate.M[elCol][elCol];
    }

    /*
     * Making an permutation matrix from an initial.
     */
    void makeMatrixPermutation(int swap1, int swap2, const Matrix &toSwap) {
        swap1--; swap2--;
        makeMatrixIdentity();
        swap(M[swap2][swap1], M[swap1][swap1]);
        swap(M[swap2][swap2], M[swap1][swap2]);
    }
};

/**
 * The class for storing square matrices.
 */
class SquareMatrix : public Matrix {
    // A class to storing square Matrices
public:
    int n;

    // Constructor with a call to Matrix constructor
    SquareMatrix(int rows) : Matrix(rows, rows) {
        n = rows;
    }
};

/**
 * The class for storing identity matrices.
 */
class IdentityMatrix : public SquareMatrix {
    // A class to store identityMatrix
public:
    int n;
    int M[100][100] = {};
    IdentityMatrix(int number) : SquareMatrix(number) {
        n = number;
        makeMatrixIdentity();
    }
};

/**
 * The class for storing elimination matrices.
 */
class EliminationMatrix : public Matrix {
    // A class to store eliminationMatrix
public:
    int n;
    int eliminateRow;
    int eliminateColumn;

    EliminationMatrix(int rows, int elRow, int elCol, Matrix toEliminate) : Matrix(rows, rows) {
        n = rows;
        eliminateRow = elRow;
        eliminateColumn = elCol;
        makeMatrixElimination(elRow, elCol, toEliminate);
    }
};

/**
 * The class for storing permutation matrices.
 */
class PermutationMatrix : public Matrix {
    // A class to store permutationMatrix
public:
    int n;

    PermutationMatrix(int rows, int swap1, int swap2, Matrix m) : Matrix(rows, rows) {
        n = rows;
        makeMatrixPermutation(swap1, swap2, m);
    }
};

/**
 * The class for representing column vectors.
 */
class ColumnVector {
public:
    double cv[100];
    int n;

    ColumnVector(int N) {
        n = N;
    }

    /*
     * Overloading '>>' operator.
     */
    friend istream & operator>>(istream&  is, ColumnVector &columnVector) {
        for (int i = 0; i < columnVector.n; ++i) {
            is >> columnVector.cv[i];
        }
        return is;
    }

    /*
     * Overloading '<<' operator.
     */
    friend ostream& operator<< (std::ostream& out, const ColumnVector &columnVector) {
        for (int i = 0; i < columnVector.n; ++i) {
            out << columnVector.cv[i] << "\n";
        }
        return out;
    }

};


/**
 * The function that for a given matrix returns inverse.
 */
Matrix inverseMatrix(const Matrix &A) {
    Matrix augmented(A.n, A.m * 2);
    for (int i = 0; i < A.n; ++i) {
        for (int j = 0; j < A.m; ++j) {
            augmented.M[i][j] = A.M[i][j];
        }
    }
    for (int i = 0; i < A.n; ++i) {
        for (int j = A.m; j < A.m * 2; ++j) {
            if (i == j - A.m)
                augmented.M[i][j] = 1;
            else augmented.M[i][j] = 0;
        }
    }
    int step = 0;
    step++;
    for (int i = 0; i < augmented.n; i++) {
        int pivotRow = i; double pivot = augmented.M[i][i];
        for (int j = i + 1; j < augmented.n; ++j) {
            if (abs(augmented.M[j][i]) > abs(pivot)) {
                pivot = augmented.M[j][i];
                pivotRow = j;
            }
        }

        if (pivotRow != i) {
            step++;
            swap(augmented.M[i], augmented.M[pivotRow]);
            pivotRow = i;
        }

        for (int row = pivotRow + 1; row < augmented.n; row++) {
            if (augmented.M[row][pivotRow] != 0) {
                step++;
                double diff = augmented.M[row][pivotRow] / pivot;
                for (int del = 0; del < augmented.m; del++){
                    augmented.M[row][del] -= diff * augmented.M[pivotRow][del];
                }
            }
        }
    }
    for (int i = augmented.n - 1; i >= 0; i--) {
        double pivot = augmented.M[i][i];
        for (int j = i - 1; j >= 0; j--) {
            if (augmented.M[j][i] != 0) {
                step++;
                double difference = augmented.M[j][i] / pivot;
                for (int k = 0; k < augmented.m; k++){
                    augmented.M[j][k] -= difference * augmented.M[i][k];
                }
            }
        }
    }

    for (int row = 0; row < augmented.n; ++row) {
        double coef = 1.0;
        if (augmented.M[row][row] != 0)
            coef = 1.0 / augmented.M[row][row];
        for (int col = 0; col < augmented.m; ++col)
            augmented.M[row][col] *= coef;
    }


    Matrix inverse(augmented.n, augmented.n);
    for (int row = 0; row < inverse.n; ++row) {
        for (int col = 0; col < inverse.m; ++col) {
            inverse.M[row][col] = augmented.M[row][col + augmented.n];
        }
    }
    return inverse;
}


/**
 * The process of solving equations Ax = b, with return value x.
 */
vector<double> solve(const Matrix &A, ColumnVector &columnVector) {
    // Solving equation Ax = b
    Matrix augmented(A.n, A.m + 1);
    for (int i = 0; i < A.n; ++i) {
        for (int j = 0; j < A.m; ++j) {
            augmented.M[i][j] = A.M[i][j];
        }
    }
    for (int i = 0; i < A.n; ++i) {
        augmented.M[i][A.m] = columnVector.cv[i];
    }
    int step = 0;
    step++;
    for (int i = 0; i < augmented.n; i++) {
        int pivotRow = i; double pivot = augmented.M[i][i];
        for (int j = i + 1; j < augmented.n; ++j) {
            if (abs(augmented.M[j][i]) > abs(pivot)) {
                pivot = augmented.M[j][i];
                pivotRow = j;
            }
        }

        if (pivotRow != i) {
            step++;
            swap(augmented.M[i], augmented.M[pivotRow]);
            pivotRow = i;
        }

        for (int row = pivotRow + 1; row < augmented.n; row++) {
            if (abs(augmented.M[row][pivotRow]) > 0.001) {
                step++;
                double diff = augmented.M[row][pivotRow] / pivot;
                for (int del = 0; del < augmented.m; del++) {
                    augmented.M[row][del] -= diff * augmented.M[pivotRow][del];
                }
            }
        }
    }
    for (int i = augmented.n - 1; i >= 0; i--) {
        double pivot = augmented.M[i][i];
        for (int j = i - 1; j >= 0; j--) {
            if (abs(augmented.M[j][i]) > 0.001) {
                step++;
                double difference = augmented.M[j][i] / pivot;
                for (int k = 0; k < augmented.m; k++) {
                    augmented.M[j][k] -= difference * augmented.M[i][k];
                }
            }
        }
    }

    // Diagonal normalization
    for (int row = 0; row < augmented.n; ++row) {
        double coef = 1;
        if (augmented.M[row][row] != 0)
            coef = 1 / augmented.M[row][row];
        for (int col = 0; col < augmented.m; ++col)
            augmented.M[row][col] *= coef;
    }


    // Constructing a solution.
    vector<double> solutionX;
    for (int i = 0; i < A.n; ++i) {
        solutionX.push_back(augmented.M[i][A.m]);
    }
    return solutionX;
}


/**
 * The implementation for a Jacobi method.
 */
void Jacobi() {
    int sizeAn;
    cin >> sizeAn;
    SquareMatrix A(sizeAn);
    cin >> A;

    int sizeB;
    cin >> sizeB;
    Matrix B(sizeB, 1);
    cin >> B;

    double e;
    cin >> e;

    bool isApplicable = true;

    for (int i = 0; i < sizeAn; ++i) {
        double aii = abs(A.M[i][i]);
        double sum = 0;
        for (int j = 0; j < sizeAn; ++j) {
            if (i != j)
                sum += abs(A.M[i][j]);
        }
        if (sum >= aii)
            isApplicable = false;
    }

    if (isApplicable) {
        SquareMatrix D(sizeAn);
        for (int i = 0; i < sizeAn; ++i) {
            D.M[i][i] = A.M[i][i];
        }

        IdentityMatrix I(sizeAn);

        Matrix alpha = I - inverseMatrix(D) * A;

        Matrix beta = inverseMatrix(D) * B;

        cout << "alpha:\n";
        cout << fixed << setprecision(4) << alpha;

        cout << "beta:\n";
        cout << fixed << setprecision(4) << beta;

        int step = 0;

        Matrix xprev(alpha.n, 1);
        xprev = beta;

        cout << "x(" << step << "):\n";
        step++;
        cout << fixed << setprecision(4) << xprev;

        Matrix xcur(alpha.n, 1);
        xcur = beta;
        double curE = 1;

        while (curE > e) {
            xprev = xcur;
            curE = 0;
            xcur = alpha * xprev + beta;
            for (int i = 0; i < xcur.n; ++i) {
                curE += pow((xcur.M[i][0] - xprev.M[i][0]), 2);
            }
            curE = sqrt(curE);

            cout << "e: ";
            cout << fixed << setprecision(4) << curE << "\n";

            cout << "x(" << step << "):\n";
            step++;
            cout << fixed << setprecision(4) << xcur;
        }
    } else {
        cout << "The method is not applicable!\n";
    }
}


/**
 * The implementation for a Seidel's method.
 */
void Seidel() {
    int sizeAn;
    cin >> sizeAn;
    SquareMatrix A(sizeAn);
    cin >> A;

    int sizeb;
    cin >> sizeb;
    Matrix b(sizeb, 1);
    cin >> b;

    double e;
    cin >> e;

    bool isApplicable = true;

    for (int i = 0; i < sizeAn; ++i) {
        double aii = abs(A.M[i][i]);
        double sum = 0;
        for (int j = 0; j < sizeAn; ++j) {
            if (i != j)
                sum += abs(A.M[i][j]);
        }
        if (sum >= aii)
            isApplicable = false;
    }

    if (isApplicable) {
        SquareMatrix D(sizeAn);
        for (int i = 0; i < sizeAn; ++i) {
            D.M[i][i] = A.M[i][i];
        }

        IdentityMatrix I(sizeAn);

        Matrix alpha = I - inverseMatrix(D) * A;

        Matrix beta = inverseMatrix(D) * b;

        cout << "beta:\n";
        cout << fixed << setprecision(4) << beta;

        cout << "alpha:\n";
        cout << fixed << setprecision(4) << alpha;

        int step = 0;

        SquareMatrix B(sizeAn);
        for (int i = 0; i < sizeAn; ++i) {
            for (int j = 0; j < sizeAn; ++j) {
                if (i > j)
                    B.M[i][j] = alpha.M[i][j];
                else
                    B.M[i][j] = 0;
            }
        }
        cout << "B:\n";
        cout << fixed << setprecision(4) << B;

        SquareMatrix C(sizeAn);
        for (int i = 0; i < sizeAn; ++i) {
            for (int j = 0; j < sizeAn; ++j) {
                if (i < j)
                    C.M[i][j] = alpha.M[i][j];
                else
                    C.M[i][j] = 0;
            }
        }
        cout << "C:\n";
        cout << fixed << setprecision(4) << C;

        Matrix IB = I - B;
        cout << "I-B:\n";
        cout << fixed << setprecision(4) << IB;

        cout << "(I-B)_-1:\n";
        Matrix IB_1 = inverseMatrix(IB);
        cout << fixed << setprecision(4) << IB_1;

        Matrix xprev(alpha.n, 1);
        xprev = beta;

        cout << "x(" << step << "):\n";
        step++;
        cout << fixed << setprecision(4) << xprev;

        Matrix xcur(alpha.n, 1);
        xcur = beta;
        double curE = 1;

        while (curE > e) {
            xprev = xcur;
            curE = 0;
            xcur = IB_1 * C * xprev + IB_1 * beta;
            for (int i = 0; i < xcur.n; ++i) {
                curE += pow((xcur.M[i][0] - xprev.M[i][0]), 2);
            }
            curE = sqrt(curE);

            cout << "e: ";
            cout << fixed << setprecision(4) << curE << "\n";

            cout << "x(" << step << "):\n";
            step++;
            cout << fixed << setprecision(4) << xcur;
        }
    } else {
        cout << "The method is not applicable!\n";
    }
}


/**
 * The implementation of least square approximation
 */
void LeastSquareApproximation() {
    int m;
    cin >> m;
    vector<double> input;
    vector<double> b;
    for (int i = 0; i < m; ++i) {
        double in;
        cin >> in;
        input.push_back(in);
        cin >> in;
        b.push_back(in);
    }
    int n; cin >> n;
    int sizeInRows = m, sizeInColumns;
    sizeInColumns = 2 + (n - 1);
    Matrix A(sizeInRows, sizeInColumns);

    for (int i = 0; i < sizeInColumns; ++i) {
        for (int j = 0; j < sizeInRows; ++j)
            i == 0 ? A.M[j][i] = 1 : A.M[j][i] = pow(input[j], i);
    }

    Matrix B(m, 1);
    for (int i = 0; i < m; ++i) {
        B.M[i][0] = b[i];
    }

    Matrix A_T = A.Transpose();
    Matrix A_TA = A_T * A;

    Matrix A_TA1 = inverseMatrix(A_TA);

    Matrix A_Tb = A_T * B;

    Matrix A_TA1_b = A_TA1 * A_Tb;
    cout << "x~:\n";
    cout << fixed << setprecision(4) << A_TA1_b;
}
