#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

const int MAX = 4;
const double EPS = 1e-9;

// 행렬 출력
void printMatrix(double matrix[MAX][MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(10) << fixed << setprecision(2)
            << matrix[i][j];
        }
        cout << endl;
    }
}

// 다음 단계로 넘어가기
void waitForEnter()
{
    cout << "\nEnter를 누르면 다음 단계로 넘어갑니다...";
    cin.ignore(10000, '\n');
    cin.get();
}

// 행렬식 계산
double determinant(double matrix[MAX][MAX], int n)
{
    if (n == 1)
        return matrix[0][0];
    if (n == 2)
    {
        return matrix[0][0] * matrix[1][1]
             - matrix[0][1] * matrix[1][0];
    }
    double det = 0;
    for (int col = 0; col < n; col++)
    {
        double minor[MAX][MAX];
        int r = 0;
        for (int i = 1; i < n; i++)
        {
            int c = 0;
            for (int j = 0; j < n; j++)
            {
                if (j == col)
                    continue;
                minor[r][c] = matrix[i][j];
                c++;
            }
            r++;
        }
        double sign = (col % 2 == 0) ? 1 : -1;
        det += sign * matrix[0][col]
             * determinant(minor, n - 1);
    }
    return det;
}

// 행렬식을 이용한 역행렬
bool inverseByDeterminant(
    double matrix[MAX][MAX],
    double inverse[MAX][MAX],
    int n)
{
    double det = determinant(matrix, n);
    if (fabs(det) < EPS)
        return false;
    if (n == 1)
    {
        inverse[0][0] = 1.0 / matrix[0][0];
        return true;
    }
    double cofactor[MAX][MAX];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double minor[MAX][MAX];
            int r = 0;
            for (int row = 0; row < n; row++)
            {
                if (row == i)
                    continue;
                int c = 0;
                for (int col = 0; col < n; col++)
                {
                    if (col == j)
                        continue;
                    minor[r][c] = matrix[row][col];
                    c++;
                }
                r++;
            }
            double sign = ((i + j) % 2 == 0) ? 1 : -1;
            cofactor[i][j]
                = sign * determinant(minor, n - 1);
        }
    }

    // 수반행렬을 이용해서 역행렬 계산
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            inverse[i][j] = cofactor[j][i] / det;
        }
    }

    return true;
}

// Gauss-Jordan으로 역행렬 계산
bool inverseByGaussJordan(
    double matrix[MAX][MAX],
    double inverse[MAX][MAX],
    int n)
{
    double augmented[MAX][MAX * 2];

    // [A | I] 만들기
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            augmented[i][j] = matrix[i][j];
            augmented[i][j + n] = (i == j) ? 1 : 0;
        }
    }

    for (int col = 0; col < n; col++)
    {
        int pivotRow = col;

        for (int row = col + 1; row < n; row++)
        {
            if (fabs(augmented[row][col])
                > fabs(augmented[pivotRow][col]))
            {
                pivotRow = row;
            }
        }

        if (fabs(augmented[pivotRow][col]) < EPS)
            return false;

        // 필요한 경우 행을 바꿈
        if (pivotRow != col)
        {
            for (int j = 0; j < 2 * n; j++)
            {
                swap(augmented[col][j],
                    augmented[pivotRow][j]);
            }
        }

        // 피벗을 1로 만들기
        double pivot = augmented[col][col];

        for (int j = 0; j < 2 * n; j++)
        {
            augmented[col][j] /= pivot;
        }

        // 다른 행의 해당 열을 0으로 만들기
        for (int row = 0; row < n; row++)
        {
            if (row == col)
                continue;

            double factor = augmented[row][col];

            for (int j = 0; j < 2 * n; j++)
            {
                augmented[row][j]
                    -= factor * augmented[col][j];
            }
        }
    }

    // 오른쪽 부분이 역행렬
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            inverse[i][j]
                = augmented[i][j + n];
        }
    }

    return true;
}

// 두 행렬 비교
bool compareMatrices(
    double A[MAX][MAX],
    double B[MAX][MAX],
    int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fabs(A[i][j] - B[i][j]) > EPS)
                return false;
        }
    }

    return true;
}

// 역행렬 확인
bool verifyInverse(
    double matrix[MAX][MAX],
    double inverse[MAX][MAX],
    int n)
{
    double result[MAX][MAX] = {0};

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j]
                    += matrix[i][k] * inverse[k][j];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double expected = (i == j) ? 1.0 : 0.0;

            if (fabs(result[i][j] - expected) > EPS)
                return false;
        }
    }

    return true;
}


int main()
{
    int n;

    // 1단계
    cout << "[1단계] 행렬 크기 입력" << endl;
    cout << "n = ";
    cin >> n;

    if (n < 1 || n > MAX)
    {
        cout << "잘못된 크기입니다." << endl;
        return 0;
    }

    waitForEnter();

    // 2단계
    double matrix[MAX][MAX];

    cout << "\n[2단계] 행렬 입력" << endl;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> matrix[i][j];
        }
    }

    cout << "\n입력한 행렬" << endl;
    printMatrix(matrix, n);

    waitForEnter();

    // 3단계
    cout << "\n[3단계] 행렬식 계산" << endl;

    double det = determinant(matrix, n);

    cout << "Det(A) = "
        << fixed << setprecision(2)
        << det << endl;

    waitForEnter();

    // 4단계
    cout << "\n[4단계] 행렬식을 이용한 역행렬" << endl;

    double inverseDet[MAX][MAX];
    bool detResult =
        inverseByDeterminant(matrix, inverseDet, n);

    if (detResult)
    {
        printMatrix(inverseDet, n);
    }
    else
    {
        cout << "행렬식이 0이므로 역행렬이 없습니다."
            << endl;
    }

    waitForEnter();

    // 5단계
    cout << "\n[5단계] Gauss-Jordan" << endl;

    double inverseGauss[MAX][MAX];

    bool gaussResult =
        inverseByGaussJordan(matrix, inverseGauss, n);

    if (gaussResult)
    {
        printMatrix(inverseGauss, n);
    }
    else
    {
        cout << "Gauss-Jordan으로 역행렬을 구할 수 없습니다."
            << endl;
    }

    waitForEnter();

    // 6단계
    cout << "\n[6단계] 두 결과 비교" << endl;

    if (detResult && gaussResult)
    {
        if (compareMatrices(
                inverseDet, inverseGauss, n))
        {
            cout << "두 결과가 같습니다." << endl;
        }
        else
        {
            cout << "두 결과가 다릅니다." << endl;
        }
    }
    else
    {
        cout << "역행렬이 없어서 비교할 수 없습니다."
            << endl;
    }

    waitForEnter();

    // 7단계
    cout << "\n[7단계] 역행렬 검증" << endl;

    if (detResult && gaussResult)
    {
        if (verifyInverse(
                matrix, inverseDet, n))
        {
            cout << "A × A^-1 = I 확인 완료."
                << endl;
        }
        else
        {
            cout << "검증에 실패했습니다."
                << endl;
        }
    }
    else
    {
        cout << "역행렬이 없어서 검증할 수 없습니다."
            << endl;
    }

    cout << "\n프로그램 종료" << endl;

    return 0;
}