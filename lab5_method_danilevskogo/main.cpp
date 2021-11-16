// Шабуня Полина (1 гр)
// лаб5 - метод Данилевского

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

// генерирует матрицу размерности n
std::vector<std::vector<float>> generate_matrix(int n, int left_border, int right_border) {
    std::vector<std::vector<float>> matrix(n, std::vector<float>(n));
    srand(time(nullptr));
    for (auto &row: matrix) {
        for (auto &element: row) {
            element = left_border + rand() % (2 * right_border + 1);
        }
    }
    return matrix;
}

// вывод матрицы
void show_matrix(std::vector<std::vector<float>> &matrix) {
    for (auto &row: matrix) {
        for (auto &element: row) {
            if (std::fabs(element) < 1e-8)
                std::cout << 0 << '\t';
            else
                std::cout << element << '\t';
        }
        std::cout << '\n';
    }
}

// проверяет отличие от нуля ведущего элемента (a[k][k-1] != 0)
bool is_valid_leading_element(std::vector<std::vector<float>> &matrix, int step) {
    int k = (int) matrix.size() - step;
    bool is_zero = std::fabs(matrix[k][k - 1]) < 1e-8;
    if (is_zero)
        std::cout << "Ведущий элемент равен 0.\n";
    return !is_zero;
}

// возвращает сгенерированную матрицу после проверки
std::vector<std::vector<float>> get_initial_matrix(int n, int left_border, int right_border) {
    std::vector<std::vector<float>> matrix;
    do {
        matrix = generate_matrix(n, left_border, right_border);
    } while (!is_valid_leading_element(matrix, 1));
    return matrix;
}

// возвращает трансформирующую матрицу M
std::vector<std::vector<float>> get_transform_matrix(std::vector<std::vector<float>> &A, int step) {
    int n = (int) A.size();
    std::vector<std::vector<float>> M(n, std::vector<float>(n));
    for (int i = 0; i < n; ++i) {
        if (i == n - 1 - step) {
            for (int j = 0; j < n; ++j) {
                if (j == n - 1 - step)
                    M[i][j] = 1 / A[n - step][n - 1 - step];
                else
                    M[i][j] = -A[n - step][j] / A[n - step][n - 1 - step];
            }
        } else
            M[i][i] = 1;
    }
    return M;
}

// возвращает обратную трансформирующую матрицу M^-1
std::vector<std::vector<float>> get_inverse_transform_matrix(std::vector<std::vector<float>> &A, int step) {
    int n = (int) A.size();
    std::vector<std::vector<float>> iM(n, std::vector<float>(n));
    for (int i = 0; i < n; ++i) {
        if (i == n - 1 - step) {
            for (int j = 0; j < n; ++j) {
                iM[i][j] = A[i + 1][j];
            }
        } else
            iM[i][i] = 1;
    }
    return iM;
}

// перемножает матрицы
std::vector<std::vector<float>> multiply(std::vector<std::vector<float>> &a,
                                         std::vector<std::vector<float>> &b) {
    int n = (int) a.size();
    std::vector<std::vector<float>> res(n, std::vector<float>(n));
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res[k][i] += a[k][j] * b[j][i];
            }
        }
    }
    return res;
}

// возвращает трансформированную матрицу после очередного шага
std::vector<std::vector<float>> get_transformed_matrix(std::vector<std::vector<float>> &A, int step) {
    int n = (int) A.size();
    std::vector<std::vector<float>> new_A(n, std::vector<float>(n));
    std::vector<std::vector<float>> M(n, std::vector<float>(n));
    new_A = get_inverse_transform_matrix(A, step);
    new_A = multiply(new_A, A);  // iM * A
    M = get_transform_matrix(A, step);
    new_A = multiply(new_A, M);  // iM * A * M
    return new_A;
}

// находит каноническую форму Фробениуса матрицы
bool get_Frobenius_form(std::vector<std::vector<float>> &A,
                        std::vector<std::vector<float>> &F) {
    int n = (int) A.size();
    std::vector<std::vector<float>> new_Ak = A;
    for (int i = 1; i < n; ++i) {
        new_Ak = get_transformed_matrix(new_Ak, i);
        if (!is_valid_leading_element(new_Ak, i))
            return false;
    }
    F = new_Ak;
    return true;
}

// возвращает коэффициент p1 матрицы Фробениуса
float get_p1(std::vector<std::vector<float>> &F) {
    return F[0][0];
}

// возвращает след матрицы
float get_trace(std::vector<std::vector<float>> &A) {
    int n = (int) A.size();
    float trace = 0;
    for (int i = 0; i < n; ++i) {
        trace += A[i][i];
    }
    return trace;
}


int main() {
    const int N = 4;  // порядок матрицы
    const int RANGE_BORDER = 50;  // для заполнения матрицы случайными числами

    bool is_finished = true;
    std::vector<std::vector<float>> A(N, std::vector<float>(N));
    std::vector<std::vector<float>> F(N, std::vector<float>(N));

    do {
        std::cout << "Матрица А:\n";
        A = get_initial_matrix(N, -RANGE_BORDER, RANGE_BORDER);
        show_matrix(A);

        std::cout << "\nКаноническая форма Фробениуса матрицы А:\n";

        if (!get_Frobenius_form(A, F)) {
            std::cout << "Нерегулярный случай.\n";
            is_finished = false;
        }
    } while (!is_finished);

    show_matrix(F);

    std::cout << "\nШаг 1:";
    std::vector<std::vector<float>> M3 = get_transform_matrix(A, 1);
    std::cout << "\nM3:\n";
    show_matrix(M3);
    std::vector<std::vector<float>> iM3 = get_inverse_transform_matrix(A, 1);
    std::cout << "\nM3^-1:\n";
    show_matrix(iM3);
    std::vector<std::vector<float>> A1 = get_transformed_matrix(A, 1);
    std::cout << "\nA1:\n";
    show_matrix(A1);

    std::cout << "\nШаг 2:";
    std::vector<std::vector<float>> M2 = get_transform_matrix(A1, 2);
    std::cout << "\nM3:\n";
    show_matrix(M2);
    std::vector<std::vector<float>> iM2 = get_inverse_transform_matrix(A1, 2);
    std::cout << "\nM3^-1:\n";
    show_matrix(iM2);
    std::vector<std::vector<float>> A2 = get_transformed_matrix(A1, 2);
    std::cout << "\nA2:\n";
    show_matrix(A2);

    std::cout << "\nШаг 3:";
    std::vector<std::vector<float>> M1 = get_transform_matrix(A2, 3);
    std::cout << "\nM3:\n";
    show_matrix(M1);
    std::vector<std::vector<float>> iM1 = get_inverse_transform_matrix(A2, 3);
    std::cout << "\nM3^-1:\n";
    show_matrix(iM1);
    std::vector<std::vector<float>> A3 = get_transformed_matrix(A2, 3);
    std::cout << "\nA3:\n";
    show_matrix(A3);

    std::cout << "\nКоэффициент p1: " << get_p1(F) << '\n' <<
              "След матрицы А (trsceА): " << get_trace(A) << '\n';

    return 0;
}
