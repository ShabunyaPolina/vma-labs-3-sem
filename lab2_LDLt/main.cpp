// Шабуня Полина (1 гр)
// лаб2 LDLt разложение

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

// ввод матрицы
void fill_matrix(std::vector<std::vector<float>> &a) {
    int n = 0;
    std::cout << "n = ";
    std::cin >> n;
    a.assign(n, std::vector<float>(n));
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << "a[" << i << ',' << j << "] = ";
            std::cin >> a[i][j];
        }
    }
}

// вывод мвтрицы
void show_matrix(std::vector<std::vector<float>> &a) {
    int n = (int)a.size();
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << a[i][j] << '\t';
        }
        std::cout << '\n';
    }
}

// вывод расширенной матрицы СЛАУ
void show_system(std::vector<std::vector<float>> &a, std::vector<float> &b) {
    int n = (int)a.size();
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << a[i][j] << '\t';
        }
        std::cout << "|\t" << b[i] << '\n';
    }
}

// LDLt разложение матрицы
void decompose(std::vector<std::vector<float>> &a) {
    int n = (int)a.size();
    std::vector<float> tmp(n*n);
    for(int k = 0; k < n - 1; ++k) {
        for(int i = k + 1; i < n; ++i) {
            tmp[i] = a[i][k];
            a[i][k] /= (float)a[k][k];
            for(int j = k + 1; j <= i; ++j) {
                a[i][j] -= a[i][k] * tmp[j];
            }
        }
    }
}

// ФУНКЦИИ ДЛЯ ПРОВЕРКИ РАЗЛОЖЕНИЯ

// выделяет матрицу L
void toL(std::vector<std::vector<float>> &a, std::vector<std::vector<float>> &l) {
    int n = (int)a.size();
    l.assign(n, std::vector<float>(n));
    for(int i = 0; i < n; ++i) {
        l[i][i] = 1;
    }
    for(int i = 1; i < n; ++i) {
        for(int j = 0; j < i; ++j) {
            l[i][j] = a[i][j];
        }
    }
}

// выделяет матрицу D
void toD(std::vector<std::vector<float>> &a, std::vector<std::vector<float>> &d) {
    int n = (int)a.size();
    d.assign(n, std::vector<float>(n));
    for(int i = 0; i < n; ++i) {
        d[i][i] = a[i][i];
    }
}

// транспонирует матрицу
void transpose(std::vector<std::vector<float>> &x) {
    int n = (int)x.size();
    for(int i = 0; i < n; ++i) {
        for(int j = i + 1; j < n; ++j) {
            std::swap(x[i][j], x[j][i]);
        }
    }
}

// перемножает матрицы
void multiply(std::vector<std::vector<float>> &a,
              std::vector<std::vector<float>> &b,
              std::vector<std::vector<float>> &res) {
    int n = (int)a.size();
    res.assign(n, std::vector<float>(n));
    for(int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res[k][i] += a[k][j] * b[j][i];
            }
        }
    }
}

// генерирует матрицу А размерности n*n
void generate_matrix(std::vector<std::vector<float>> &a, int n) {
    a.assign(n, std::vector<float>(n));
    srand ( time(nullptr) );
    const int MIN_RAND = -4;
    const int MAX_RAND = 0;
    // генерация недиагональных элементов матрицы А
    for(int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            a[i][j] = a[j][i] = MIN_RAND + rand() % (MAX_RAND - MIN_RAND + 1);
        }
    }
    // вычисление диагональных элементов матрицы А
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if(j == i) continue;
            a[i][i] -= a[i][j];
        }
    }
}

// вычисляет a[0][0] в зависимости от k
void set_a00(std::vector<std::vector<float>> &a, int k) {
    a[0][0] += (float)pow(10, -1 * k);
}

// вычисляет столбец b СЛАУ
void calculate_b(std::vector<std::vector<float>> &a,
                 std::vector<float> &b,
                 std::vector<float> &x, int m) {
    // генерация вектора X
    int n = (int)a.size();
    x.resize(n);
    b.resize(n);
    for(int i = 0; i < n; ++i) {
        x[i] = (float)(m + i);
    }
    // вычисление столбца b
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += a[i][j] * x[j];
        }
    }
}

// решает систему
void solve_system(std::vector<std::vector<float>> &a,
                  std::vector<float> &b,
                  std::vector<float> &x) {
    int n = (int)a.size();
    x.resize(n);
    std::vector<float> tmp(n);
    // решение системы L*tmp=b, где tmp = DLt
    for(int i = 0; i < n; ++i) {
        tmp[i] = b[i];
        for(int j = 0; j < i; ++j) {
            tmp[i] -= a[i][j] * tmp[j];
        }
    }
    // решение системы DLt*x=tmp
    for(int i = n - 1; i >= 0; --i) {
        x[i] = tmp[i] / a[i][i];
        for(int j = n - 1; j > i; --j) {
            x[i] -= a[j][i] * x[j];
        }
    }
}

// вычисляет относительную погрешность в кубической норме
float calculate_error(std::vector<float> &x, std::vector<float> &true_x) {
    float tmp = 0, max_delta = 0, max_true_x = 0;
    for(int i = 0; i < x.size(); ++i) {
        tmp = std::fabs(x[i] - true_x[i]);
        if(tmp > max_delta)
            max_delta = tmp;
        if(std::fabs(x[i]) > max_true_x)
            max_true_x = std::fabs(x[i]);
    }
    return max_delta / max_true_x;
}

int main() {
//    std::vector<std::vector<float>> a;
//    std::vector<std::vector<float>> a1;
//    std::vector<float> b;
//    std::vector<float> b1;
//    std::vector<float> x;
//    std::vector<float> x1;
//    std::vector<float> x2;
//
//    generate_matrix(a, 12);
//    a1 = a;
//
//    set_a00(a,0);
//    set_a00(a1,1);
//    calculate_b(a, b, x, 10);
//    calculate_b(a1, b1, x, 10);
//    std::cout << "СЛАУ:\n";
//    std::cout << "\nk = 0\n";
//    show_system(a, b);
//    std::cout << "\nk = 1\n";
//    show_system(a1, b1);
//    std::cout << "\nПреобразованная матрица А (нижний треугольник - L "
//                 "(за исключением главной диагонали), диагональ - D):\n";
//
//    std::cout << "\nk = 0\n";
//    decompose(a);
//    show_matrix(a);
//    std::cout << "\nk = 1\n";
//    decompose(a1);
//    show_matrix(a1);
//
//    solve_system(a,b,x1);
//    solve_system(a1,b1,x2);
//
//    std::cout << "\nВектор решений систем:\n";
//    for(auto &i : x) {
//        std::cout << i << ' ';
//    }
//    std::cout << '\n';
//
//    std::cout << "\nВекторы решений системы, полученные решением системы на основе LDLt разложения:";
//    std::cout << "\nk = 0\n";
//    for(auto &i : x1) {
//        std::cout << i << ' ';
//    }
//    std::cout << "\nk = 1\n";
//    for(auto &i : x2) {
//        std::cout << i << ' ';
//    }
//    std::cout << '\n';
//
//    float er1 = calculate_error(x1, x);
//    float er2 = calculate_error(x2, x);
//    std::cout << "\nПогрешность метода:\n";
//    std::cout << "k = 0:  " << er1;
//    std::cout << "\nk = 1:  " << er2;

//    std::vector<std::vector<float>> a{{1,-1,1},{-1,0.001, 0.001},{1,0.001,0.001}};
    std::vector<std::vector<float>> a{{5,7,9},{5,3,7},{9,5,5}};
    std::vector<std::vector<float>> l;
    std::vector<std::vector<float>> d;
    show_matrix(a);
    decompose(a);
    toL(a, l);
    toD(a, d);
    std::cout << std::endl;
    show_matrix(l);
    //show_matrix(d);


    return 0;
}
