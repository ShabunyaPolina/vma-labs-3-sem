// Шабуня Полина (1 гр)
// лаб1 метод Гаусса

#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

// ввод коэффициентов уравнений СЛАУ
void fill_system(std::vector<std::vector<double>> &a, std::vector<double> &b) {
    int n = 0;
    std::cout << "Количество уравнений: ";
    std::cin >> n;
    double tmp = 0;
    for(int i = 0; i < n; ++i) {
        std::vector<double> vec;
        for (int j = 0; j < n; ++j) {
            std::cout << "a[" << i << ',' << j << "] = ";
            std::cin >> tmp;
            vec.push_back(tmp);
        }
        a.push_back(vec);
        std::cout << "b[" << i << "] = ";
        std::cin >> tmp;
        b.push_back(tmp);
    }
}

// вывод расширенной матрицы СЛАУ
void show_system(std::vector<std::vector<double>> &a, std::vector<double> &b) {
    int n = (int)a.size();
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << a[i][j] << '\t';
        }
        std::cout << "|\t" << b[i] << '\n';
    }
}

// метод Гаусса без выбора ведущего элемента
void solve(std::vector<std::vector<double>> &a,
           std::vector<double> &b, std::vector<double> &x) {
    int n = (int)a.size();
    x.resize(n);
    // прямой ход
    double l = 0;
    for(int k = 0; k < n - 1; ++k) {  // шаги
        for(int i = k + 1; i < n; ++i) {
            l = a[i][k] / a[k][k];
            std::cout << "l = " << l << '\n';
            for(int j = k; j < n; ++j) {
                a[i][j] -= l * a[k][j];
            }
            b[i] -= l * b[k];
        }
    }
    // обратный ход
    double s, xi;
    for(int i = n - 1; i >= 0; --i) {
        s = 0;
        for (int j = i; j < n; ++j) {
            s += a[i][j] * x[j];
            for(double q : x) {
                std::cout << q << ' ';
            }
            std::cout <<"s = " << s << '\n';
        }
        xi = (b[i] - s) / a[i][i];
        std::cout << xi << '\n';
        x[i] = xi;
    }
}

// генерирует СЛАУ размерности 10
void generate_system(std::vector<std::vector<double>> &a, std::vector<double> &f, int n) {
    std::vector<double> x(n);
    a.assign(n, std::vector<double>(n));
    f.resize(n);
    srand ( time(nullptr) );
    const int MIN_RAND = -100;
    const int MAX_RAND = 100;
    // генерируем матрицу А и вектор Х
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = MIN_RAND + rand() % (MAX_RAND - MIN_RAND + 1);
        }
        x[i] = MIN_RAND + rand() % (MAX_RAND - MIN_RAND + 1);
    }
    // высчитываем столбец f
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            f[i] += a[i][j] * x[j];
        }
    }
}

int main() {
    std::vector<std::vector<double>> a;
    std::vector<double> b;
    std::vector<double> x;

//    fill_system(a, b);
//    show_system(a, b);
//    solve(a, b, x);
//    show_system(a, b);
//    for(double i : x) {
//        std::cout << i << ' ';
//    }

    generate_system(a, b, 10);
    show_system(a,b);
    return 0;
}
