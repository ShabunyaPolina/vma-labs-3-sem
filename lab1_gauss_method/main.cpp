// Шабуня Полина (1 гр)
// лаб1 метод Гаусса

#include <iostream>
#include <vector>
#include <algorithm>

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

// вывод расширенной матрицы А СЛАУ
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
    // прямой ход
    double l = 0;
    for(int k = 0; k < n - 1; ++k) {  // шаги
        for(int i = k + 1; i < n; ++i) {
            l = a[i][k] / a[k][k];
            for(int j = k; j < n; ++j) {
                a[i][j] -= l * a[k][j];
            }
            b[i] -= l * b[k];
        }
    }
    // обратный ход
    double s, xi;
    x.push_back(0);
    int v = 0;
    for(int i = n - 1; i >= 0; --i) {
        s = 0;
        for (int j = i; j < n; ++j) {
            s += a[i][j] * x[v];
            std::cout << s << '\n';
            ++v;
        }
        xi = (b[i] - s) / a[i][i];
        std::cout << xi << '\n';
        x.push_back(xi);
    }
    std::reverse(x.begin(), x.end());
    x.pop_back();
}

int main() {
    std::vector<std::vector<double>> a;
    std::vector<double> b;
    std::vector<double> x;
    fill_system(a, b);
    show_system(a, b);
    solve(a, b, x);
    show_system(a, b);
    for(double i : x) {
        std::cout << i << ' ';
    }
    return 0;
}
