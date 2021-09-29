// Шабуня Полина (1 гр)
// лаб1 метод Гаусса

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>

// ввод расширенной матрицы СЛАУ
void fill_system(std::vector<std::vector<float>> &a, std::vector<float> &b) {
    int n = 0;
    std::cout << "Количество уравнений: ";
    std::cin >> n;
    a.assign(n, std::vector<float>(n));
    b.resize(n);
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << "a[" << i << ',' << j << "] = ";
            std::cin >> a[i][j];
        }
        std::cout << "b[" << i << "] = ";
        std::cin >> b[i];
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

// метод Гаусса без выбора ведущего элемента
void solve(std::vector<std::vector<float>> &a,
           std::vector<float> &b, std::vector<float> &x) {
    int n = (int)a.size();
    x.resize(n);
    const float FLOAT_ZERO = 1e-07;  // машинный эквивалент нуля для типа float
    // прямой ход
    float l = 0;
    for(int k = 0; k < n - 1; ++k) {
        for(int i = k + 1; i < n; ++i) {
            if(std::fabs(a[k][k]) < FLOAT_ZERO)
                throw std::logic_error("Деление на ноль");
            l = a[i][k] / a[k][k];
            for(int j = k + 1; j < n; ++j) {
                a[i][j] -= l * a[k][j];
            }
            b[i] -= l * b[k];
        }
    }
    // обратный ход
    float s;
    for(int i = n - 1; i >= 0; --i) {
        s = 0;
        for (int j = n - 1; j > i; --j) {
            s += a[i][j] * x[j];
        }
        x[i] = (b[i] - s) / a[i][i];
    }
}

// метод Гаусса с выбором ведущего элемента по строке
void solve_with_item_selection(std::vector<std::vector<float>> &a,
           std::vector<float> &b, std::vector<float> &x) {
    int n = (int) a.size();
    x.resize(n);
    const float FLOAT_ZERO = 1e-07;  // машинный эквивалент нуля для типа float
    std::vector<std::pair<int, int>> changed_indexes;  // пары индексов, элементы которых меняются местами
    float l = 0, max = 0, tmp = 0;
    int max_index = 0;
    // прямой ход
    for (int k = 0; k < n - 1; ++k) {
            // выбор максимального элемента в строке
            max = std::fabs(a[k][k]);
            max_index = k;
            for (int j = k + 1; j < n; ++j) {
                if (std::fabs(a[k][j]) > max) {
                    max = a[k][j];
                    max_index = j;
                }
            }
            if(std::fabs(max) < FLOAT_ZERO)
                throw std::logic_error("Неопределенная система");
            // перестановка переменных
            if(max != a[k][k]) {
                changed_indexes.emplace_back(k, max_index);  // сохраняет последовательность перестановок
                for (int i = 0; i < n; ++i) {
                    tmp = a[i][k];
                    a[i][k] = a[i][max_index];
                    a[i][max_index] = tmp;
                }
            }
        for (int i = k + 1; i < n; ++i) {
            if(std::fabs(a[k][k]) < FLOAT_ZERO)
                throw std::logic_error("Деление на ноль");
            l = a[i][k] / a[k][k];
            for (int j = k + 1; j < n; ++j) {
                a[i][j] -= l * a[k][j];
            }
            b[i] -= l * b[k];
        }
    }
    // обратный ход
    float s;
    for(int i = n - 1; i >= 0; --i) {
        s = 0;
        for (int j = n - 1; j > i; --j) {
            s += a[i][j] * x[j];
        }
        x[i] = (b[i] - s) / a[i][i];
    }
    // восстановление правильного порядка неизвестных
    for(int i = (int)changed_indexes.size(); i >=0; --i) {
        tmp = x[changed_indexes[i].first];
        x[changed_indexes[i].first] = x[changed_indexes[i].second];
        x[changed_indexes[i].second] = tmp;
    }
}

// генерирует СЛАУ размерности n
void generate_system(std::vector<std::vector<float>> &a, std::vector<float> &f,
                     std::vector<float> &x, int n) {
    a.assign(n, std::vector<float>(n));
    f.resize(n);
    x.resize(n);
    srand ( time(nullptr) );
    const int MIN_RAND = -100;
    const int MAX_RAND = 100;
    // генерация матрицы А и вектора Х
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = MIN_RAND + rand() % (MAX_RAND - MIN_RAND + 1);
        }
        x[i] = MIN_RAND + rand() % (MAX_RAND - MIN_RAND + 1);
    }
    // вычисление столбца f
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            f[i] += a[i][j] * x[j];
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
    std::vector<std::vector<float>> a1;
    std::vector<float> b1;
    std::vector<float> x;
    std::vector<float> x1;
    std::vector<float> x2;

    generate_system(a1, b1, x, 10);
    std::vector<std::vector<float>> a2 = a1;
    std::vector<float> b2 = b1;

    std::cout << "Расширенная матрица системы:\n";
    show_system(a1, b1);
    solve(a1, b1, x1);

    std::cout << "\nВектор решений системы:\n";
    for(auto &i : x) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;

    std::cout << "\nВектор решений системы методом Гаусса "
                 "без выбора главного элемента:\n";
    for(auto &i : x1) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;

    solve_with_item_selection(a2, b2, x2);
    std::cout << "\nВектор решений системы методом Гаусса "
                 "с выбором главного элемента по строке:\n";
    for(float i : x2) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;

    float error1 = calculate_error(x1, x);
    float error2 = calculate_error(x2, x);
    std::cout << "\nПогрешность метода Гаусса без выбора главного элемента: "
        << error1
        << "\nПогрешность метода Гаусса с выбором главного элемента по строке: "
        << error2;

    return 0;
}
