// Шабуня Полина (1 гр)
// лаб2 LDLt разложение

#include <iostream>
#include <vector>


// ввод матрицы
void fill_matrix(std::vector<std::vector<float>> &a) {
    int n = 0;
    std::cout << "n: ";
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

void decompose(std::vector<std::vector<float>> &a) {
    int n = (int)a.size();
    std::vector<float> tmp(n*n);
    for(int k = 0; k < n - 1; ++k) {
        for(int i = k + 1; i < n; ++i) {
            tmp[i] = a[i][k];
            a[i][k] /= a[k][k];
            for(int j = k + 1; j <= i; ++j) {
                a[i][j] -= a[i][k] * tmp[j];
            }
        }
    }
}

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

void toD(std::vector<std::vector<float>> &a, std::vector<std::vector<float>> &d) {
    int n = (int)a.size();
    d.assign(n, std::vector<float>(n));
    for(int i = 0; i < n; ++i) {
        d[i][i] = a[i][i];
    }
}

void transpose(std::vector<std::vector<float>> &x) {
    int n = (int)x.size();
    for(int i = 0; i < n; ++i) {
        for(int j = i + 1; j < n; ++j) {
            std::swap(x[i][j], x[j][i]);
        }
    }
}

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

int main() {
    std::vector<std::vector<float>> a{{4,2,3},{2,4,1},{3,1,3}};
    std::vector<std::vector<float>> l;
    std::vector<std::vector<float>> lt;
    std::vector<std::vector<float>> d;
    std::vector<std::vector<float>> res1;
    std::vector<std::vector<float>> res2;
    //fill_matrix(a);
    decompose(a);
    toL(a,l);
    lt = l;
    show_matrix(a);
    std::cout << "\nL:\n";
    show_matrix(l);
    transpose(lt);
    std::cout << "\nLt:\n";
    show_matrix(lt);
    toD(a, d);
    std::cout << "\nD:\n";
    show_matrix(d);
    std::cout << "\n";
    multiply(l,d, res1);
    multiply(res1, lt, res2);
    show_matrix(res2);
    return 0;
}
