// Шабуня Полина (1 гр)
// лаб 6 - степенной метод

import java.util.Arrays;
import java.util.Random;

public class Program {

    static final int m = 10;  // номер в списке группы
    static final int n = 10;  // размер матриц

    static double[] u0 = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    public static void main(String[] args) {

        double[][] A = generateMatrix(n);
        System.out.println("Матрица А:");
        showMatrix(A);

        System.out.println("\nПЕРВОЕ СОБСТВЕННОЕ ЗНАЧЕНИЕ:");
        double firstEigenvalue1 = 0;
        double firstEigenvalue2 = 0;
        double[] uk = new double[0];
        System.out.println("Приближенное значение наибольшего " +
                "по модулю собственного значения на ");
        for (int k = 46; k <= 50; ++k) {
            System.out.println(k + "-ой итерации:  ");
            firstEigenvalue1 = calculateFirstEigenvalue1(A, u0, k);
            System.out.print("(1 сп.) " + firstEigenvalue1);
            firstEigenvalue2 = calculateFirstEigenvalue2(A, u0, k);
            System.out.println("  (2 сп.) " + firstEigenvalue2);
            uk = calculateU(calcFinalV(A, u0, k));
            System.out.println("Собственный вектор:  " + Arrays.toString(uk));
        }

        double[] vk = calculateV(A, uk);  // V^(k+1)
        double[] vectorToCheck1 = calculateVectorToCheck(vk, uk, firstEigenvalue1);
        double[] vectorToCheck2 = calculateVectorToCheck(vk, uk, firstEigenvalue2);
        System.out.println("\nВекторы для проверки: ");
        System.out.println("(1) " + Arrays.toString(vectorToCheck1));
        System.out.println("Кубическая норма:  " + calcVectCubicNorm(vectorToCheck1));
        System.out.println("(2) " + Arrays.toString(vectorToCheck2));
        System.out.println("Кубическая норма:  " + calcVectCubicNorm(vectorToCheck2));

        System.out.println("\nВТОРОЕ СОБСТВЕННОЕ ЗНАЧЕНИЕ:");
        System.out.println("Приближенное значение второго по величине модуля собственное значение:");
        System.out.println("m=30 (1):");
        double secondEigenvalue = calculateSecondEigenvalue(A, u0, firstEigenvalue1, 30);
        System.out.println(secondEigenvalue);
        System.out.println("m=50 (1):");
        secondEigenvalue = calculateSecondEigenvalue(A, u0, firstEigenvalue1, 50);
        System.out.println(secondEigenvalue);
        System.out.println("m=50 (2):");
        secondEigenvalue = calculateSecondEigenvalue(A, u0, firstEigenvalue2, 50);
        System.out.println(secondEigenvalue);

        double[] v1 = calculateVectorToCheck(vk, uk, firstEigenvalue1);
        double[] v2 = calculateVectorToCheck(vk, uk, firstEigenvalue2);
        System.out.println("\nСобственные векторы:  " + "\n(1) " +
                Arrays.toString(v1) +
                "\n(2) " + Arrays.toString(v2));

        double[] vec1 = calculateVectorToCheck2(A, v1, secondEigenvalue);
        double[] vec2 = calculateVectorToCheck2(A, v2, secondEigenvalue);
        double[] vec3 = calculateVectorToCheck2(A, v2, secondEigenvalue);
        System.out.println("\nВекторы для проверки: ");
        System.out.println("(1) " + Arrays.toString(vec1));
        System.out.println("Кубическая норма:  " + calcVectCubicNorm(vec1));
        System.out.println("(2) " + Arrays.toString(vec2));
        System.out.println("Кубическая норма:  " + calcVectCubicNorm(vec2));
    }

    // вывод матрицы
    private static void showMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double element : row) {
                System.out.print(element + "  ");
            }
            System.out.println();
        }
    }

    // генерирует симметрическую матрицу А размерности n*n
    private static double[][] generateMatrix(int n) {
        double[][] matrix = new double[n][n];
        final int MIN_RAND = -4;
        final int MAX_RAND = 0;
        Random random = new Random();
        // генерация недиагональных элементов матрицы А
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                matrix[i][j] = matrix[j][i] = MIN_RAND +
                        random.nextInt(MAX_RAND - MIN_RAND + 1);
            }
        }
        // вычисление диагональных элементов матрицы А
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                matrix[i][i] -= matrix[i][j];
            }
        }
        return matrix;
    }

    // перемножает матрицу и вектор
    private static double[] multiply(double[][] matrix, double[] vector) {
        int n = vector.length;
        double[] res = new double[n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res[i] += matrix[i][j] * vector[j];
            }
        }
        return res;
    }

    // умножает вектор на число
    private static double[] multiply(double[] vector, double number) {
        double[] res = new double[vector.length];
        for (int i = 0; i < vector.length; ++i) {
            res[i] = vector[i] * number;
        }
        return res;
    }

    // делит вектор на число
    private static double[] divide(double[] vector, double number) {
        double[] res = new double[vector.length];
        for (int i = 0; i < vector.length; ++i) {
            res[i] = vector[i] / number;
        }
        return res;
    }

    // вычисляет разность векторов
    private static double[] subtract(double[] a, double[] b) {
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            res[i] = a[i] - b[i];
        }
        return res;
    }

    // вычисляет скалярное произведение векторов
    private static double calculateScalarProduct(double[] a, double[] b) {
        double res = 0;
        for (int i = 0; i < a.length; ++i) {
            res += a[i] * b[i];
        }
        return res;
    }

    // вычисляет кубическую норму вектора
    private static double calcVectCubicNorm(double[] vector) {
        double max = 0;
        for (double element : vector) {
            if (element > max)
                max = element;
        }
        return max;
    }

    // вычисляет вектор V на следующей итерации
    private static double[] calculateV(double[][] A, double[] u) {
        return multiply(A, u);
    }

    // вычисляет вектор U на следующей итерации
    private static double[] calculateU(double[] v) {
        return divide(v, calcVectCubicNorm(v));
    }

    // вычисляет V^(k) (V на k-ой итерации)
    private static double[] calcFinalV(double[][] A, double[] u0, int k) {
        double[] u = u0;
        while (k > 1) {
            u = calculateU(calculateV(A, u));
            --k;
        }
        return calculateV(A, u);
    }

    // возвращает индекс наибольшей по модулю компоненты вектора
    private static int chooseComponent(double[] vk) {
        int maxIndex = 0;
        for (int i = 1; i < vk.length; ++i) {
            if (Math.abs(vk[i]) > vk[maxIndex])
                maxIndex = i;
        }
        return maxIndex;
    }

    // вычисляет наибольшее по модулю собственное значение матрицы А (1 способ)
    private static double calculateFirstEigenvalue1(double[][] A, double[] u0, int k) {
        double[] vk = calcFinalV(A, u0, k);  // V^(k)
        int i = chooseComponent(vk);
        double[] uk = calculateU(vk);  // U^(k)
        vk = calculateV(A, uk);  // V^(k+1)
        return vk[i] * Math.signum(uk[i]);
    }

    // вычисляет наибольшее по модулю собственное значение матрицы А (2 способ)
    private static double calculateFirstEigenvalue2(double[][] A, double[] u0, int k) {
        double[] vk = calcFinalV(A, u0, k);  // V^(k)
        double[] uk = calculateU(vk);  // U^(k)
        vk = calculateV(A, uk);  // V^(k+1)
        return calculateScalarProduct(vk, uk) / calculateScalarProduct(uk, uk);
    }

    // вычисляет вектор V^(k+1)-lambda1*U^(k)
    private static double[] calculateVectorToCheck(double[] vk, double[] uk, double lambda) {
        return subtract(vk, multiply(uk, lambda));
    }

    // вычисляет второе по величине молудя СЗ матрицы А
    private static double calculateSecondEigenvalue(double[][] A, double[] u0,
                                                    double firstEigenvalue, int m) {
        double[] vm = calcFinalV(A, u0, m - 1);  // V^(m-1)
        double[] um = calculateU(vm);  // U^(m-1)
        int i = chooseComponent(calculateVectorToCheck(vm, um, firstEigenvalue));
        vm = calculateV(A, um);  // V^(m)
        double[] vm1 = calculateV(A, calculateU(vm));  // V^(m+1)
        return ((vm1[i] * calcVectCubicNorm(vm) - firstEigenvalue * vm[i])
                / (vm[i] - firstEigenvalue * um[i]));
    }

    private static double[] calculateVectorToCheck2(double[][] A, double[] vec,
                                                    double secondEigenvalue) {
        return subtract(multiply(A, vec), multiply(vec, secondEigenvalue));
    }
}