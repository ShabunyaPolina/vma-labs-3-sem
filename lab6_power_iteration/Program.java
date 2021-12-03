// Шабуня Полина (1 гр)
// лаб 6 - степенной метод

import java.util.Arrays;
import java.util.Random;

public class Program {

    static final int m = 10;  // номер в списке группы
    static final int n = 10;  // размер матриц

    static double[] u0 = {1, 0, 0};

    public static void main(String[] args) {
        double[][] A = {
                {1, 2, 3},
                {2, 1, 2},
                {3, 2, 1}};

        showMatrix(A);

        System.out.println();
        double[] v1 = calculateV(A, u0);
        System.out.println(Arrays.toString(v1));
        System.out.println();
        System.out.println(calcVectCubicNorm(v1));
        System.out.println();
        double[] u1 = calculateU(v1);
        System.out.println(Arrays.toString(u1));
        System.out.println();

        double[] v2 = calculateV(A, u1);
        System.out.println(Arrays.toString(v2));
        System.out.println();
        System.out.println(calcVectCubicNorm(v2));
        System.out.println();
        double[] u2 = calculateU(v2);
        System.out.println(Arrays.toString(u2));
        System.out.println();

        double[] finalV = calcFinalV(A, u0, 4);
        System.out.println(Arrays.toString(finalV));
        double[] finalU = calculateU(finalV);
        System.out.println(Arrays.toString(finalU));
        double[] finalfinalV = calculateV(A, finalU);
        System.out.println(Arrays.toString(finalfinalV));
        System.out.println();

        System.out.println(calculateFirstEigenvalue1(A, u0, 4));
        System.out.println(calculateFirstEigenvalue2(A, u0, 4));

        System.out.println("Приближенный собственный вектор u^(k): " +
                Arrays.toString(finalU));
        System.out.println();

        System.out.println(Arrays.toString(calculateVectorToCheck(finalfinalV,
                finalU, calculateFirstEigenvalue1(A, u0, 4))));
        System.out.println(Arrays.toString(calculateVectorToCheck(finalfinalV,
                finalU, calculateFirstEigenvalue2(A, u0, 4))));

        System.out.println(calcVectCubicNorm(calculateVectorToCheck(finalfinalV,
                finalU, calculateFirstEigenvalue1(A, u0, 4))));
        System.out.println(calcVectCubicNorm(calculateVectorToCheck(finalfinalV,
                finalU, calculateFirstEigenvalue2(A, u0, 4))));
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
        int i = chooseComponent(vk);
        double[] uk = calculateU(vk);  // U^(k)
        vk = calculateV(A, uk);  // V^(k+1)
        return calculateScalarProduct(vk, uk) / calculateScalarProduct(uk, uk);
    }

    // выбирает номер компоненты векторов V и U для вычисления первого СЗ (max|V^(k)i|)
    private static int chooseComponent(double[] vk) {
        int maxIndex = 0;
        for (int i = 1; i < vk.length; ++i) {
            if (vk[i] > vk[maxIndex])
                maxIndex = i;
        }
        return maxIndex;
    }

    // вычисляет скалярное произведение векторов
    private static double calculateScalarProduct(double[] a, double[] b) {
        double res = 0;
        for(int i = 0; i < a.length; ++i) {
            res += a[i] * b[i];
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

    // умножает вектор на число
    private static double[] multiply(double[] vector, double number) {
        double[] res = new double[vector.length];
        for (int i = 0; i < vector.length; ++i) {
            res[i] = vector[i] * number;
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

    // вычисляет вектор V^(k+1)-lambda1*U^(k)
    private static double[] calculateVectorToCheck(double[] vk, double[] uk, double lambda) {
        return subtract(vk, multiply(uk, lambda));
    }

    // вычисляет разность векторов
    private static double[] subtract(double[] a, double[] b) {
        double[] res = new double[a.length];
        for(int i = 0; i < a.length; ++i) {
            res[i] = a[i] - b[i];
        }
        return res;
    }

    // вычисляет второе по величине молудя СЗ матрицы А
    private static double calculateSecondEigenvalue(double[][] A, double[] u0, int k) {
//        double[] vk = calcFinalV(A, u0, k);  // V^(k)
//        int i = chooseComponent(vk);
//        double[] uk = calculateU(vk);  // U^(k)
//        vk = calculateV(A, uk);  // V^(k+1)
//        return vk[i] * Math.signum(uk[i]);
        return 1.1;
    }
}