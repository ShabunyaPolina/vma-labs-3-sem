// Шабуня Полина (1 гр)
// лыб 4 итерационные методы

import java.util.Arrays;

public class Program {

    static final int m = 10;  // номер в списке группы
    static final int p = (m % 4) + 1;
    static final int q = m % 3;
    static final int n = 10;  // размер матриц
    static final double e = 0.0001;
    static final int kMax = 1000;

    public static void main(String[] args) {

        float[][] A = generateMatrix(n);
        System.out.println("Матрица А:");
        showMatrix(A);

        float[] x = generateX();

        float[] b = calculateB(A, x);
        System.out.println("Вектор B:\n" + Arrays.toString(b));

        System.out.println("\nРасширенная матрица:");
        showExtendedMatrix(A, b);

        float[][] newA = new float[n][n];
        float[] newB = new float[n];
        convertToIterationView(A, b, newA, newB);
        System.out.println("\nВид матрицы, пригодный для итераций:");
        showExtendedMatrix(newA, newB);

        //--------------------------------------------------------

        System.out.println("\nЗАДАНИЕ 1");

        System.out.println("\nПроверка строгого диагонального преобладания: "
                + isDiagonalDominance(A));

        float[] x0 = new float[x.length];
        System.arraycopy(newB, 0, x0, 0, n);
        System.out.println("\nВектор начального приближения:\n"
                + Arrays.toString(x0));

        float[] X1 = new float[x.length];
        int k = solveSimpleIterationMethod(newA, newB, x0, X1);

        System.out.println("\nПравильный вектор решений:\n"
                + Arrays.toString(x));
        System.out.println("\nВектор решений, полученный методом " +
                "простых итераций:\n" + Arrays.toString(X1));

        System.out.println("\nКоличество итераций, необходимых для " +
                "достижения точности в кубичесткой норме " + e + ":  "
                + calcNumberOfIterations(newA, newB));

        System.out.println("\nКоличество совершенных итераций:  " + k);

        System.out.println("\nОтносительная погрешность метода:  "
                + calculateError(X1, x));

        //--------------------------------------------------------

        System.out.println("\nЗАДАНИЕ 2");

        float[] X2 = new float[x.length];

        int k2 = solveRelaxationMethod(newA, newB, x0, X2, 0.5F);
        System.out.println("\nВектор решений, полученный методом " +
                "релаксации при w = 0.5:\n" + Arrays.toString(X2));
        System.out.println("Количество совершенных итераций:  " + k2);
        System.out.println("Относительная погрешность:  "
                + calculateError(X2, x));

        System.arraycopy(newB, 0, x0, 0, n);
        int k3 = solveRelaxationMethod(newA, newB, x0, X2, 1F);
        System.out.println("\nВектор решений, полученный методом релаксации при w = 1:\n" + Arrays.toString(X2));
        System.out.println("Количество совершенных итераций:  " + k3);
        System.out.println("Относительная погрешность:  "
                + calculateError(X2, x));

        System.arraycopy(newB, 0, x0, 0, n);
        int k4 = solveRelaxationMethod(newA, newB, x0, X2, 1.5F);
        System.out.println("\nВектор решений, полученный методом релаксации при w = 1.5:\n" + Arrays.toString(X2));
        System.out.println("Количество совершенных итераций:  " + k4);
        System.out.println("Относительная погрешность:  "
                + calculateError(X2, x));
    }

    // вывод матрицы
    public static void showMatrix(float[][] matrix) {
        for (float[] row : matrix) {
            for (float element : row) {
                System.out.print(element + "  ");
            }
            System.out.println();
        }
    }

    // вывод расширенной матрицы
    public static void showExtendedMatrix(float[][] a, float[] b) {
        int n = a.length;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                System.out.print(a[i][j] + "  ");
            }
            System.out.print("   |   " + b[i] + '\n');
        }
    }

    // генерирует матрицу А размерности n*n
    public static float[][] generateMatrix(int n) {
        float[][] matrix = new float[n][n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    matrix[i][j] = (float) (5 * Math.pow(i + 1, p / 2.));
                else
                    matrix[i][j] = (float) (0.01 * (Math.pow(i + 1, p / 2.)
                            + Math.pow(j + 2, q / 2.)));
            }
        }
        return matrix;
    }

    // генерирует вектор x
    private static float[] generateX() {
        float[] x = new float[n];
        for (int i = 0; i < n; ++i) {
            x[i] = m + i;
        }
        return x;
    }

    // вычисляет столбец b СЛАУ
    public static float[] calculateB(float[][] a, float[] x) {
        float[] b = new float[x.length];
        for (int i = 0; i < x.length; ++i) {
            for (int j = 0; j < x.length; ++j) {
                b[i] += a[i][j] * x[j];
            }
        }
        return b;
    }

    // проверка строгого диагонального преобладания в матрице по строкам
    private static boolean isDiagonalDominance(float[][] a) {
        int n = a.length;
        boolean res = true;
        float sum;
        for (int i = 0; i < n; ++i) {
            sum = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j)
                    sum += Math.abs(a[i][j]);
            }
            if (Math.abs(a[i][i]) <= sum)
                res = false;
        }
        return res;
    }

    // приводит матрицу к виду, пригодному для итераций
    private static void convertToIterationView(float[][] a, float[] b,
                                               float[][] newA, float[] newB) {
        int n = a.length;
        for (int i = 0; i < n; ++i) {
            newA[i][i] = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j)
                    newA[i][j] = -a[i][j] / a[i][i];
            }
            newB[i] = b[i] / a[i][i];
        }
    }

    //  реализует метод простых итераций для решения СЛАУ
    public static int solveSimpleIterationMethod(float[][] a, float[] b, float[] xk, float[] res) {
        int n = a.length;
        int k = 0;
        System.arraycopy(xk, 0, res, 0, n);
        do {
            ++k;
            System.arraycopy(res, 0, xk, 0, n);
            for (int i = 0; i < n; ++i) {
                res[i] = 0;
                for (int j = 0; j < n; ++j) {
                    res[i] += xk[j] * a[i][j];
                }
                res[i] += b[i];
            }
        } while (isSolutionOver(xk, res));
        return k;
    }

    // проверяет условие выхода из итерационного цикла
    private static boolean isSolutionOver(float[] x0, float[] x1) {
        float[] tmp = new float[x0.length];
        for (int i = 0; i < x0.length; ++i) {
            tmp[i] = Math.abs(x1[i] - x0[i]);
        }
        return calcVectCubicNorm(tmp) > e;
    }

    // вычисляет количество итераций, необходимых для достижения точности е
    public static int calcNumberOfIterations(float[][] B, float[] x0) {
        return (int) Math.ceil((-(Math.log(e * (1 - calcMatrixCubicNorm(B)))) /
                Math.log(calcVectCubicNorm(x0))) - 1);
    }

    public static int solveRelaxationMethod(float[][] a, float[] b, float[] xk,
                                            float[] res, float w) {
        int n = a.length;
        System.arraycopy(xk, 0, res, 0, n);
        int k = 0;
        float tmp;
        do {
            ++k;
            System.arraycopy(res, 0, xk, 0, n);
            for (int i = 0; i < n; ++i) {
                tmp = 0;
                res[i] = (1 - w) * xk[i];
                for (int j = 0; j < i; ++j) {
                    tmp += res[j] * a[i][j];
                }
                for (int j = i + 1; j < n; ++j) {
                    tmp += xk[i] * a[i][j];
                }
                res[i] += w * (b[i] + tmp);
            }
        } while (isSolutionOver(xk, res) && k < kMax);
        if(k >= kMax) {
            System.out.println("Превышено допустимое максимальное " +
                    "количество итераций.");
        }
        return k;
    }

    // вычисляет кубическую норму вектора
    public static float calcVectCubicNorm(float[] v) {
        float max = 0;
        for (float element : v) {
            if (element > max)
                max = element;
        }
        return max;
    }

    // вычисляет кубическую норму матрицы
    public static float calcMatrixCubicNorm(float[][] m) {
        float sum = 0, maxSum = 0;
        for (float[] row : m) {
            sum = 0;
            for (float element : row) {
                sum += element;
            }
            if (sum > maxSum)
                maxSum = sum;
        }
        return maxSum;
    }

    // вычисляет относительную погрешность в кубической норме
    public static float calculateError(float[] x, float[] trueX) {
        float[] tmp = new float[x.length];
        for (int i = 0; i < x.length; ++i) {
            tmp[i] = Math.abs(x[i] - trueX[i]);
        }
        return calcVectCubicNorm(tmp) / calcVectCubicNorm(trueX);
    }
}