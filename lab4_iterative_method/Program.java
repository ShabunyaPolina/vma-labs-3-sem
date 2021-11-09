// Шабуня Полина (1 гр)
// лыб 4 итерационные методы

import java.util.Arrays;

public class Program {

    static final int m = 10;  // номер в списке группы
    static final int p = (m % 4) + 1;
    static final int q = m % 3;
    static final int n = 10;  // размер матриц
    static final double e = 0.0001;

    public static void main(String[] args) {
        float[][] A = generateMatrix(n);
        showMatrix(A);
        float[] x = generateX();
        System.out.println(Arrays.toString(x));
        float[] b = calculateB(A, x);
        System.out.println(Arrays.toString(b));

        showExtendedMatrix(A, b);

        System.out.println(isDiagonalDominance(A));

        float[][] newA = new float[n][n];
        float[] newB = new float[n];
        convertToIterationView(A, b, newA, newB);
        showExtendedMatrix(newA, newB);
        System.out.println(Arrays.toString(newB));

//        float[][] a = {{10, 1, 1}, {2, -10, 1}, {2, 2, 10}};
//        float[] b = {12, 13, 14};
//        showExtendedMatrix(a, b);
//        System.out.println(isDiagonalDominance(a));
//        float[][] newA = new float[3][3];
//        float[] newB = new float[3];
//        convertToIterationView(a, b, newA, newB);
//        showExtendedMatrix(newA, newB);

        System.out.println(Arrays.toString(solve(newA, newB, b)));
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
                            + Math.pow(i + 2, q / 2.)));
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

    public static float[] solve(float[][] a, float[] b, float[] xk) {
        int n = a.length;
        float[] nextX = new float[xk.length];
        int k = 0;
        while(isSolutionOver(xk, nextX)) {
            ++k;
            System.arraycopy(nextX, 0, xk, 0, n);
        for (int i = 0; i < n; ++i) {
            nextX[i] = 0;
            for (int j = 0; j < n; ++j) {
                nextX[i] += xk[j] * a[i][j];
               // System.out.print(xk[j] + " * " + a[i][j] + " + ");
            }
            nextX[i] += b[i];
            //System.out.println(b[i] + " = " + nextX[i]);
        }
           // System.out.println();
         }
        System.out.println(k);
        return nextX;
    }

    private static boolean isSolutionOver(float[] x0, float[] x1) {
        float tmp, max = 0;
        for (int i = 0; i < x0.length; ++i) {
            tmp = Math.abs(x1[i] - x0[i]);
            if (tmp > max)
                max = tmp;
        }
        return max > e;
    }

    // вычисляет относительную погрешность в кубической норме
    public static float calculateError(float[] x, float[] trueX) {
        float tmp = 0, maxDelta = 0, maxTrueX = 0;
        for (int i = 0; i < x.length; ++i) {
            tmp = Math.abs(x[i] - trueX[i]);
            if (tmp > maxDelta)
                maxDelta = tmp;
            if (Math.abs(x[i]) > maxTrueX)
                maxTrueX = Math.abs(x[i]);
        }
        return maxDelta / maxTrueX;
    }
}
