// Шабуня Полина (1 гр)
// лыб 3 метод прогонки (левая)
// 10 номер в списке подгруппы

import java.util.Arrays;

public class Program {

    static final int m = 10;
    static final int p = (m % 4) + 1;
    static final int q = m % 3;
    static final int n = 10;

    public static void main(String[] args) {
        System.out.println(m + " " + p + " " + q);

        float[][] A = generateMatrix(n);
        showMatrix(A);

        float[] x = generateX(n);
        System.out.println(Arrays.toString(x));

        float[] b = calculateB(A, x, m);
        System.out.println(Arrays.toString(b));

        System.out.println(isStable(A));
    }

    // вывод мвтрицы
    public static void showMatrix(float[][] matrix) {
        int n = matrix.length;
        for (float[] floats : matrix) {
            for (int j = 0; j < n; ++j) {
                System.out.print(floats[j] + "  ");
            }
            System.out.println();
        }
    }

    // генерирует трехдиагональную матрицу А размерности n*n
    public static float[][] generateMatrix(int n) {
        float[][] matrix = new float[n][n];
        for (int i = 0; i < n; ++i) {
            matrix[i][i] = (float) (5 * Math.pow(i + 1, p / 2.));
            if (i != 0)
                matrix[i][i - 1] = (float) (Math.pow(i + 1, p / 2.) + Math.pow(i, q / 2.));
            if (i != n - 1)
                matrix[i][i + 1] = (float) (Math.pow(i + 1, p / 2.) + Math.pow(i + 2, q / 2.));
        }
        return matrix;
    }

    // генерирует вектор x
    public static float[] generateX(int n) {
        float[] x = new float[n];
        for (int i = 0; i < n; ++i) {
            x[i] = m + i;
        }
        return x;
    }

    // вычисляет столбец b СЛАУ
    public static float[] calculateB(float[][] a, float[] x, int m) {
        if (a.length != x.length)
            return null;
        for (float[] elem : a) {
            if (elem.length != x.length)
                return null;
        }

        float[] b = new float[x.length];
        for (int i = 0; i < x.length; ++i) {
            for (int j = 0; j < x.length; ++j) {
                b[i] += a[i][j] * x[j];
            }
        }
        return b;
    }

    // проверка устойчивости метода
    public static boolean isStable(float[][] a) {
        boolean isThereStrict = false;
        int cmp = 0;
        for (int i = 0; i < n; ++i) {
            if (i == 0)
                cmp = Float.compare(Math.abs(a[i][i]), Math.abs(a[i][i + 1]));
            else if (i == n - 1)
                cmp = Float.compare(Math.abs(a[i][i]), Math.abs(a[i][i - 1]));
            else
                cmp = Float.compare(Math.abs(a[i][i]), Math.abs(a[i][i - 1]) + Math.abs(a[i][i + 1]));

            if(cmp > 0)
                isThereStrict =  true;
            else if(cmp < 0)
                return false;
        }
        return isThereStrict;
    }

    // реализация метода прогонки
    public static float[] solve(float[][] a, float[] newA, float[] b, float[] newB) {
        float[] x = new float[b.length];
        float[] alfa = new float[b.length];
        float[] beta = new float[b.length];
        int n = b.length;
        float tmp = 0;

        // прямой ход
        tmp = a[n-1][n-1];
        alfa[n-1] = -a[n -1][n - 2] / tmp;
        beta[n-1] = b[n - 1] / tmp;

        for(int i = n-2; i > 0; --i) {
            tmp = 
        }
    }
}
