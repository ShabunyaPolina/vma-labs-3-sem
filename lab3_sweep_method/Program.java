// Шабуня Полина (1 гр)
// лыб 3 метод прогонки (левая)
// 10 номер в списке подгруппы

import java.util.Arrays;

public class Program {

    static final int m = 10;  // номер в списке группы
    static final int p = (m % 4) + 1;
    static final int q = m % 3;
    static final int n = 10; // размер матрицы

    public static void main(String[] args) {
        float[][] A = generateMatrix(n);
        System.out.println("Матрица A:");
        showMatrix(A);

        float[] x0 = generateX();

        float[] b = calculateB(A, x0);
        System.out.println("Вектор B:\n" + Arrays.toString(b));

        System.out.println("\nРасширенная матрица:");
        showExtendedMatrix(A, b);

        System.out.println("\nПроверка устойчивости метода: " + isStable(A));

        float[][] newA = new float[A.length][A[0].length];
        float[] newB = new float[b.length];

        float[] x = solve(A, newA, b, newB);

        System.out.println("\nПреобразованная матрица:");
        showExtendedMatrix(newA, newB);

        System.out.println("\nПравильный вектор решений:\n" + Arrays.toString(x0));
        System.out.println("\nВектор решений, полученный методом левой прогонки:\n" + Arrays.toString(x));

        System.out.println("\nОтносительная погрешность метода:  " + calculateError(x, x0));
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

            if (cmp > 0)
                isThereStrict = true;
            else if (cmp < 0)
                return false;
        }
        return isThereStrict;
    }

    // реализация метода прогонки
    public static float[] solve(float[][] a, float[][] newA, float[] b, float[] newB) {
        float[] x = new float[b.length];
        int n = b.length;
        float tmp = 0;

        // прямой ход
        for(int i = 0; i < n; ++i) {
            newA[i][i] = 1;
        }
        tmp = a[n - 1][n - 1];
        newA[n-1][n-2] = a[n - 1][n - 2] / tmp;
        newB[n-1] = b[n - 1] / tmp;
        for (int i = n - 2; i > 0; --i) {
            tmp = a[i][i] - a[i][i+1] * newA[i+1][i];
            newA[i][i - 1] = a[i][i-1] / tmp;
            newB[i] = (b[i] - a[i][i+1] * newB[i+1]) / tmp;
        }
        newB[0] = (b[0] - a[0][1] * newB[1]) / (a[0][0] - a[0][1] * newA[1][0]);

        // обратный ход
        x[0] = newB[0];
        for(int i = 1; i < n; ++i) {
            x[i] = -newA[i][i-1] * x[i-1] + newB[i];
        }
        return x;
    }

    // вычисляет относительную погрешность в кубической норме
    public static float calculateError(float[] x, float[] trueX) {
        float tmp = 0, maxDelta = 0, maxTrueX = 0;
        for(int i = 0; i < x.length; ++i) {
            tmp = Math.abs(x[i] - trueX[i]);
            if(tmp > maxDelta)
                maxDelta = tmp;
            if(Math.abs(x[i]) > maxTrueX)
                maxTrueX = Math.abs(x[i]);
        }
        return maxDelta / maxTrueX;
    }
}
