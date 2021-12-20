// Шабуня Полина (1 гр)
// лыб 7 Нахождение корней собственного многочлена

public class Program {

    static final double e = 0.0001;  // точность вычисления корней
    // интервал, в котором находится корень многочлена
    // (получен в результате "ручного" отделения корней)
    static final double l = -22047;
    static final double r = -0.9993;

    static double x0;  // начальное приближение

    public static void main(String[] args) {

        double[][] F = {
                {-8, 244, -22047, 726036},
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0}};

        System.out.println("Каноническая форма Фробениуса матрицы из ЛР5:");
        showMatrix(F);
        double[] p = {-F[0][3], -F[0][2], -F[0][1], -F[0][0], 1};
        System.out.println("Собственный многочлен:");
        System.out.print("P(L) = ");
        System.out.println(polinomToString(p));
        System.out.println("\nПроизводная:");
        System.out.print("P`(L) = ");
        double[] dp = calculateDerivative(p);
        System.out.println(polinomToString(dp));
        System.out.println("\nРешение уравнения P`(L) = 0 методом Ньютона:");
        System.out.println("После предварительного отделения корней производной многочлена P(L) получаем," +
                "\nчто корень производной находится в промежутке [" + l + "; " + r + "]");
        System.out.println("Выбор начального приближения:");
        System.out.print(l + " - ");
        System.out.println(isInitialApproximation(dp, l));
        System.out.print(r + " - ");
        System.out.println(isInitialApproximation(dp, r));
        System.out.println("Начальное приближение = " + r);
        x0 = r;
        System.out.println("Корень уравнения P`(L) = 0, полученный методом Ньютона:");
        double x = calculateRootByNewtonsMethod(x0, dp);
        System.out.println(x);
        double x01 = -23, x02 = -22;
        System.out.println("Корни уравнения P(L) = 0, полученные методом Ньютона:");
        double x1 = calculateRootByNewtonsMethod(x01, p);
        System.out.println("Первое собственное значение: " + x1);
        double x2 = calculateRootByNewtonsMethod(x02, p);
        System.out.println("Второе собственное значение: " + x2);
        System.out.println("Значение P(L) при первом собственном значении : "
                + calculatePolinom(p, x1));
        System.out.println("Значение P(L) при втором собственном значении : "
                + calculatePolinom(p, x2));
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

    // вывод многочлена
    private static String polinomToString(double[] p) {
        StringBuilder sb = new StringBuilder();
        for (int i = p.length - 1; i >= 1; --i) {
            sb.append(p[i]).append("*").append("L^").append(i).append(" + ");
        }
        sb.append(p[0]);
        return sb.toString();
    }

    // вычисляет производную
    private static double[] calculateDerivative(double[] p) {
        if (p.length <= 1)
            return new double[0];
        double[] derivative = new double[p.length - 1];
        for (int i = 1; i < p.length; ++i) {
            derivative[i - 1] = p[i] * i;
        }
        return derivative;
    }

    // проверяет условие выбора начального приближения
    private static boolean isInitialApproximation(double[] p, double x) {
        return calculatePolinom(p, x) *
                calculatePolinom(calculateDerivative(calculateDerivative(p)), x) > 0;
    }

    // вычисляется значение многочлена в точке х
    private static double calculatePolinom(double[] p, double x) {
        double res = 0;
        //System.out.println(Arrays.toString(p));
        for (int i = 0; i < p.length; ++i) {
            res += p[i] * Math.pow(x, i);
            //System.out.print(p[i] + "*" + x + "^" + i + " + ");
        }
        return res;
    }

    // совершает одну итерацию метода Ньютона
    private static double iterate(double xk, double[] p) {
        double[] dp = calculateDerivative(p);
        double res = xk;
        double tmp = 0;
        for (int i = 0; i < dp.length; ++i) {
            res += p[i] * Math.pow(xk, i);
            tmp += dp[i] * Math.pow(xk, i);
        }
        res += p[p.length - 1] * Math.pow(xk, p.length - 1);
        res /= tmp;
        res = -res + xk;
        return res;
    }

    // находит корень методом Ньютона
    private static double calculateRootByNewtonsMethod(double x0, double[] p) {
        double xk = x0, xkk = x0;
        do {
            xk = xkk;
            xkk = iterate(xk, p);
        } while (Math.abs(xkk - xk) >= e);
        return xkk;
    }
}
