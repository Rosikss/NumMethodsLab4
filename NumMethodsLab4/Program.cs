using System;

class Program
{
    static double FLagrange(double x)
    {
        return 2 * Math.Pow(x, 7) + 3 * Math.Pow(x, 6) + 3 * Math.Pow(x, 4) - 3;
    }

    static double LagrangeInterpolation(double x, double[] nodes, double[] values)
    {
        double result = 0;
        int n = nodes.Length;

        for (int i = 0; i < n; i++)
        {
            double term = values[i];
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    term *= (x - nodes[j]) / (nodes[i] - nodes[j]);
                }
            }
            result += term;
        }

        return result;
    }
    //////////////Lagrange
    static double FSquare(double x)
    {
        return 2 * Math.Pow(x, 7) + 3 * Math.Pow(x, 6) + 3 * Math.Pow(x, 4) - 3;
    }

    static double DF(double x)
    {
        return 14 * Math.Pow(x, 6) + 18 * Math.Pow(x, 5) + 12 * Math.Pow(x, 3);
    }

    static void Gauss(double[,] m, double[] ans, int N)
    {
        double q, S;
        for (int k = 1; k <= N; ++k)
        {
            q = m[k, k];
            for (int i = 1; i <= N + 1; ++i) m[k, i] /= q;
            for (int i = k + 1; i <= N; ++i)
            {
                q = m[i, k];
                for (int j = 1; j <= N + 1; ++j) m[i, j] -= m[k, j] * q;
            }
        }
        ans[N] = m[N, N + 1];
        for (int i = N - 1; i >= 1; --i)
        {
            S = 0;
            for (int j = i + 1; j <= N; ++j) S += m[i, j] * ans[j];
            ans[i] = m[i, N + 1] - S;
        }
    }

    static void SquareSpline(double[] ans)
    {
        double[,] m = new double[101, 101];
        m[1, 1] = 4;
        m[1, 2] = 2;
        m[1, 6] = FSquare(5) - FSquare(3);
        m[2, 1] = 4;
        m[2, 2] = 1;
        m[2, 3] = -4;
        m[2, 4] = -1;
        m[3, 3] = 4;
        m[3, 4] = 2;
        m[3, 5] = 1;
        m[3, 6] = FSquare(5);
        m[4, 3] = 8;
        m[4, 4] = 6;
        m[4, 6] = DF(7);
        m[5, 3] = 16;
        m[5, 4] = 4;
        m[5, 5] = 1;
        m[5, 6] = FSquare(7);
        Gauss(m, ans, 5);
    }

    /////////// square spline
    static List<Tuple<double, double>> splineRes = new List<Tuple<double, double>>();

    static double FLinear(double x)
    {
        return 2 * Math.Pow(x, 7) + 3 * Math.Pow(x, 6) + 3 * Math.Pow(x, 4) - 3;
    }

    static void LinearSpline(double start, double step, double finish)
    {
        double prev = start;
        double next = start + step;
        double a, b;
        while (next <= finish)
        {
            a = (FLinear(next) - FLinear(prev)) / (next - prev);
            b = FLinear(prev) - a * prev;
            splineRes.Add(Tuple.Create(a, b));
            prev = next;
            next += step;
        }
    }
    static void Main()
    {
        double[] nodes = { 3, 4, 5, 6, 7 };
        double[] values = new double[nodes.Length];

        for (int i = 0; i < nodes.Length; i++)
        {
            values[i] = FLagrange(nodes[i]);
        }

        Console.WriteLine("Interpolation polynomial of Lagrange for function 2x^7 + 3x^6 + 3x^4 - 3 on interval [3..7]:");
        for (double x = 3; x <= 7; x += 0.5)
        {
            Console.WriteLine($"F({x}) = {LagrangeInterpolation(x, nodes, values)}");
        }
        ///////Lagrange
        Console.WriteLine("This program finds the square spline for 2x^7 + 3x^6 + 3x^4 - 3");
        Console.WriteLine("The points are 3, 5, 7");
        double[] ans = new double[101];
        SquareSpline(ans);
        Console.WriteLine("The square spline polynomials are:");
        Console.WriteLine($"{ans[1]}x^2 + {ans[2]}x + {FSquare(3)}");
        Console.WriteLine($"{ans[3]}x^2 + {ans[4]}x + {ans[5]}");
        //////////Square spline
        Console.WriteLine("This program finds the linear spline for 2x^7 + 3x^6 + 3x^4 - 3");
        Console.WriteLine("Please enter the borders of your segment and the step size");
        double start = double.Parse(Console.ReadLine());
        double finish = double.Parse(Console.ReadLine());
        double step = double.Parse(Console.ReadLine());

        Console.WriteLine("The functions for linear spline are:");
        LinearSpline(start, step, finish);
        foreach (var res in splineRes)
        {
            Console.Write($"{res.Item1}x");
            Console.WriteLine(res.Item2 >= 0 ? $"+{res.Item2}" : $"{res.Item2}");
        }
    }
}