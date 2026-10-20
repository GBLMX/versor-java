package vsr;

import java.util.Random;

/** Output-driven dense GP variants: fewer result writes, split accumulators, FMA. */
public final class GpKernelBench {
    public static void main(String[] args) {
        Random rnd = new Random(1234);
        for (int n = 8; n <= 10; n++) {
            Signature sig = Signature.euclidean(n);
            int size = sig.size;
            double[] a = new double[size];
            double[] b = new double[size];
            for (int i = 0; i < size; i++) {
                a[i] = rnd.nextDouble() * 2 - 1;
                b[i] = rnd.nextDouble() * 2 - 1;
            }
            Multivector ma = new Multivector(sig, a);
            Multivector mb = new Multivector(sig, b);
            int[][] sign = sig.signTable();
            Multivector ref = ma.gp(mb);
            double base = GaBench.ns(() -> {
                Multivector r = ma.gp(mb);
                return r.get(0) + r.get(size - 1);
            });
            double out1 = GaBench.ns(() -> {
                double[] r = gpOut(a, b, sign, size, 1, false);
                return r[0] + r[size - 1];
            });
            double out4 = GaBench.ns(() -> {
                double[] r = gpOut(a, b, sign, size, 4, false);
                return r[0] + r[size - 1];
            });
            double out4f = GaBench.ns(() -> {
                double[] r = gpOut(a, b, sign, size, 4, true);
                return r[0] + r[size - 1];
            });
            double err = 0;
            for (int v = 0; v < 3; v++) {
                double[] r = gpOut(a, b, sign, size, v == 2 ? 4 : 1, v == 2);
                for (int i = 0; i < size; i++) {
                    err = Math.max(err, Math.abs(r[i] - ref.get(i)));
                }
            }
            System.out.printf("GPK,%d,base,%.1f%n", n, base);
            System.out.printf("GPK,%d,out1,%.1f%n", n, out1);
            System.out.printf("GPK,%d,out4,%.1f%n", n, out4);
            System.out.printf("GPK,%d,out4fma,%.1f%n", n, out4f);
            System.out.printf("GPK,%d,err,%.2e%n", n, err);
        }
    }

    private static double[] gpOut(double[] a, double[] b, int[][] sign, int size, int unroll,
                                  boolean fma) {
        double[] r = new double[size];
        for (int k = 0; k < size; k++) {
            double s0 = 0;
            double s1 = 0;
            double s2 = 0;
            double s3 = 0;
            int x = 0;
            if (unroll == 4) {
                for (; x + 3 < size; x += 4) {
                    int b0 = x ^ k;
                    int b1 = (x + 1) ^ k;
                    int b2 = (x + 2) ^ k;
                    int b3 = (x + 3) ^ k;
                    double a0 = a[x];
                    double a1 = a[x + 1];
                    double a2 = a[x + 2];
                    double a3 = a[x + 3];
                    if (a0 != 0 && b[b0] != 0) {
                        s0 = fma ? Math.fma(sign[x][b0] * a0, b[b0], s0)
                                : s0 + sign[x][b0] * a0 * b[b0];
                    }
                    if (a1 != 0 && b[b1] != 0) {
                        s1 = fma ? Math.fma(sign[x + 1][b1] * a1, b[b1], s1)
                                : s1 + sign[x + 1][b1] * a1 * b[b1];
                    }
                    if (a2 != 0 && b[b2] != 0) {
                        s2 = fma ? Math.fma(sign[x + 2][b2] * a2, b[b2], s2)
                                : s2 + sign[x + 2][b2] * a2 * b[b2];
                    }
                    if (a3 != 0 && b[b3] != 0) {
                        s3 = fma ? Math.fma(sign[x + 3][b3] * a3, b[b3], s3)
                                : s3 + sign[x + 3][b3] * a3 * b[b3];
                    }
                }
            }
            for (; x < size; x++) {
                int bx = x ^ k;
                double ax = a[x];
                if (ax != 0 && b[bx] != 0) {
                    s0 = fma ? Math.fma(sign[x][bx] * ax, b[bx], s0)
                            : s0 + sign[x][bx] * ax * b[bx];
                }
            }
            r[k] = (s0 + s1) + (s2 + s3);
        }
        return r;
    }
}
