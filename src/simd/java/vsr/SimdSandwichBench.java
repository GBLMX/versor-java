package vsr;

import java.util.Random;

/** Scalar vs Vector-API comparison for the batch sandwich kernels. */
public final class SimdSandwichBench {
    private static final int N = 1 << 20;

    public static void main(String[] args) {
        Random rnd = new Random(42);
        double[] x = new double[N];
        double[] y = new double[N];
        double[] z = new double[N];
        double[] p7 = new double[N];
        double[] p11 = new double[N];
        double[] p13 = new double[N];
        double[] p14 = new double[N];
        for (int i = 0; i < N; i++) {
            x[i] = rnd.nextDouble() * 4 - 2;
            y[i] = rnd.nextDouble() * 4 - 2;
            z[i] = rnd.nextDouble() * 4 - 2;
            p7[i] = rnd.nextDouble() * 4 - 2;
            p11[i] = rnd.nextDouble() * 4 - 2;
            p13[i] = rnd.nextDouble() * 4 - 2;
            p14[i] = rnd.nextDouble() * 4 - 2;
        }
        double s = 0.8, b1 = 0.1, b2 = -0.2, b3 = 0.3;
        double b01 = 0.2, b02 = -0.1, b12 = 0.4, b03 = 0.3, b13 = -0.2, b23 = 0.1, ps = 0.05;

        double[] sx = new double[N], sy = new double[N], sz = new double[N];
        double[] vx = new double[N], vy = new double[N], vz = new double[N];
        scalarE3(s, b1, b2, b3, x, y, z, sx, sy, sz);
        SimdSandwich.e3(s, b1, b2, b3, x, y, z, vx, vy, vz);
        System.out.printf("VERIFY,E3,%.3e%n", maxDiff(sx, vx, sy, vy, sz, vz));

        double[] s7 = new double[N], s11 = new double[N], s13 = new double[N], s14 = new double[N];
        double[] t7 = new double[N], t11 = new double[N], t13 = new double[N], t14 = new double[N];
        scalarPga(s, b01, b02, b12, b03, b13, b23, ps,
                p7, p11, p13, p14, s7, s11, s13, s14);
        SimdSandwich.pga(s, b01, b02, b12, b03, b13, b23, ps,
                p7, p11, p13, p14, t7, t11, t13, t14);
        System.out.printf("VERIFY,PGA,%.3e%n", maxDiff(s7, t7, s11, t11, s13, t13, s14, t14));

        for (int i = 0; i < 5; i++) {
            scalarE3(s, b1, b2, b3, x, y, z, sx, sy, sz);
            SimdSandwich.e3(s, b1, b2, b3, x, y, z, vx, vy, vz);
            scalarPga(s, b01, b02, b12, b03, b13, b23, ps,
                    p7, p11, p13, p14, s7, s11, s13, s14);
            SimdSandwich.pga(s, b01, b02, b12, b03, b13, b23, ps,
                    p7, p11, p13, p14, t7, t11, t13, t14);
        }

        double scalarE3 = time(() -> scalarE3(s, b1, b2, b3, x, y, z, sx, sy, sz));
        double vectorE3 = time(() -> SimdSandwich.e3(s, b1, b2, b3, x, y, z, vx, vy, vz));
        double scalarPga = time(() -> scalarPga(s, b01, b02, b12, b03, b13, b23, ps,
                p7, p11, p13, p14, s7, s11, s13, s14));
        double vectorPga = time(() -> SimdSandwich.pga(s, b01, b02, b12, b03, b13, b23, ps,
                p7, p11, p13, p14, t7, t11, t13, t14));

        System.out.printf("SIMD,E3,scalar,%.3f%n", scalarE3);
        System.out.printf("SIMD,E3,vector,%.3f%n", vectorE3);
        System.out.printf("SIMD,PGA,scalar,%.3f%n", scalarPga);
        System.out.printf("SIMD,PGA,vector,%.3f%n", vectorPga);
        System.out.printf("SPEEDUP,E3,%.2f%n", scalarE3 / vectorE3);
        System.out.printf("SPEEDUP,PGA,%.2f%n", scalarPga / vectorPga);
    }

    private static double time(Runnable task) {
        double best = Double.MAX_VALUE;
        for (int round = 0; round < 7; round++) {
            long t0 = System.nanoTime();
            task.run();
            best = Math.min(best, (System.nanoTime() - t0) / (double) N);
        }
        return best;
    }

    private static double maxDiff(double[]... arrays) {
        double max = 0;
        for (int a = 0; a < arrays.length; a += 2) {
            double[] left = arrays[a];
            double[] right = arrays[a + 1];
            for (int i = 0; i < left.length; i++) {
                max = Math.max(max, Math.abs(left[i] - right[i]));
            }
        }
        return max;
    }

    private static void scalarE3(double s, double b1, double b2, double b3,
                                 double[] x, double[] y, double[] z,
                                 double[] rx, double[] ry, double[] rz) {
        for (int i = 0; i < x.length; i++) {
            double vx = x[i];
            double vy = y[i];
            double vz = z[i];
            rx[i] = (s * s - b1 * b1 - b2 * b2 + b3 * b3) * vx
                    + 2 * (b1 * (s * vy + b3 * vz) + b2 * (s * vz - b3 * vy));
            ry[i] = (s * s - b1 * b1 + b2 * b2 - b3 * b3) * vy
                    + 2 * (-vx * (b1 * s + b2 * b3) + vz * (b3 * s - b1 * b2));
            rz[i] = (s * s + b1 * b1 - b2 * b2 - b3 * b3) * vz
                    + 2 * (vx * (b1 * b3 - b2 * s) - vy * (b1 * b2 + b3 * s));
        }
    }

    private static void scalarPga(double s, double b01, double b02, double b12,
                                  double b03, double b13, double b23, double ps,
                                  double[] p7, double[] p11, double[] p13, double[] p14,
                                  double[] r7, double[] r11, double[] r13, double[] r14) {
        for (int i = 0; i < p7.length; i++) {
            double a7 = p7[i];
            double a11 = p11[i];
            double a13 = p13[i];
            double a14 = p14[i];
            r7[i] = (s * s + b12 * b12 - b13 * b13 - b23 * b23) * a7
                    + 2 * (b12 * (b13 * a11 + b23 * a13) + s * (b23 * a11 - b13 * a13)
                    + a14 * (ps * b12 - b01 * b13 - b02 * b23 + b03 * s));
            r11[i] = (s * s - b12 * b12 + b13 * b13 - b23 * b23) * a11
                    + 2 * (b12 * (b13 * a7 + s * a13) + b13 * b23 * a13 - s * b23 * a7
                    + a14 * (ps * b13 + b01 * b12 - b02 * s - b03 * b23));
            r13[i] = (s * s - b12 * b12 - b13 * b13 + b23 * b23) * a13
                    + 2 * (b12 * (b23 * a7 - s * a11) + b13 * (b23 * a11 + s * a7)
                    + a14 * (ps * b23 + b01 * s + b02 * b12 + b03 * b13));
            r14[i] = (s * s + b12 * b12 + b13 * b13 + b23 * b23) * a14;
        }
    }
}
