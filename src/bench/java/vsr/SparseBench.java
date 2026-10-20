package vsr;

import java.util.Random;

/**
 * High-dimensional sparse/COO stress: does sparsity break the dense 4^n cost?
 * Prints CSV lines for the report.
 */
public final class SparseBench {
    public static void main(String[] args) {
        System.out.println("## sparse / high-dimensional stress");
        int[] densities = {2, 8, 10, 32, 128, 512};
        for (int n = 8; n <= 14; n++) {
            Signature sig = Signature.euclidean(n);
            int size = sig.size;
            long t0 = System.nanoTime();
            int[][] table = sig.signTable();
            long tableMs = (System.nanoTime() - t0) / 1_000_000;
            long tableBytes = table == null ? 0 : 4L * size * size;
            System.out.printf("TABLE,%d,%d,%d,%d%n", n, tableMs, tableBytes, table == null ? 0 : 1);

            for (int k : densities) {
                if (k > size) continue;
                Random rnd = new Random(1000L * n + k);
                int[] ma = distinctMasks(rnd, size, k);
                int[] mb = distinctMasks(rnd, size, k);
                double[] va = new double[k];
                double[] vb = new double[k];
                for (int i = 0; i < k; i++) {
                    va[i] = rnd.nextDouble() * 2 - 1;
                    vb[i] = rnd.nextDouble() * 2 - 1;
                }
                SparseMultivector sa = SparseMultivector.of(sig, ma, va);
                SparseMultivector sb = SparseMultivector.of(sig, mb, vb);
                long pairs = (long) sa.nnz() * sb.nnz();
                double gp = GaBench.ns(() -> sa.gp(sb).nnz());
                double sand = GaBench.ns(() -> sa.sandwich(sb).nnz());
                double bytes = GaBench.bytesPerOp(() -> sa.gp(sb).nnz(), 20);
                System.out.printf("SPARSE,%d,%d,gp,%.2f,%d,%.3f,%.0f%n",
                        n, k, gp, pairs, gp / pairs, bytes);
                System.out.printf("SPARSE,%d,%d,sandwich,%.2f,%d,%.3f,%.0f%n",
                        n, k, sand, pairs, sand / pairs, 0.0);
            }

            if (n <= 10) {
                Random rnd = new Random(2000L + n);
                Multivector da = dense(sig, rnd);
                Multivector db = dense(sig, rnd);
                double denseGp = GaBench.ns(() -> da.gp(db).get(0));
                double denseSand = GaBench.ns(() -> da.sandwich(db).get(0));
                System.out.printf("DENSE,%d,gp,%.2f%n", n, denseGp);
                System.out.printf("DENSE,%d,sandwich,%.2f%n", n, denseSand);
            }
        }
    }

    private static int[] distinctMasks(Random rnd, int size, int k) {
        int[] masks = new int[k];
        java.util.HashSet<Integer> seen = new java.util.HashSet<>();
        int p = 0;
        while (p < k) {
            int m = rnd.nextInt(size);
            if (seen.add(m)) masks[p++] = m;
        }
        return masks;
    }

    private static Multivector dense(Signature sig, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < c.length; i++) c[i] = rnd.nextDouble() * 2 - 1;
        return new Multivector(sig, c);
    }
}
