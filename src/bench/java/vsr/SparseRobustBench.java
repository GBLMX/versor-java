package vsr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

/**
 * Sparse robustness: structured inputs, mixed densities, strategy crossover,
 * on-demand sign cost, and first-call/table behaviour.
 */
public final class SparseRobustBench {
    private static final int[] NS = {10, 11, 12, 13, 14};

    public static void main(String[] args) throws Exception {
        System.out.println("## structured sparse cases");
        Random rnd = new Random(42);
        for (int n : NS) {
            Signature sig = Signature.euclidean(n);
            int[] ks = n == 10 ? new int[]{8, 32, 128, 512} : new int[]{8, 32, 128, 512};
            for (int k : ks) {
                if (k > sig.size / 2) continue;
                for (String kind : new String[]{"uniform", "versor", "pgaPoint", "cgaPoint"}) {
                    SparseMultivector a = structured(sig, kind, k, rnd);
                    SparseMultivector b = structured(sig, kind, k, rnd);
                    double gp = GaBench.ns(() -> a.gp(b).nnz());
                    double p99 = p99(() -> a.gp(b).nnz(), gp);
                    double err = -1;
                    if (n <= 10) {
                        err = a.gp(b).toDense()
                                .distanceTo(a.toDense().gp(b.toDense()))
                                / Math.max(1, a.toDense().maxAbs() * b.toDense().maxAbs());
                    }
                    System.out.printf("ROBUST,%d,%s,%d,gp,%.1f,%.1f,%d,%.2e%n",
                            n, kind, k, gp, p99, a.gp(b).nnz(), err);
                    if (k <= 32) {
                        double sand = GaBench.ns(() -> a.sandwich(b).nnz());
                        System.out.printf("ROBUST,%d,%s,%d,sandwich,%.1f,%.1f,%d,%.2e%n",
                                n, kind, k, sand, p99(() -> a.sandwich(b).nnz(), sand),
                                a.sandwich(b).nnz(), err);
                    }
                }
            }
        }

        System.out.println("## mixed density (1% dense + 99% sparse), 1000 ops");
        for (int n : new int[]{10, 12, 14}) {
            Signature sig = Signature.euclidean(n);
            List<SparseMultivector[]> ops = new ArrayList<>();
            int denseK = Math.min(sig.size, 1024);
            for (int i = 0; i < 1000; i++) {
                int k = i < 10 ? denseK : 10;
                ops.add(new SparseMultivector[]{
                        uniform(sig, Math.min(k, sig.size), rnd),
                        uniform(sig, Math.min(k, sig.size), rnd)});
            }
            long t0 = System.nanoTime();
            double sink = 0;
            for (SparseMultivector[] p : ops) sink += p[0].gp(p[1]).nnz();
            long dt = System.nanoTime() - t0;
            sink += sink == 0 ? 1 : 0;
            System.out.printf("MIXED,%d,%dpct,%.1f,%.1f%n", n,
                    (int) Math.round(100.0 * denseK / sig.size), dt / 1000.0, sink);
            double sparseOnly = GaBench.ns(() -> uniform(sig, 10, new Random(1))
                    .gp(uniform(sig, 10, new Random(2))).nnz());
            System.out.printf("MIXED,%d,sparse10,%.1f,0%n", n, sparseOnly);
        }

        System.out.println("## auto threshold: forced hash vs dense vs auto");
        for (int n : new int[]{10, 12, 14}) {
            Signature sig = Signature.euclidean(n);
            for (int k : new int[]{2, 8, 16, 32, 64, 128, 256, 512}) {
                if (k > sig.size / 4) continue;
                SparseMultivector a = uniform(sig, k, new Random(n * 31L + k));
                SparseMultivector b = uniform(sig, k, new Random(n * 17L + k));
                double hash = GaBench.ns(() -> a.gpHashed(b).nnz());
                double dense = GaBench.ns(() -> a.gpDense(b).nnz());
                double auto = GaBench.ns(() -> a.gp(b).nnz());
                boolean hashChosen = (long) k * k * (n <= 10 ? 4 : 16) < sig.size;
                boolean correct = hashChosen == (hash <= dense);
                System.out.printf("THRESH,%d,%d,%.1f,%.1f,%.1f,%s,%s%n",
                        n, k, hash, dense, auto, hashChosen ? "hash" : "dense",
                        correct ? "ok" : "miss");
            }
        }

        System.out.println("## dynamic density growth");
        for (int n : new int[]{10, 12}) {
            Signature sig = Signature.euclidean(n);
            SparseMultivector x = structured(sig, "versor", 8, rnd);
            SparseMultivector g = structured(sig, "versor", 8, rnd);
            for (int step = 0; step < 12; step++) {
                SparseMultivector cur = x;
                double ns = GaBench.ns(() -> cur.gp(g).nnz());
                x = cur.gp(g);
                System.out.printf("DYNAMIC,%d,%d,%d,%.1f%n", n, step, x.nnz(), ns);
            }
        }

        System.out.println("## sign table: cold first product, rebuilds, concurrency");
        for (int n = 8; n <= 14; n++) {
            final int dim = n;
            Signature sig = Signature.euclidean(dim);
            long t0 = System.nanoTime();
            coldProduct(sig);
            long first = System.nanoTime() - t0;
            long bytes = dim <= 10 ? 4L * sig.size * sig.size : 0;
            int rounds = 20;
            t0 = System.nanoTime();
            for (int i = 0; i < rounds; i++) coldProduct(Signature.euclidean(dim));
            long rebuild = (System.nanoTime() - t0) / rounds;
            System.out.printf("TABLE,%d,1,%d,%d,%d%n", dim, first / 1000, rebuild / 1000, bytes);
            for (int threads : new int[]{4, 8}) {
                CountDownLatch start = new CountDownLatch(1);
                Thread[] ts = new Thread[threads];
                long[] times = new long[threads];
                for (int t = 0; t < threads; t++) {
                    final int idx = t;
                    final Signature fresh = Signature.euclidean(dim);
                    ts[t] = new Thread(() -> {
                        try {
                            start.await();
                        } catch (InterruptedException e) {
                            return;
                        }
                        long a = System.nanoTime();
                        coldProduct(fresh);
                        times[idx] = System.nanoTime() - a;
                    });
                    ts[t].start();
                }
                t0 = System.nanoTime();
                start.countDown();
                for (Thread t : ts) t.join();
                long wall = System.nanoTime() - t0;
                Arrays.sort(times);
                System.out.printf("TABLE,%d,%d,%d,%d,%d%n", dim, threads, wall / 1000,
                        times[times.length / 2] / 1000, bytes);
            }
        }

        System.out.println("## 3x3 PGA versor matrix via sparse entries");
        Random mrand = new Random(5);
        Multivector[][] entries = new Multivector[3][3];
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                entries[r][c] = P3.motor(P3.rotor(mrand.nextDouble(), mrand.nextDouble(),
                        mrand.nextDouble(), mrand.nextDouble()),
                        P3.translator(mrand.nextDouble(), mrand.nextDouble(), mrand.nextDouble()));
            }
        }
        MultivectorMatrix denseMatrix = new MultivectorMatrix(P3.SIG, 3, 3);
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) denseMatrix.set(r, c, entries[r][c]);
        }
        double denseNs = GaBench.ns(() -> denseMatrix.inverse().get(0, 0).get(0));
        double convertNs = GaBench.ns(() -> {
            MultivectorMatrix m = new MultivectorMatrix(P3.SIG, 3, 3);
            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 3; c++) {
                    m.set(r, c, SparseMultivector.fromDense(entries[r][c]).toDense());
                }
            }
            return m.inverse().get(0, 0).get(0);
        });
        System.out.printf("MATRIX,3x3PGA,denseInverse,%.1f,0%n", denseNs);
        System.out.printf("MATRIX,3x3PGA,sparseConvert,%.1f,0%n", convertNs);
    }

    private static double coldProduct(Signature sig) {
        SparseMultivector a = SparseMultivector.of(sig, new int[]{1}, new double[]{1});
        SparseMultivector b = SparseMultivector.of(sig, new int[]{2}, new double[]{1});
        return a.gp(b).nnz();
    }

    private static SparseMultivector structured(Signature sig, String kind, int k, Random rnd) {
        List<Integer> pool = new ArrayList<>();
        for (int mask = 1; mask < sig.size; mask++) {
            int g = Integer.bitCount(mask);
            boolean take = switch (kind) {
                case "versor" -> (g == 2 || g == 4);
                case "pgaPoint" -> g == 3;
                case "cgaPoint" -> g == 1;
                default -> true;
            };
            if (take) pool.add(mask);
        }
        if (pool.isEmpty()) pool.add(1);
        Collections.shuffle(pool, rnd);
        int count = Math.min(k, pool.size());
        int[] masks = new int[count];
        double[] vals = new double[count];
        for (int i = 0; i < count; i++) {
            masks[i] = pool.get(i);
            vals[i] = rnd.nextDouble() * 2 - 1;
        }
        return SparseMultivector.of(sig, masks, vals);
    }

    private static SparseMultivector uniform(Signature sig, int k, Random rnd) {
        int count = Math.min(k, sig.size);
        int[] masks = new int[count];
        double[] vals = new double[count];
        for (int i = 0; i < count; i++) {
            masks[i] = rnd.nextInt(sig.size);
            vals[i] = rnd.nextDouble() * 2 - 1;
        }
        return SparseMultivector.of(sig, masks, vals);
    }

    private static double p99(GaBench.Op op, double measuredNs) {
        int samples = measuredNs > 1_000_000 ? 5 : 50;
        long[] times = new long[samples];
        for (int i = 0; i < 10; i++) op.apply();
        for (int i = 0; i < samples; i++) {
            long t0 = System.nanoTime();
            op.apply();
            times[i] = System.nanoTime() - t0;
        }
        Arrays.sort(times);
        return times[(int) Math.min(samples - 1, samples * 0.99)];
    }
}
