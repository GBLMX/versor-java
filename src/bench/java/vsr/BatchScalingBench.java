package vsr;

import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Batch memory behaviour: allocation vs output reuse, block reuse, thread
 * scaling, and the cost of misaligned block boundaries.
 */
public final class BatchScalingBench {
    public static void main(String[] args) throws Exception {
        System.out.println("## single thread: allocation vs reuse (PGA points)");
        for (int n : new int[]{100_000, 1_000_000, 10_000_000}) {
            final MultivectorBatch batch = points(n);
            final Multivector motor = motor(0.4);
            final MultivectorBatch out = new MultivectorBatch(P3.SIG, n);
            final MultivectorBatch out2 = new MultivectorBatch(P3.SIG, n);
            int reps = Math.max(3, 2_000_000 / n);
            variant("alloc", n, reps, () -> batch.sandwich(motor));
            variant("reuse", n, reps, () -> batch.sandwich(motor, out));
            boolean[] flip = {false};
            variant("double", n, reps, () -> {
                flip[0] = !flip[0];
                return batch.sandwich(motor, flip[0] ? out : out2);
            });
            if (n <= 1_000_000) {
                Multivector[] pts = batch.toMultivectors();
                variant("aos", n, reps, () -> {
                    double sum = 0;
                    for (Multivector p : pts) sum += motor.sandwich(p).get(14);
                    return sum;
                });
            }
            System.gc();
        }

        System.out.println("## thread scaling with per-thread output reuse");
        int n = 1_000_000;
        MultivectorBatch all = points(n);
        for (int threads : new int[]{1, 2, 4, 8, 14}) {
            int chunk = (n + threads - 1) / threads;
            final List<MultivectorBatch> slices = new ArrayList<>();
            final List<MultivectorBatch> outs = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                int from = t * chunk;
                int to = Math.min(n, from + chunk);
                slices.add(all.slice(from, to));
                outs.add(new MultivectorBatch(P3.SIG, to - from));
            }
            Multivector motor = motor(0.4);
            ExecutorService pool = Executors.newFixedThreadPool(threads);
            int reps = 10;
            long t0 = System.nanoTime();
            for (int r = 0; r < reps; r++) {
                List<Future<Object>> fs = new ArrayList<>();
                for (int t = 0; t < threads; t++) {
                    final int idx = t;
                    fs.add(pool.submit((Callable<Object>) () -> {
                        slices.get(idx).sandwich(motor, outs.get(idx));
                        return null;
                    }));
                }
                for (Future<Object> f : fs) f.get();
            }
            long dt = System.nanoTime() - t0;
            pool.shutdown();
            double nsPerPoint = dt / (double) n / reps;
            System.out.printf("SCALE,%d,%d,reuse,%.1f,%.2f%n", threads, 0,
                    1e9 / nsPerPoint, nsPerPoint);
            System.gc();
        }

        System.out.println("## block size and alignment (8 threads, 1M points)");
        Multivector motor = motor(0.4);
        for (int block : new int[]{1_000, 10_000, 100_000, 10_001}) {
            ExecutorService pool = Executors.newFixedThreadPool(8);
            int blocks = (n + block - 1) / block;
            long t0 = System.nanoTime();
            List<Future<Object>> fs = new ArrayList<>();
            for (int b = 0; b < blocks; b++) {
                final int from = b * block;
                final int len = Math.min(block, n - from);
                fs.add(pool.submit((Callable<Object>) () -> {
                    MultivectorBatch slice = all.slice(from, from + len);
                    MultivectorBatch out = new MultivectorBatch(P3.SIG, len);
                    slice.sandwich(motor, out);
                    return null;
                }));
            }
            for (Future<Object> f : fs) f.get();
            long dt = System.nanoTime() - t0;
            pool.shutdown();
            double ns = dt / (double) n;
            int aligned = (block * 8) % 64 == 0 ? 1 : 0;
            System.out.printf("BLOCK,%d,%d,%d,%.2f%n", block, aligned, blocks, ns);
        }
    }

    private static MultivectorBatch points(int n) {
        MultivectorBatch batch = new MultivectorBatch(P3.SIG, n);
        Multivector p = P3.point(1.5, -2.0, 0.75);
        for (int m = 0; m < P3.SIG.size; m++) {
            double v = p.get(m);
            if (v != 0) {
                for (int i = 0; i < n; i++) batch.set(m, i, v);
            }
        }
        return batch;
    }

    private static Multivector motor(double angle) {
        return P3.motor(P3.rotor(0, 0, 1, angle), P3.translator(1, 2, 3));
    }

    private static void variant(String name, int n, int reps, Callable<Object> body) throws Exception {
        for (int i = 0; i < Math.min(3, reps); i++) body.call();
        long gc0 = gcCount();
        long gct0 = gcTime();
        long a0 = allocated();
        long t0 = System.nanoTime();
        for (int i = 0; i < reps; i++) body.call();
        long dt = System.nanoTime() - t0;
        long a1 = allocated();
        long gc1 = gcCount();
        long gct1 = gcTime();
        double ns = dt / (double) n / reps;
        double bytes = a1 < 0 ? -1 : (a1 - a0) / (double) n / reps;
        double gbs = 256.0 / ns; // 16 doubles in + 16 doubles out per point
        System.out.printf("REUSE,%d,%s,%.2f,%.1f,%.1f,%d,%d%n",
                n, name, ns, bytes, gbs, gc1 - gc0, gct1 - gct0);
    }

    private static long gcCount() {
        long c = 0;
        for (GarbageCollectorMXBean b : ManagementFactory.getGarbageCollectorMXBeans()) c += b.getCollectionCount();
        return c;
    }

    private static long gcTime() {
        long c = 0;
        for (GarbageCollectorMXBean b : ManagementFactory.getGarbageCollectorMXBeans()) c += b.getCollectionTime();
        return c;
    }

    private static long allocated() {
        try {
            return ((com.sun.management.ThreadMXBean) ManagementFactory.getThreadMXBean())
                    .getCurrentThreadAllocatedBytes();
        } catch (Throwable t) {
            return -1;
        }
    }
}
