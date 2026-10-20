package vsr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/** Multicore scaling for dense GP, O(3^n) wedge, scalar and batched sandwich. */
public final class MultiCoreBench {
    private static final int[] THREADS = {1, 2, 4, 8, 14, 20};
    private static final int REPS = 5;

    public static void main(String[] args) throws Exception {
        System.out.println("cpus=" + Runtime.getRuntime().availableProcessors());
        Signature sig = Signature.euclidean(10);
        Random rnd = new Random(2026);
        Multivector a = dense(sig, rnd);
        Multivector b = dense(sig, rnd);
        a.gp(b); // build the shared sign table before timing
        a.op(b);
        for (int threads : THREADS) {
            scale("gp_n10", threads, 512, () -> a.gp(b).get(0));
        }
        for (int threads : THREADS) {
            scale("op_n10", threads, 256, () -> a.op(b).get(0));
        }

        int points = 2_000_000;
        MultivectorBatch all = points(points);
        Multivector motor = P3.motor(P3.rotor(0, 0, 1, 0.4), P3.translator(1, 2, 3));
        for (int threads : THREADS) {
            batchScale(threads, all, motor);
        }
        for (int threads : THREADS) {
            scalarScale(threads, points, motor);
        }
    }

    private static void scale(String name, int threads, int opsPerThread,
                              Callable<Double> work) throws Exception {
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        warmup(pool, threads, opsPerThread, work);
        double[] times = new double[REPS];
        for (int r = 0; r < REPS; r++) {
            CountDownLatch start = new CountDownLatch(1);
            List<Future<Double>> futures = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                futures.add(pool.submit(() -> {
                    start.await();
                    double sink = 0;
                    for (int i = 0; i < opsPerThread; i++) sink += work.call();
                    return sink;
                }));
            }
            long t0 = System.nanoTime();
            start.countDown();
            for (Future<Double> f : futures) f.get();
            times[r] = System.nanoTime() - t0;
        }
        pool.shutdown();
        Arrays.sort(times);
        long totalOps = (long) threads * opsPerThread;
        double nsPerOp = times[REPS / 2] / totalOps;
        double base = reference(name);
        System.out.printf("MC,%s,%d,%d,%.1f,%.2f,%.2f%n",
                name, threads, totalOps, nsPerOp,
                base > 0 ? base / nsPerOp : 0, base > 0 ? base / nsPerOp / threads : 0);
    }

    private static void warmup(ExecutorService pool, int threads, int ops,
                               Callable<Double> work) throws Exception {
        List<Future<Double>> fs = new ArrayList<>();
        for (int t = 0; t < threads; t++) {
            fs.add(pool.submit(() -> {
                double sink = 0;
                for (int i = 0; i < Math.max(1, ops / 8); i++) sink += work.call();
                return sink;
            }));
        }
        for (Future<Double> f : fs) f.get();
    }

    private static double reference(String name) {
        return switch (name) {
            case "gp_n10" -> 850_000;
            case "op_n10" -> 77_000;
            case "batch_sandwich" -> 19.4;
            case "scalar_sandwich" -> 31.0;
            default -> 0;
        };
    }

    private static void scalarScale(int threads, int points, Multivector motor)
            throws Exception {
        int chunk = (points + threads - 1) / threads;
        List<Multivector> shared = new ArrayList<>();
        for (int i = 0; i < chunk; i++) shared.add(P3.point(i, -i, 0.5));
        shared = List.copyOf(shared);
        final List<Multivector> pts = shared;
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        for (int w = 0; w < 2; w++) {
            List<Future<Double>> fs = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                fs.add(pool.submit(() -> {
                    double sum = 0;
                    for (Multivector p : pts) sum += motor.sandwich(p).get(14);
                    return sum;
                }));
            }
            for (Future<Double> f : fs) f.get();
        }
        double[] times = new double[REPS];
        for (int r = 0; r < REPS; r++) {
            CountDownLatch start = new CountDownLatch(1);
            List<Future<Double>> fs = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                fs.add(pool.submit(() -> {
                    start.await();
                    double sum = 0;
                    for (Multivector p : pts) sum += motor.sandwich(p).get(14);
                    return sum;
                }));
            }
            long t0 = System.nanoTime();
            start.countDown();
            for (Future<Double> f : fs) f.get();
            times[r] = System.nanoTime() - t0;
        }
        pool.shutdown();
        Arrays.sort(times);
        double nsPerPoint = times[REPS / 2] / (double) points;
        double base = reference("scalar_sandwich");
        System.out.printf("MC,scalar_sandwich,%d,%d,%.1f,%.2f,%.2f%n",
                threads, points, nsPerPoint, base / nsPerPoint,
                base / nsPerPoint / threads);
    }

    private static void batchScale(int threads, MultivectorBatch all, Multivector motor)
            throws Exception {
        int n = all.size();
        int chunk = (n + threads - 1) / threads;
        List<MultivectorBatch> slices = new ArrayList<>();
        List<MultivectorBatch> outs = new ArrayList<>();
        for (int t = 0; t < threads; t++) {
            int from = t * chunk;
            int to = Math.min(n, from + chunk);
            slices.add(all.slice(from, to));
            outs.add(new MultivectorBatch(P3.SIG, to - from));
        }
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        double sink = 0;
        for (int r = 0; r < 2; r++) {
            CountDownLatch start = new CountDownLatch(1);
            List<Future<Double>> fs = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                final int idx = t;
                fs.add(pool.submit(() -> {
                    start.await();
                    slices.get(idx).sandwich(motor, outs.get(idx));
                    return outs.get(idx).get(14, 0);
                }));
            }
            start.countDown();
            for (Future<Double> f : fs) sink += f.get();
        }
        double[] times = new double[REPS];
        for (int r = 0; r < REPS; r++) {
            CountDownLatch start = new CountDownLatch(1);
            List<Future<Double>> fs = new ArrayList<>();
            for (int t = 0; t < threads; t++) {
                final int idx = t;
                fs.add(pool.submit(() -> {
                    start.await();
                    slices.get(idx).sandwich(motor, outs.get(idx));
                    return outs.get(idx).get(14, 0);
                }));
            }
            long t0 = System.nanoTime();
            start.countDown();
            for (Future<Double> f : fs) sink += f.get();
            times[r] = System.nanoTime() - t0;
        }
        pool.shutdown();
        Arrays.sort(times);
        double nsPerPoint = times[REPS / 2] / (double) n;
        double base = reference("batch_sandwich");
        System.out.printf("MC,batch_sandwich,%d,%d,%.1f,%.2f,%.2f,%.1f%n",
                threads, n, nsPerPoint, base / nsPerPoint, base / nsPerPoint / threads,
                2.0 * 16 * 8 / nsPerPoint);
        if (sink == 42) System.out.print("");
    }

    private static Multivector dense(Signature sig, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < c.length; i++) c[i] = rnd.nextDouble() * 2 - 1;
        return new Multivector(sig, c);
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
}
