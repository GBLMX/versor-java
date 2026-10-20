package vsr;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.joml.Matrix4f;
import org.joml.Vector3f;

import vsr.blaze3d.Blaze3D;

/**
 * End-to-end hot paths: per-frame transforms (GA, batch, cached Matrix4f,
 * adapter), thread scaling, and the robot-style 3x3 multivector inverse.
 */
public final class RenderBench {
    public static void main(String[] args) throws Exception {
        int n = args.length > 0 ? Integer.parseInt(args[0]) : 200_000;
        int frames = args.length > 1 ? Integer.parseInt(args[1]) : 300;
        Random rnd = new Random(11);
        Multivector[] points = new Multivector[n];
        Vector3f[] vectors = new Vector3f[n];
        for (int i = 0; i < n; i++) {
            double x = rnd.nextDouble() * 100 - 50;
            double y = rnd.nextDouble() * 100 - 50;
            double z = rnd.nextDouble() * 100 - 50;
            points[i] = P3.point(x, y, z);
            vectors[i] = new Vector3f((float) x, (float) y, (float) z);
        }
        MultivectorBatch batch = MultivectorBatch.of(points);
        Multivector base = P3.motor(P3.rotor(0, 0, 1, 0.4), P3.translator(1, 2, 3));
        Vector3f scratch = new Vector3f();

        System.out.println("## end-to-end render paths");
        frame("ga", frames, n, f -> {
            Multivector motor = spin(base, f, frames);
            double sum = 0;
            for (int i = 0; i < n; i++) sum += motor.sandwich(points[i]).get(14);
            return sum;
        });
        frame("batch", frames, n, f -> batch.sandwich(spin(base, f, frames)).get(14, n - 1));
        frame("matrix", frames, n, f -> {
            Matrix4f m = Blaze3D.matrix(spin(base, f, frames));
            double sum = 0;
            for (int i = 0; i < n; i++) {
                Vector3f p = vectors[i];
                sum += m.m00() * p.x + m.m10() * p.y + m.m20() * p.z + m.m30();
            }
            return sum;
        });
        int adapterN = Math.min(n, 100_000);
        frame("adapter", frames, adapterN, f -> {
            Multivector motor = spin(base, f, frames);
            double sum = 0;
            for (int i = 0; i < adapterN; i++) {
                sum += Blaze3D.transform(motor, vectors[i], scratch).x;
            }
            return sum;
        });

        System.out.println("## thread scaling (batch sandwich, N=200000)");
        for (int threads : new int[]{1, 2, 4, 8, 14}) {
            Multivector[][] slices = new Multivector[threads][];
            int chunk = (n + threads - 1) / threads;
            for (int t = 0; t < threads; t++) {
                int from = t * chunk;
                int to = Math.min(n, from + chunk);
                slices[t] = new Multivector[Math.max(0, to - from)];
                System.arraycopy(points, from, slices[t], 0, slices[t].length);
            }
            MultivectorBatch[] batches = new MultivectorBatch[threads];
            for (int t = 0; t < threads; t++) batches[t] = MultivectorBatch.of(slices[t]);
            ExecutorService pool = Executors.newFixedThreadPool(threads);
            Multivector motor = base;
            int repeats = 30;
            long t0 = System.nanoTime();
            for (int r = 0; r < repeats; r++) {
                java.util.List<Future<Double>> fs = new java.util.ArrayList<>(threads);
                for (int t = 0; t < threads; t++) {
                    final int idx = t;
                    fs.add(pool.submit((Callable<Double>) () ->
                            batches[idx].sandwich(motor).get(14, batches[idx].size() - 1)));
                }
                for (Future<Double> f : fs) f.get();
            }
            long dt = System.nanoTime() - t0;
            pool.shutdown();
            double pointsPerSec = (double) n * repeats / (dt / 1e9);
            System.out.printf("SCALE,%d,%.0f,%.2f%n", threads, pointsPerSec,
                    dt / (double) n / repeats);
        }

        System.out.println("## robot-style solver");
        MultivectorMatrix m = new MultivectorMatrix(P3.SIG, 3, 3);
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                m.set(r, c, P3.motor(P3.rotor(rnd.nextDouble(), rnd.nextDouble(),
                        rnd.nextDouble(), rnd.nextDouble()), P3.translator(
                        rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble())));
            }
        }
        int iters = 200;
        long t0 = System.nanoTime();
        for (int i = 0; i < iters; i++) m.inverse().get(0, 0).get(0);
        long t1 = System.nanoTime();
        System.out.printf("ROBOT,matrix3Inverse,%.2f%n", (t1 - t0) / (double) iters);

        Multivector[] chain = {
                P3.motor(P3.rotor(0, 0, 1, 0.1), P3.translator(0.1, 0, 0)),
                P3.motor(P3.rotor(0, 1, 0, 0.2), P3.translator(0, 0.1, 0)),
                P3.motor(P3.rotor(1, 0, 0, 0.3), P3.translator(0, 0, 0.1)),
                P3.motor(P3.rotor(0, 0, 1, -0.2), P3.translator(0.05, 0, 0))
        };
        int solverIters = 20_000;
        for (int i = 0; i < 5_000; i++) {
            Multivector p = P3.point(0, 0, 0);
            for (Multivector motor : chain) p = motor.sandwich(p);
        }
        t0 = System.nanoTime();
        for (int i = 0; i < solverIters; i++) {
            Multivector p = P3.point(0, 0, 0);
            for (Multivector motor : chain) p = motor.sandwich(p);
        }
        t1 = System.nanoTime();
        System.out.printf("ROBOT,chain4Forward,%.2f%n", (t1 - t0) / (double) solverIters);
    }

    private interface FrameBody {
        double run(int frame);
    }

    private static void frame(String name, int frames, int n, FrameBody body) {
        for (int i = 0; i < 20; i++) body.run(i);
        double[] times = new double[frames];
        for (int f = 0; f < frames; f++) {
            long t0 = System.nanoTime();
            body.run(f);
            times[f] = System.nanoTime() - t0;
        }
        Arrays.sort(times);
        double sum = 0;
        for (double t : times) sum += t;
        long a0 = allocated();
        body.run(0);
        long a1 = allocated();
        System.out.printf("RENDER,%s,%d,%d,%.0f,%.0f,%.2f,%.0f%n",
                name, n, frames, sum / frames, times[(int) (frames * 0.99)],
                sum / frames / n, a1 < 0 || a0 < 0 ? -1 : (double) (a1 - a0));
    }

    private static Multivector spin(Multivector base, int frame, int frames) {
        double angle = 0.4 + 2 * Math.PI * frame / frames;
        return P3.motor(P3.rotor(0, 0, 1, angle), P3.translator(1, 2, 3));
    }

    private static long allocated() {
        try {
            com.sun.management.ThreadMXBean bean =
                    (com.sun.management.ThreadMXBean) java.lang.management.ManagementFactory.getThreadMXBean();
            return bean.getCurrentThreadAllocatedBytes();
        } catch (Throwable t) {
            return -1;
        }
    }
}
