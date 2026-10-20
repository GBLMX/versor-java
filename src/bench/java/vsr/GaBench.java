package vsr;

import java.lang.management.ManagementFactory;
import java.util.Random;

/**
 * Standalone performance harness for the geometric-algebra core (no game or
 * Gradle dependencies). Every measurement calibrates its own iteration count,
 * warms up the JIT, then reports the best of five rounds.
 *
 * <pre>
 * java -cp out-bench vsr.GaBench            # full sweep
 * java -cp out-bench vsr.GaBench jit        # key operations only (for JIT flags)
 * java -cp out-bench vsr.GaBench batch      # batch scaling only
 * java -cp out-bench vsr.GaBench stress     # matrix / high-dimensional stress
 * </pre>
 */
public final class GaBench {
    private static final int WARMUP = 3;
    private static final int ROUNDS = 5;
    private static final long TARGET_NANOS = 40_000_000L;
    private static final long MAX_ITERS = 1L << 24;
    private static volatile double sink;

    interface Op {
        double apply();
    }

    public static void main(String[] args) {
        String mode = args.length == 0 ? "all" : args[0];
        if (mode.equals("all") || mode.equals("sweep")) sweep();
        if (mode.equals("all") || mode.equals("batch")) batch();
        if (mode.equals("all") || mode.equals("stress")) stress();
        if (mode.equals("all") || mode.equals("jit")) jit();
    }

    // ---------------------------------------------------------------- sweep

    private static void sweep() {
        System.out.println("## dimension sweep (Euclidean Cl(n,0), dense inputs)");
        int[] dims = {3, 4, 5, 6, 7, 8, 9, 10};
        for (int n : dims) {
            Signature sig = Signature.euclidean(n);
            Random rnd = new Random(1000 + n);
            Multivector a = dense(sig, rnd);
            Multivector b = dense(sig, rnd);
            Multivector v = Multivector.scalar(sig, 1.7)
                    .add(Multivector.blade(sig, 1, 0.3))
                    .add(Multivector.blade(sig, 2, -0.4))
                    .add(Multivector.blade(sig, 4, 0.8));
            double size = sig.size;
            System.out.printf("n=%d size=%.0f bytesPerMv=%.0f%n", n, size, 16 + 16 + 8 * size);
            line(n, "add", ns(() -> a.add(b).get(0)));
            line(n, "gp", ns(() -> a.gp(b).get(0)));
            line(n, "op", ns(() -> a.op(b).get(0)));
            line(n, "ip", ns(() -> a.ip(b).get(0)));
            line(n, "rightContraction", ns(() -> a.rightContraction(b).get(0)));
            line(n, "vee", ns(() -> a.vee(b).get(0)));
            line(n, "commutator", ns(() -> a.commutator(b).get(0)));
            line(n, "reverse", ns(() -> a.reverse().get(0)));
            line(n, "conjugate", ns(() -> a.conjugate().get(0)));
            line(n, "involution", ns(() -> a.involution().get(0)));
            line(n, "dual", ns(() -> a.dual().get(0)));
            line(n, "grade1", ns(() -> a.grade(1).get(0)));
            line(n, "norm2", ns(() -> a.norm2()));
            line(n, "unit", ns(() -> a.unit().get(0)));
            line(n, "inverse", ns(() -> v.inverse().get(0)));
            line(n, "sandwich", ns(() -> a.sandwich(b).get(0)));
            if (n <= 6) {
                line(n, "inverseGeneral", ns(() -> v.inverseGeneral().get(0)));
                Multivector small = v.scale(0.1);
                line(n, "exp", ns(() -> small.exp().get(0)));
            }
            System.out.println();
        }
    }

    private static Multivector dense(Signature sig, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < c.length; i++) c[i] = rnd.nextDouble() * 2 - 1;
        return new Multivector(sig, c);
    }

    private static void line(int n, String op, double nanos) {
        System.out.printf("DIM,%d,%s,%.2f%n", n, op, nanos);
    }

    // ---------------------------------------------------------------- batch

    private static void batch() {
        System.out.println("## batch scaling: one P3 motor over N points");
        Random rnd = new Random(5);
        int max = 100_000;
        Multivector[] points = new Multivector[max];
        for (int i = 0; i < max; i++) {
            points[i] = P3.point(rnd.nextDouble() * 4 - 2,
                    rnd.nextDouble() * 4 - 2, rnd.nextDouble() * 4 - 2);
        }
        Multivector motor = P3.motor(P3.rotor(0.3, -0.7, 0.2, 0.9),
                P3.translator(1.5, -2.5, 0.75));
        for (int n : new int[]{1, 10, 100, 1_000, 10_000, 100_000}) {
            Multivector[] slice = new Multivector[n];
            System.arraycopy(points, 0, slice, 0, n);
            MultivectorBatch batch = MultivectorBatch.of(slice);
            final int count = n;
            double loop = ns(() -> {
                double s = 0;
                for (int i = 0; i < count; i++) s += motor.sandwich(points[i]).get(0);
                return s;
            });
            double batched = ns(() -> batch.sandwich(motor).get(0, count - 1));
            Multivector[] samplePoints = new Multivector[1000];
            System.arraycopy(points, 0, samplePoints, 0, 1000);
            MultivectorBatch sample = MultivectorBatch.of(samplePoints);
            System.out.printf("BATCH,%d,loop=%.2f,batch=%.2f,loopB=%.0f,batchB=%.0f%n",
                    n, loop, batched, bytesPerOp(() -> motor.sandwich(points[0]).get(0), 1000),
                    bytesPerOp(() -> sample.sandwich(motor).get(0, 0), 10) / 1000);
        }

        System.out.println("## dense vs sparse sandwich, Cl(4)");
        Signature sig = Signature.euclidean(4);
        Random r2 = new Random(9);
        Multivector denseA = dense(sig, r2);
        Multivector denseB = dense(sig, r2);
        Multivector motorLike = Multivector.scalar(sig, 0.9)
                .add(Multivector.blade(sig, 3, 0.2)).add(Multivector.blade(sig, 6, -0.3));
        Multivector pointLike = Multivector.blade(sig, 14, 1).add(Multivector.blade(sig, 7, -0.4));
        System.out.printf("SPARSE,denseSandwich=%.2f,motorSandwich=%.2f,motorGp=%.2f%n",
                ns(() -> denseA.sandwich(denseB).get(0)),
                ns(() -> motorLike.sandwich(pointLike).get(0)),
                ns(() -> motorLike.gp(pointLike).get(0)));
        Multivector e3r = E3.rotor(1, 2, 3, 0.7);
        Multivector e3v = E3.vec(1, 0, 0);
        Multivector p3m = P3.motor(P3.rotor(0, 0, 1, 0.6), P3.translator(1, 2, 3));
        Multivector p3p = P3.point(1, 2, 3);
        System.out.printf("FASTPATH,e3Sandwich=%.2f,e3Generic=%.2f,p3Sandwich=%.2f,p3Generic=%.2f%n",
                ns(() -> e3r.sandwich(e3v).get(0)),
                ns(() -> e3r.gp(e3v).gp(e3r.reverse()).get(0)),
                ns(() -> p3m.sandwich(p3p).get(0)),
                ns(() -> p3m.gp(p3p).gp(p3m.reverse()).get(0)));
    }

    // ---------------------------------------------------------------- stress

    private static void stress() {
        System.out.println("## stress: high-dimensional dense sandwich / inverse");
        for (int n = 4; n <= 8; n++) {
            Signature sig = Signature.euclidean(n);
            Random rnd = new Random(2000 + n);
            Multivector a = dense(sig, rnd);
            Multivector b = dense(sig, rnd);
            Multivector v = dense(sig, rnd).scale(0.1).add(Multivector.scalar(sig, 1.5));
            line(n, "denseSandwich", ns(() -> a.sandwich(b).get(0)));
            line(n, "inverseGeneral", ns(() -> v.inverseGeneral().get(0)));
            if (n <= 6) {
                MultivectorMatrix m = new MultivectorMatrix(P3.SIG, 3, 3);
                Random rm = new Random(3000 + n);
                for (int r = 0; r < 3; r++) {
                    for (int c = 0; c < 3; c++) {
                        m.set(r, c, P3.motor(
                                P3.rotor(rm.nextDouble() - 0.5, rm.nextDouble() - 0.5,
                                        rm.nextDouble() - 0.5, rm.nextDouble() * 2),
                                P3.translator(rm.nextDouble(), rm.nextDouble(), rm.nextDouble())));
                    }
                }
                line(n, "matrix3inverse", ns(() -> m.inverse().get(0, 0).get(0)));
            }
        }
    }

    // ---------------------------------------------------------------- jit

    private static void jit() {
        System.out.println("## jit: key operations (run with different JVM flags)");
        Signature sig = Signature.euclidean(5);
        Random rnd = new Random(11);
        Multivector a = dense(sig, rnd);
        Multivector b = dense(sig, rnd);
        line(5, "gp", ns(() -> a.gp(b).get(0)));
        line(5, "inverse", ns(() -> a.inverse().get(0)));
        line(5, "sandwichDense", ns(() -> a.sandwich(b).get(0)));
        Multivector motor = P3.motor(P3.rotor(0, 0, 1, 0.6), P3.translator(1, 2, 3));
        Multivector point = P3.point(1, 2, 3);
        line(5, "sandwichP3", ns(() -> motor.sandwich(point).get(0)));
    }

    // ---------------------------------------------------------------- timing

    static double ns(Op op) {
        long iters = 1;
        double dummy = 0;
        while (true) {
            long t0 = System.nanoTime();
            for (long i = 0; i < iters; i++) dummy += op.apply();
            long dt = System.nanoTime() - t0;
            if (dt >= TARGET_NANOS || iters >= MAX_ITERS) break;
            iters <<= 1;
        }
        sink += dummy;
        for (int w = 0; w < WARMUP; w++) run(op, iters);
        double best = Double.MAX_VALUE;
        for (int r = 0; r < ROUNDS; r++) {
            double value = run(op, iters) / (double) iters;
            if (value < best) best = value;
        }
        return best;
    }

    private static double run(Op op, long iters) {
        double sum = 0;
        long t0 = System.nanoTime();
        for (long i = 0; i < iters; i++) sum += op.apply();
        long dt = System.nanoTime() - t0;
        sink += sum;
        return dt;
    }

    static double bytesPerOp(Op op, int iterations) {
        long iters = iterations;
        for (int w = 0; w < 2; w++) run(op, iters);
        long a0 = allocatedBytes();
        if (a0 < 0) return -1;
        run(op, iters);
        long a1 = allocatedBytes();
        return (a1 - a0) / (double) iters;
    }

    private static long allocatedBytes() {
        try {
            com.sun.management.ThreadMXBean bean =
                    (com.sun.management.ThreadMXBean) ManagementFactory.getThreadMXBean();
            return bean.getCurrentThreadAllocatedBytes();
        } catch (Throwable t) {
            return -1;
        }
    }
}
