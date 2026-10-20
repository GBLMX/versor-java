package vsr;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/** Fast paths, batch vs AoS layout, and inverseGeneral. */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 400, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 4, time = 400, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Benchmark)
public class JmhFast {
    static final int N = 10_000;

    private Multivector rotor;
    private Multivector vec;
    private Multivector motor;
    private Multivector point;
    private Multivector[] points;
    private MultivectorBatch batch;
    private Multivector general;

    @Setup
    public void setup() {
        rotor = E3.rotor(1, 2, 3, 0.7);
        vec = E3.vec(1.5, -2.0, 0.25);
        motor = P3.motor(P3.rotor(0, 0, 1, 0.6), P3.translator(1, 2, 3));
        point = P3.point(1.5, 2.5, -0.5);
        Random rnd = new Random(7);
        points = new Multivector[N];
        for (int i = 0; i < N; i++) {
            points[i] = P3.point(rnd.nextDouble() * 4 - 2, rnd.nextDouble() * 4 - 2,
                    rnd.nextDouble() * 4 - 2);
        }
        batch = MultivectorBatch.of(points);
        Signature sig4 = Signature.euclidean(4);
        general = Multivector.scalar(sig4, 1.3).add(Multivector.blade(sig4, 1, 0.4));
    }

    @Benchmark
    public double e3Rotate() {
        return rotor.sandwich(vec).get(1);
    }

    @Benchmark
    public double p3Transform() {
        return motor.sandwich(point).get(14);
    }

    @Benchmark
    @OperationsPerInvocation(N)
    public double loopTransform() {
        double sum = 0;
        for (int i = 0; i < N; i++) sum += motor.sandwich(points[i]).get(14);
        return sum;
    }

    @Benchmark
    @OperationsPerInvocation(N)
    public double batchTransform() {
        return batch.sandwich(motor).get(14, N - 1);
    }

    @Benchmark
    public double inverseGeneral() {
        return general.inverseGeneral().get(0);
    }
}
