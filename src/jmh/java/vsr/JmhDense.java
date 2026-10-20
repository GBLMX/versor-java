package vsr;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/** Dense operation matrix over dimensions; run with -prof gc for allocation. */
@BenchmarkMode({Mode.AverageTime, Mode.SampleTime})
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 4, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
public class JmhDense {
    @Param({"3", "4", "5", "6", "7", "8", "9", "10"})
    public int n;

    private Multivector a;
    private Multivector b;
    private Multivector versor;

    @Setup
    public void setup() {
        Signature sig = Signature.euclidean(n);
        Random rnd = new Random(100L + n);
        a = random(sig, rnd);
        b = random(sig, rnd);
        double[] c = new double[sig.size];
        c[0] = 1.2;
        c[3] = 0.4;
        if (sig.size > 5) c[5] = -0.3;
        versor = new Multivector(sig, c);
    }

    private static Multivector random(Signature sig, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < c.length; i++) c[i] = rnd.nextDouble() * 2 - 1;
        return new Multivector(sig, c);
    }

    @Benchmark
    public double gp() {
        return a.gp(b).get(0);
    }

    @Benchmark
    public double sandwich() {
        return a.sandwich(b).get(0);
    }

    @Benchmark
    public double inverse() {
        return versor.inverse().get(0);
    }
}
