package vsr;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import java.util.Random;

/**
 * Emits a Google-Benchmark-shaped JSON consumable by ga-benchmark's gabm.py.
 * Coverage: Euclidean D=3/4, Homogeneous D=3, Conformal D=3 (all grades).
 * Minkowski is omitted; it is not a signature used by this library.
 */
public final class GaBenchmarkJson {
    private static final double TOL = 1e-8;

    public static void main(String[] args) throws IOException {
        Path out = Path.of(args.length > 0 ? args[0] : "GABenchmark_versor-java.json");
        StringBuilder json = new StringBuilder();
        json.append("{\n  \"context\": {")
                .append("\"library_build_type\": \"release\", ")
                .append("\"num_cpus\": ").append(Runtime.getRuntime().availableProcessors()).append(", ")
                .append("\"mhz_per_cpu\": 0, \"cpu_scaling_enabled\": false, ")
                .append("\"load_avg\": [0.0, 0.0, 0.0]},\n  \"benchmarks\": [\n");
        boolean first = true;
        boolean[] firstRef = {true};
        for (String[] md : new String[][]{{"Euclidean", "3"}, {"Euclidean", "4"},
                {"Homogeneous", "3"}, {"Conformal", "3"}}) {
            String model = md[0];
            int d = Integer.parseInt(md[1]);
            {
                Signature sig = signature(model, d);
                int modelId = switch (model) {
                    case "Conformal" -> 1;
                    case "Euclidean" -> 2;
                    case "Homogeneous" -> 3;
                    default -> 4;
                };
                Random rnd = new Random(1000 + d + model.hashCode());
                for (int grade = 0; grade <= sig.n; grade++) {
                    Multivector a = randomBlade(sig, grade, rnd);
                    Multivector rhs = randomBlade(sig, grade, rnd);
                    for (String op : new String[]{"CliffordConjugation", "Dualization",
                            "GradeInvolution", "Inversion", "Normalization", "Reversion",
                            "SquaredReverseNorm", "UnaryMinus", "UnaryPlus"}) {
                        Multivector arg = a;
                        if (op.equals("Inversion") || op.equals("Normalization")) {
                            arg = invertible(sig, grade, rnd);
                            if (arg == null) continue; // null blades have no inverse
                        }
                        final Multivector argFinal = arg;
                        double ns = time(() -> applyUnary(op, argFinal).get(0));
                        first = emit(json, first, "UnaryOperation_" + op, modelId, d,
                                "Grade:" + grade, ns);
                    }
                    for (String op : new String[]{"Addition", "Subtraction",
                            "GeometricProduct", "OuterProduct", "LeftContraction",
                            "RightContraction", "RegressiveProduct", "CommutatorProduct",
                            "ScalarProduct", "HestenesInnerProduct",
                            "InverseGeometricProduct"}) {
                        Multivector l = a;
                        Multivector r = rhs;
                        if (op.equals("InverseGeometricProduct")) {
                            r = invertible(sig, grade, rnd);
                            if (r == null) continue;
                        }
                        final Multivector rFinal = r;
                        double ns = time(() -> applyBinary(op, l, rFinal).get(0));
                        first = emit(json, first, "BinaryOperation_" + op, modelId, d,
                                "LeftGrade:" + grade + "/RightGrade:" + grade, ns);
                    }
                }
            }
        }
        json.append("\n  ]\n}\n");
        Files.writeString(out, json.toString());
        System.out.println("wrote " + out + " (" + out.toFile().length() + " bytes)");
    }

    /** Short adaptive timing: warm up, then five ~20 ms rounds, keep the minimum. */
    private static double time(GaBench.Op op) {
        for (int i = 0; i < 3; i++) op.apply();
        long t0 = System.nanoTime();
        op.apply();
        double est = Math.max(1, System.nanoTime() - t0);
        long iters = Math.max(3, Math.min(200_000, (long) (20_000_000 / est)));
        for (int i = 0; i < iters; i++) op.apply();
        double best = Double.MAX_VALUE;
        for (int round = 0; round < 5; round++) {
            t0 = System.nanoTime();
            for (long i = 0; i < iters; i++) op.apply();
            best = Math.min(best, (System.nanoTime() - t0) / (double) iters);
        }
        return best;
    }

    private static boolean emit(StringBuilder b, boolean first, String op, int model, int d,
                                String args, double ns) {
        if (!first) b.append(",\n");
        b.append("    {\"name\": \"GABM_").append(op)
                .append("/SystemName:Linux/SystemVersion:unknown/CompilerID:GNU/CompilerVersion:25")
                .append("/Solution:versor-java/Model:").append(model)
                .append("/D:").append(d).append('/').append(args)
                .append("\", \"iterations\": 1, \"cpu_time\": ")
                .append(String.format(Locale.ROOT, "%.6f", ns / 1e6))
                .append(", \"time_unit\": \"ms\", \"aggregate_name\": \"mean\"}");
        return false;
    }

    private static Signature signature(String model, int d) {
        return switch (model) {
            case "Euclidean" -> Signature.euclidean(d);
            case "Homogeneous" -> {
                int[] m = new int[d + 1];
                m[0] = 0;
                for (int i = 1; i <= d; i++) m[i] = 1;
                yield Signature.of(m);
            }
            case "Conformal" -> {
                int[] m = new int[d + 2];
                for (int i = 0; i < d + 1; i++) m[i] = 1;
                m[d + 1] = -1;
                yield Signature.of(m);
            }
            default -> throw new IllegalArgumentException(model);
        };
    }

    private static Multivector randomBlade(Signature sig, int grade, Random rnd) {
        Multivector m = Multivector.scalar(sig, rnd.nextDouble() * 2 - 1);
        for (int i = 0; i < grade; i++) {
            double[] coords = new double[sig.n];
            coords[i] = 1;
            for (int j = grade; j < sig.n; j++) coords[j] = rnd.nextDouble() * 2 - 1;
            m = m.op(Multivector.vector(sig, coords));
        }
        return m;
    }

    /** Returns null when every blade of this grade is null (e.g. contains e0). */
    private static Multivector invertible(Signature sig, int grade, Random rnd) {
        for (int attempt = 0; attempt < 50; attempt++) {
            Multivector m = randomBlade(sig, grade, rnd);
            if (Math.abs(m.norm2()) > TOL) return m;
        }
        return null;
    }

    private static Multivector applyUnary(String op, Multivector a) {
        return switch (op) {
            case "CliffordConjugation" -> a.conjugate();
            case "Dualization" -> a.dual();
            case "GradeInvolution" -> a.involution();
            case "Inversion" -> a.inverse();
            case "Normalization" -> a.unit();
            case "Reversion" -> a.reverse();
            case "SquaredReverseNorm" -> Multivector.scalar(a.sig, a.norm2());
            case "UnaryMinus" -> a.neg();
            case "UnaryPlus" -> a;
            default -> throw new IllegalArgumentException(op);
        };
    }

    private static Multivector applyBinary(String op, Multivector a, Multivector b) {
        return switch (op) {
            case "Addition" -> a.add(b);
            case "Subtraction" -> a.sub(b);
            case "GeometricProduct" -> a.gp(b);
            case "OuterProduct" -> a.op(b);
            case "LeftContraction" -> a.ip(b);
            case "RightContraction" -> a.rightContraction(b);
            case "RegressiveProduct" -> a.vee(b);
            case "CommutatorProduct" -> a.commutator(b);
            case "ScalarProduct" -> Multivector.scalar(a.sig, a.gp(b).scalarPart());
            case "HestenesInnerProduct" -> hestenes(a, b);
            case "InverseGeometricProduct" -> a.gp(b.inverse());
            default -> throw new IllegalArgumentException(op);
        };
    }

    private static Multivector hestenes(Multivector a, Multivector b) {
        Multivector sum = Multivector.zero(a.sig);
        for (int r = 1; r <= a.sig.n; r++) {
            Multivector ar = a.grade(r);
            if (ar.maxAbs() == 0) continue;
            for (int s = 1; s <= b.sig.n; s++) {
                Multivector bs = b.grade(s);
                if (bs.maxAbs() == 0) continue;
                sum = sum.add(ar.gp(bs).grade(Math.abs(r - s)));
            }
        }
        return sum;
    }
}
