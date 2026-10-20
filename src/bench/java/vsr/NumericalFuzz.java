package vsr;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.joml.Quaternionf;
import org.joml.Vector3f;

/** Numerical fuzz: dense vs sparse, ordering, conditioning, degeneracy, float refs. */
public final class NumericalFuzz {
    public static void main(String[] args) {
        Random rnd = new Random(20261021);
        System.out.println("## dense vs sparse");
        for (int n = 4; n <= 10; n++) {
            Signature sig = Signature.euclidean(n);
            double maxGp = 0;
            double maxSand = 0;
            int failures = 0;
            for (int trial = 0; trial < 100; trial++) {
                int k = 1 + rnd.nextInt(Math.min(64, sig.size));
                SparseMultivector a = uniform(sig, k, rnd);
                SparseMultivector b = uniform(sig, k, rnd);
                Multivector da = a.toDense();
                Multivector db = b.toDense();
                double scale = Math.max(1e-30, da.maxAbs() * db.maxAbs());
                double e1 = a.gp(b).toDense().distanceTo(da.gp(db)) / scale;
                double e2 = a.sandwich(b).toDense().distanceTo(da.sandwich(db)) / scale;
                if (!(e1 < 1e-12) || !(e2 < 1e-10)) failures++;
                maxGp = Math.max(maxGp, e1);
                maxSand = Math.max(maxSand, e2);
            }
            System.out.printf("FUZZ,denseSparse,%d,gp,%.2e,%d%n", n, maxGp, failures);
            System.out.printf("FUZZ,denseSparse,%d,sandwich,%.2e,%d%n", n, maxSand, failures);
        }

        System.out.println("## accumulation order (input order shuffled)");
        for (int n : new int[]{8, 10, 12}) {
            Signature sig = Signature.euclidean(n);
            int[] masks = new int[64];
            double[] vals = new double[64];
            for (int i = 0; i < 64; i++) {
                masks[i] = rnd.nextInt(sig.size);
                vals[i] = rnd.nextDouble() * 2 - 1;
            }
            SparseMultivector a = SparseMultivector.of(sig, masks, vals);
            SparseMultivector b = SparseMultivector.of(sig, masks.clone(), vals.clone());
            List<Integer> order = new ArrayList<>();
            for (int i = 0; i < masks.length; i++) order.add(i);
            Collections.shuffle(order, rnd);
            int[] sm = new int[64];
            double[] sv = new double[64];
            for (int i = 0; i < 64; i++) {
                sm[i] = masks[order.get(i)];
                sv[i] = vals[order.get(i)];
            }
            SparseMultivector shuffled = SparseMultivector.of(sig, sm, sv);
            double diff = a.gp(b).toDense().distanceTo(shuffled.gp(b).toDense());
            System.out.printf("FUZZ,order,%d,gp,%.2e,0%n", n, diff);
        }

        System.out.println("## inverseGeneral conditioning (Cl(4), singular base + eps*noise)");
        Signature sig = Signature.euclidean(4);
        for (double eps : new double[]{1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12}) {
            double worst = 0;
            double worstRelative = 0;
            double worstCondition = 0;
            int throwsCount = 0;
            for (int trial = 0; trial < 20; trial++) {
                Multivector a = Multivector.scalar(sig, 0.5)
                        .add(Multivector.blade(sig, 1, 0.5)); // idempotent, norm2 = 0
                for (int m = 1; m < sig.size; m++) {
                    a = a.add(Multivector.blade(sig, m, eps * (rnd.nextDouble() * 2 - 1)));
                }
                Multivector.GeneralInverse result = a.inverseGeneralResult();
                worst = Math.max(worst, result.residual());
                worstRelative = Math.max(worstRelative, result.relativeResidual());
                worstCondition = Math.max(worstCondition, result.condition());
                try {
                    a.inverseGeneral();
                } catch (ArithmeticException e) {
                    throwsCount++;
                }
            }
            System.out.printf("COND,4,%.0e,%.2e,%d,%.2e,%.2e%n",
                    eps, worst, throwsCount, worstCondition, worstRelative);
        }

        System.out.println("## degenerate signature and zero norm");
        Signature pga = Signature.pga3();
        Multivector e0 = Multivector.blade(pga, 1, 1);
        System.out.printf("DEGEN,nullSquare,%s,0%n",
                e0.gp(e0).maxAbs() == 0 ? "zero" : "nonzero");
        System.out.printf("DEGEN,nullNorm,%.1f,0%n", e0.norm2());
        boolean threw = false;
        try {
            e0.inverse();
        } catch (ArithmeticException e) {
            threw = true;
        }
        System.out.printf("DEGEN,nullInverse,%s,0%n", threw ? "throws" : "wrong");

        System.out.println("## double summation vs BigDecimal reference (Cl(8))");
        Signature sig8 = Signature.euclidean(8);
        double worstDense = 0;
        for (int trial = 0; trial < 20; trial++) {
            Multivector a = big(sig8, rnd);
            Multivector b = big(sig8, rnd);
            Multivector p = a.gp(b);
            for (int mask = 0; mask < sig8.size; mask++) {
                BigDecimal ref = refCoefficient(a, b, mask);
                double exact = ref.doubleValue();
                double got = p.get(mask);
                double err = Math.abs(got - exact) / Math.max(1e-30, Math.abs(exact));
                worstDense = Math.max(worstDense, err);
            }
        }
        System.out.printf("SUM,8,double,%.2e,0%n", worstDense);

        System.out.println("## E3 double rotor vs JOML float quaternion");
        double maxErr = 0;
        for (int i = 0; i < 10_000; i++) {
            double ax = rnd.nextDouble() * 2 - 1;
            double ay = rnd.nextDouble() * 2 - 1;
            double az = rnd.nextDouble() * 2 - 1;
            double angle = rnd.nextDouble() * Math.PI * 2;
            Multivector rotor = E3.rotor(ax, ay, az, angle);
            Quaternionf q = new Quaternionf().rotateAxis((float) angle,
                    (float) ax, (float) ay, (float) az).normalize();
            Vector3f v = new Vector3f((float) (rnd.nextDouble() * 4 - 2),
                    (float) (rnd.nextDouble() * 4 - 2), (float) (rnd.nextDouble() * 4 - 2));
            Multivector r = rotor.sandwich(E3.vec(v.x, v.y, v.z));
            Vector3f qv = q.transform(new Vector3f(v));
            double err = Math.sqrt(Math.pow(r.get(1) - qv.x, 2)
                    + Math.pow(r.get(2) - qv.y, 2) + Math.pow(r.get(4) - qv.z, 2));
            maxErr = Math.max(maxErr, err);
        }
        System.out.printf("JOML,e3FloatVsDouble,%.2e,0%n", maxErr);
    }

    private static SparseMultivector uniform(Signature sig, int k, Random rnd) {
        int[] masks = new int[k];
        double[] vals = new double[k];
        for (int i = 0; i < k; i++) {
            masks[i] = rnd.nextInt(sig.size);
            vals[i] = rnd.nextDouble() * 2 - 1;
        }
        return SparseMultivector.of(sig, masks, vals);
    }

    private static Multivector big(Signature sig, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < c.length; i++) c[i] = (rnd.nextDouble() * 2 - 1) * 1e8;
        return new Multivector(sig, c);
    }

    private static BigDecimal refCoefficient(Multivector a, Multivector b, int mask) {
        BigDecimal sum = BigDecimal.ZERO;
        MathContext mc = MathContext.DECIMAL128;
        for (int x = 0; x < a.sig.size; x++) {
            double ca = a.get(x);
            if (ca == 0) continue;
            for (int y = 0; y < b.sig.size; y++) {
                int sign = a.sig.bladeProductSign(x, y);
                if (sign == 0 || (x ^ y) != mask) continue;
                double cb = b.get(y);
                if (cb == 0) continue;
                BigDecimal term = BigDecimal.valueOf(ca).multiply(BigDecimal.valueOf(cb), mc)
                        .multiply(BigDecimal.valueOf(sign), mc);
                sum = sum.add(term, mc);
            }
        }
        return sum;
    }
}
