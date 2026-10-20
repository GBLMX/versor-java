package vsr;

/**
 * Rotor and motor constructors, following versor's {@code Gen} helpers:
 * rotors are exponentials of bivectors, so no angle bookkeeping is hidden
 * in the type.
 */
public final class Gen {
    private Gen() {}

    /** exp(B) for a bivector B (simple or null), closed form. */
    public static Multivector expBivector(Multivector b) {
        double b2 = b.gp(b).scalarPart();
        if (b2 < -1e-14) {
            double r = Math.sqrt(-b2);
            return Multivector.scalar(b.sig, Math.cos(r)).add(b.scale(Math.sin(r) / r));
        }
        if (b2 > 1e-14) {
            double r = Math.sqrt(b2);
            return Multivector.scalar(b.sig, Math.cosh(r)).add(b.scale(Math.sinh(r) / r));
        }
        // null bivector (PGA translators) or zero
        return Multivector.scalar(b.sig, 1).add(b);
    }

    /** Versor-style rotor: exp(B). */
    public static Multivector rotor(Multivector bivector) {
        return expBivector(bivector);
    }

    /** Generator B of a rotor R = exp(B); B is the grade-2 part for simple rotors. */
    public static Multivector log(Multivector rotor) {
        Multivector b = rotor.grade(2);
        double s = rotor.scalarPart();
        double n = b.norm();
        if (n < 1e-14) {
            // identity, or the antipodal rotor -1 whose logarithm is pi times
            // any Euclidean unit bivector (path is ambiguous at 180 degrees)
            return s < 0 ? piBivector(rotor.sig) : b;
        }
        double angle = Math.atan2(n, s);
        return b.scale(angle / n);
    }

    private static Multivector piBivector(Signature sig) {
        for (int m = 0; m < sig.size; m++) {
            if (sig.grade(m) != 2) continue;
            int lo = Integer.numberOfTrailingZeros(m);
            int hi = 31 - Integer.numberOfLeadingZeros(m);
            if (sig.metric[lo] == 1 && sig.metric[hi] == 1) {
                return Multivector.blade(sig, m, Math.PI);
            }
        }
        throw new ArithmeticException("no Euclidean bivector in this algebra");
    }

    /** Shortest-path interpolation of two unit rotors. */
    public static Multivector slerp(Multivector a, Multivector b, double t) {
        Multivector relative = a.reverse().gp(b);
        return a.gp(expBivector(log(relative).scale(t)));
    }

    /** R x ~R */
    public static Multivector apply(Multivector versor, Multivector x) {
        return versor.sandwich(x);
    }
}
