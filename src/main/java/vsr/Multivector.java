package vsr;

/**
 * Immutable multivector: a dense array of 2^n coefficients over the basis
 * blades of a {@link Signature}.
 */
public final class Multivector {
    public final Signature sig;
    private final double[] c;

    private Multivector(Signature sig, double[] c) {
        this.sig = sig;
        this.c = c;
    }

    public static Multivector zero(Signature sig) {
        return new Multivector(sig, new double[sig.size]);
    }

    public static Multivector scalar(Signature sig, double v) {
        double[] c = new double[sig.size];
        c[0] = v;
        return new Multivector(sig, c);
    }

    public static Multivector blade(Signature sig, int mask, double v) {
        double[] c = new double[sig.size];
        c[mask] = v;
        return new Multivector(sig, c);
    }

    /** One coordinate per generator, e.g. vector(sig, x, y, z). */
    public static Multivector vector(Signature sig, double... coords) {
        if (coords.length != sig.n) {
            throw new IllegalArgumentException("expected " + sig.n + " coordinates");
        }
        double[] c = new double[sig.size];
        for (int i = 0; i < sig.n; i++) c[1 << i] = coords[i];
        return new Multivector(sig, c);
    }

    public double get(int mask) {
        return c[mask];
    }

    public double[] toArray() {
        return c.clone();
    }

    public double scalarPart() {
        return c[0];
    }

    public Multivector add(Multivector o) {
        check(o);
        double[] r = c.clone();
        for (int i = 0; i < r.length; i++) r[i] += o.c[i];
        return new Multivector(sig, r);
    }

    public Multivector add(double s) {
        double[] r = c.clone();
        r[0] += s;
        return new Multivector(sig, r);
    }

    public Multivector sub(Multivector o) {
        check(o);
        double[] r = c.clone();
        for (int i = 0; i < r.length; i++) r[i] -= o.c[i];
        return new Multivector(sig, r);
    }

    public Multivector sub(double s) {
        return add(-s);
    }

    public Multivector neg() {
        double[] r = new double[c.length];
        for (int i = 0; i < r.length; i++) r[i] = -c[i];
        return new Multivector(sig, r);
    }

    public Multivector scale(double s) {
        double[] r = new double[c.length];
        for (int i = 0; i < r.length; i++) r[i] = c[i] * s;
        return new Multivector(sig, r);
    }

    public Multivector divide(double s) {
        return scale(1.0 / s);
    }

    /** Geometric product. */
    public Multivector gp(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0) continue;
                int sign = sig.bladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Outer (wedge) product. */
    public Multivector op(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0 || (a & b) != 0) continue;
                int sign = sig.bladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Left contraction: A ⌋ B, selecting blades with basis(A) ⊆ basis(B). */
    public Multivector ip(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0 || (a & b) != a) continue;
                int sign = sig.bladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Regressive product, built from the complement dual. */
    public Multivector vee(Multivector o) {
        return dual().op(o.dual()).dual();
    }

    public Multivector reverse() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            int k = sig.grade(m);
            double s = ((k * (k - 1) / 2) & 1) == 0 ? 1 : -1;
            r[m] = c[m] * s;
        }
        return new Multivector(sig, r);
    }

    /** Clifford conjugate: reversion times grade involution. */
    public Multivector conjugate() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            int k = sig.grade(m);
            double s = ((k * (k + 1) / 2) & 1) == 0 ? 1 : -1;
            r[m] = c[m] * s;
        }
        return new Multivector(sig, r);
    }

    /** Grade involution: (-1)^k on grade-k blades. */
    public Multivector involution() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            r[m] = (sig.grade(m) & 1) == 0 ? c[m] : -c[m];
        }
        return new Multivector(sig, r);
    }

    public Multivector grade(int k) {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            if (sig.grade(m) == k) r[m] = c[m];
        }
        return new Multivector(sig, r);
    }

    public Multivector even() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            if ((sig.grade(m) & 1) == 0) r[m] = c[m];
        }
        return new Multivector(sig, r);
    }

    public Multivector odd() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            if ((sig.grade(m) & 1) != 0) r[m] = c[m];
        }
        return new Multivector(sig, r);
    }

    /** Poincare dual: complement blade with the sign making blade ∧ dual(blade) = ±I. */
    public Multivector dual() {
        int full = sig.fullMask();
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            if (c[m] == 0) continue;
            int comp = full ^ m;
            int swaps = 0;
            for (int i = 0; i < sig.n; i++) {
                if ((m & (1 << i)) != 0) swaps += Integer.bitCount(comp & ((1 << i) - 1));
            }
            r[comp] += ((swaps & 1) == 0 ? 1 : -1) * c[m];
        }
        return new Multivector(sig, r);
    }

    /** Scalar part of A ~A; negative for timelike/null-signature blades. */
    public double norm2() {
        return gp(reverse()).scalarPart();
    }

    public double norm() {
        return Math.sqrt(Math.abs(norm2()));
    }

    public Multivector unit() {
        double n = norm();
        return n == 0 ? this : divide(n);
    }

    /** Inverse for elements with nonzero ~A A. */
    public Multivector inverse() {
        double n2 = norm2();
        if (n2 == 0) throw new ArithmeticException("element is not invertible");
        return reverse().divide(n2);
    }

    /** Sandwich product M x ~M (a transformation by the versor M). */
    public Multivector sandwich(Multivector x) {
        return gp(x).gp(reverse());
    }

    /**
     * Exponential, summed with scaling and squaring. Works for rotors,
     * translators and general motors (including null bivectors).
     */
    public Multivector exp() {
        double max = maxAbs();
        if (max == 0) return scalar(sig, 1);
        int squares = Math.max(0, (int) Math.ceil(Math.log(max) / Math.log(2)));
        Multivector a = divide(1L << squares);
        Multivector term = scalar(sig, 1);
        Multivector sum = scalar(sig, 1);
        for (int k = 1; k <= 32; k++) {
            term = term.gp(a).divide(k);
            sum = sum.add(term);
            if (term.maxAbs() < 1e-18) break;
        }
        for (int i = 0; i < squares; i++) sum = sum.gp(sum);
        return sum;
    }

    public double maxAbs() {
        double m = 0;
        for (double v : c) m = Math.max(m, Math.abs(v));
        return m;
    }

    public boolean isZero(double eps) {
        return maxAbs() < eps;
    }

    public double distanceTo(Multivector o) {
        check(o);
        double m = 0;
        for (int i = 0; i < c.length; i++) m = Math.max(m, Math.abs(c[i] - o.c[i]));
        return m;
    }

    private void check(Multivector o) {
        if (o.sig != sig) throw new IllegalArgumentException("different algebras");
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int m = 0; m < c.length; m++) {
            if (Math.abs(c[m]) < 1e-12) continue;
            if (sb.length() > 0) sb.append(c[m] < 0 ? " - " : " + ");
            else if (c[m] < 0) sb.append('-');
            sb.append(String.format("%.6g", Math.abs(c[m])));
            if (m != 0) sb.append(sig.bladeName(m));
        }
        return sb.length() == 0 ? "0" : sb.toString();
    }
}
