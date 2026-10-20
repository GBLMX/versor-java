package vsr;

/**
 * Immutable multivector: a dense array of 2^n coefficients over the basis
 * blades of a {@link Signature}.
 */
public final class Multivector {
    public final Signature sig;
    private final double[] c;

    Multivector(Signature sig, double[] c) {
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

    /** Direct coefficient array for in-package consumers; must not be mutated. */
    double[] coefficients() {
        return c;
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
        int[][] signs = sig.signTable();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0) continue;
                int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Outer (wedge) product. */
    public Multivector op(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0 || (a & b) != 0) continue;
                int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Left contraction: A ⌋ B, selecting blades with basis(A) ⊆ basis(B). */
    public Multivector ip(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0 || (a & b) != a) continue;
                int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Regressive product, built from the complement dual. */
    public Multivector vee(Multivector o) {
        return dual().op(o.dual()).dual();
    }

    /** Right contraction A ⌊ B: selects blades with basis(B) ⊆ basis(A). */
    public Multivector rightContraction(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0 || (a & b) != b) continue;
                int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                if (sign != 0) r[a ^ b] += sign * ca * cb;
            }
        }
        return new Multivector(sig, r);
    }

    /** Commutator [A, B] = (AB - BA) / 2. */
    public Multivector commutator(Multivector o) {
        return gp(o).sub(o.gp(this)).scale(0.5);
    }

    /** Anticommutator {A, B} = (AB + BA) / 2. */
    public Multivector anticommutator(Multivector o) {
        return gp(o).add(o.gp(this)).scale(0.5);
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

    /**
     * Fast inverse for versors and blades, where {@code A ~A} is the scalar
     * {@link #norm2()}. For arbitrary invertible multivectors use
     * {@link #inverseGeneral()}.
     */
    public Multivector inverse() {
        double n2 = norm2();
        if (n2 == 0) throw new ArithmeticException("element is not invertible");
        return reverse().divide(n2);
    }

    /**
     * Inverse for an arbitrary invertible multivector, including elements
     * whose reverse product is not scalar. Solves the linear system
     * {@code A x = 1} on the coefficient space, costing {@code O(size^3)}
     * scalar operations; {@link #inverse()} is the {@code O(size^2)} fast path
     * for versors.
     */
    public Multivector inverseGeneral() {
        int n = c.length;
        double[][] m = new double[n][n + 1];
        for (int col = 0; col < n; col++) {
            Multivector product = gp(new Multivector(sig, unitMask(col)));
            for (int row = 0; row < n; row++) {
                m[row][col] = product.get(row);
            }
        }
        for (int row = 0; row < n; row++) {
            m[row][n] = row == 0 ? 1 : 0;
        }
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(m[row][col]) > Math.abs(m[pivot][col])) pivot = row;
            }
            if (Math.abs(m[pivot][col]) < 1e-14) {
                throw new ArithmeticException("element is not invertible");
            }
            double[] swap = m[col];
            m[col] = m[pivot];
            m[pivot] = swap;
            double p = m[col][col];
            for (int j = col; j <= n; j++) m[col][j] /= p;
            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double f = m[row][col];
                if (f == 0) continue;
                for (int j = col; j <= n; j++) m[row][j] -= f * m[col][j];
            }
        }
        double[] r = new double[n];
        for (int row = 0; row < n; row++) r[row] = m[row][n];
        return new Multivector(sig, r);
    }

    private double[] unitMask(int mask) {
        double[] v = new double[sig.size];
        v[mask] = 1;
        return v;
    }

    /** Sandwich product M x ~M (a transformation by the versor M). */
    public Multivector sandwich(Multivector x) {
        check(x);
        if (isE3() && onlyEvenGrades() && x.onlyBlades(1, 2, 4)) {
            return sandwichE3Vector(x);
        }
        if (isPga3() && onlyEvenGrades() && x.onlyBlades(7, 11, 13, 14)) {
            return sandwichPgaPoint(x);
        }
        return gp(x).gp(reverse());
    }

    private boolean isE3() {
        return sig.n == 3 && sig.metric[0] == 1 && sig.metric[1] == 1 && sig.metric[2] == 1;
    }

    private boolean isPga3() {
        return sig.n == 4 && sig.metric[0] == 0 && sig.metric[1] == 1
                && sig.metric[2] == 1 && sig.metric[3] == 1;
    }

    /** True when every odd-grade coefficient is zero. */
    private boolean onlyEvenGrades() {
        double[] arr = c;
        for (int m = 1; m < arr.length; m++) {
            if (arr[m] != 0 && (Integer.bitCount(m) & 1) != 0) return false;
        }
        return true;
    }

    private boolean onlyBlades(int m0, int m1, int m2) {
        double[] arr = c;
        for (int m = 0; m < arr.length; m++) {
            if (m != m0 && m != m1 && m != m2 && arr[m] != 0) return false;
        }
        return true;
    }

    private boolean onlyBlades(int m0, int m1, int m2, int m3) {
        double[] arr = c;
        for (int m = 0; m < arr.length; m++) {
            if (m != m0 && m != m1 && m != m2 && m != m3 && arr[m] != 0) return false;
        }
        return true;
    }

    /**
     * Even E3 versor acting on a vector, unrolled from the symbolic sandwich.
     * {@code b1}, {@code b2}, {@code b3} are the coefficients of e12, e13 and e23
     * (masks 3, 5, 6), i.e. the implementation's ascending-blade basis.
     */
    private Multivector sandwichE3Vector(Multivector x) {
        double s = c[0];
        double b1 = c[3];
        double b2 = c[5];
        double b3 = c[6];
        double vx = x.c[1];
        double vy = x.c[2];
        double vz = x.c[4];
        double rx = (s * s - b1 * b1 - b2 * b2 + b3 * b3) * vx
                + 2 * (b1 * (s * vy + b3 * vz) + b2 * (s * vz - b3 * vy));
        double ry = (s * s - b1 * b1 + b2 * b2 - b3 * b3) * vy
                + 2 * (-vx * (b1 * s + b2 * b3) + vz * (b3 * s - b1 * b2));
        double rz = (s * s + b1 * b1 - b2 * b2 - b3 * b3) * vz
                + 2 * (vx * (b1 * b3 - b2 * s) - vy * (b1 * b2 + b3 * s));
        double[] r = new double[8];
        r[1] = rx;
        r[2] = ry;
        r[4] = rz;
        return new Multivector(sig, r);
    }

    /**
     * Even PGA3 motor acting on a point {@code P = p7 e012 + p11 e013 + p13 e023 + p14 e123},
     * unrolled from the symbolic sandwich.
     */
    private Multivector sandwichPgaPoint(Multivector x) {
        double s = c[0];
        double b01 = c[3];
        double b02 = c[5];
        double b12 = c[6];
        double b03 = c[9];
        double b13 = c[10];
        double b23 = c[12];
        double i = c[15];
        double p7 = x.c[7];
        double p11 = x.c[11];
        double p13 = x.c[13];
        double p14 = x.c[14];
        double r7 = (s * s + b12 * b12 - b13 * b13 - b23 * b23) * p7
                + 2 * (b12 * (b13 * p11 + b23 * p13) + s * (b23 * p11 - b13 * p13)
                        + p14 * (i * b12 - b01 * b13 - b02 * b23 + b03 * s));
        double r11 = (s * s - b12 * b12 + b13 * b13 - b23 * b23) * p11
                + 2 * (b12 * (b13 * p7 + s * p13) + b13 * b23 * p13 - s * b23 * p7
                        + p14 * (i * b13 + b01 * b12 - b02 * s - b03 * b23));
        double r13 = (s * s - b12 * b12 - b13 * b13 + b23 * b23) * p13
                + 2 * (b12 * (b23 * p7 - s * p11) + b13 * (b23 * p11 + s * p7)
                        + p14 * (i * b23 + b01 * s + b02 * b12 + b03 * b13));
        double r14 = (s * s + b12 * b12 + b13 * b13 + b23 * b23) * p14;
        double[] r = new double[16];
        r[7] = r7;
        r[11] = r11;
        r[13] = r13;
        r[14] = r14;
        return new Multivector(sig, r);
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
