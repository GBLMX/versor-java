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

    /**
     * Right division {@code A / B = A B^-1}. The divisor must be a versor or
     * blade for the fast inverse; for an arbitrary invertible divisor use
     * {@code a.gp(b.inverseGeneral())}.
     */
    public Multivector divide(Multivector divisor) {
        return gp(divisor.inverse());
    }

    /** Scalar product {@code A . B = scalar(AB)}, the grade-0 part of the gp. */
    public double scalarProduct(Multivector o) {
        return gp(o).scalarPart();
    }

    /** Readable alias for {@link #gp(Multivector)}. */
    public Multivector mul(Multivector o) {
        return gp(o);
    }

    /** Scalar multiplication, matching the {@link #mul(Multivector)} spelling. */
    public Multivector mul(double s) {
        return scale(s);
    }

    /** Geometric product. See {@link #mul(Multivector)} for the readable alias. */
    public Multivector gp(Multivector o) {
        double[] r = new double[c.length];
        gpAccumulate(o, r);
        return new Multivector(sig, r);
    }

    /**
     * Accumulates {@code this * o} into {@code out}; the caller must supply a
     * zeroed coefficient array of the right length. Internal scheduling hook
     * used by matrix products to avoid one intermediate multivector per term.
     */
    void gpAccumulate(Multivector o, double[] out) {
        gpAccumulate(o, out, 1.0);
    }

    /**
     * Accumulates {@code scale * (this * o)} into {@code out}. A scale of
     * {@code -1} is used by matrix elimination to avoid a temporary product
     * and subtraction.
     */
    void gpAccumulate(Multivector o, double[] out, double scale) {
        if (scale == 0) return;
        check(o);
        if (out.length != c.length) {
            throw new IllegalArgumentException("output length mismatch");
        }
        int[][] signs = sig.signTable();
        int aCount = nonZeroCount(c);
        int bCount = nonZeroCount(o.c);
        if ((long) aCount * bCount * 4 < (long) c.length * c.length) {
            int[] ai = nonZeroIndices(c, aCount);
            int[] bi = nonZeroIndices(o.c, bCount);
            for (int ia : ai) {
                double ca = scale * c[ia];
                int[] row = signs != null ? signs[ia] : null;
                for (int ib : bi) {
                    int sign = row != null ? row[ib] : sig.computeBladeProductSign(ia, ib);
                    if (sign != 0) out[ia ^ ib] += sign * ca * o.c[ib];
                }
            }
            return;
        }
        for (int a = 0; a < c.length; a++) {
            double ca = scale * c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = 0; b < c.length; b++) {
                double cb = o.c[b];
                if (cb == 0) continue;
                int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                if (sign != 0) out[a ^ b] += sign * ca * cb;
            }
        }
    }

    private static int nonZeroCount(double[] values) {
        int nz = 0;
        for (double v : values) {
            if (v != 0) nz++;
        }
        return nz;
    }

    private static int[] nonZeroIndices(double[] values, int count) {
        int[] idx = new int[count];
        int p = 0;
        for (int i = 0; i < values.length; i++) {
            if (values[i] != 0) idx[p++] = i;
        }
        return idx;
    }

    /** Outer (wedge) product. Enumerates disjoint blades: O(3^n) pairs. See {@link #wedge(Multivector)}. */
    public Multivector op(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        int full = sig.fullMask();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            int comp = full ^ a;
            for (int b = comp; ; b = (b - 1) & comp) {
                double cb = o.c[b];
                if (cb != 0) {
                    int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                    if (sign != 0) r[a ^ b] += sign * ca * cb;
                }
                if (b == 0) break;
            }
        }
        return new Multivector(sig, r);
    }

    /** Readable alias for {@link #op(Multivector)}. */
    public Multivector wedge(Multivector o) {
        return op(o);
    }

    /** Left contraction A ⌋ B. Enumerates supersets of A: O(3^n) pairs. See {@link #contract(Multivector)}. */
    public Multivector ip(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        int full = sig.fullMask();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            int comp = full ^ a;
            for (int sub = comp; ; sub = (sub - 1) & comp) {
                int b = a | sub;
                double cb = o.c[b];
                if (cb != 0) {
                    int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                    if (sign != 0) r[sub] += sign * ca * cb;
                }
                if (sub == 0) break;
            }
        }
        return new Multivector(sig, r);
    }

    /** Readable alias for {@link #ip(Multivector)}. */
    public Multivector contract(Multivector o) {
        return ip(o);
    }

    /** Regressive product, built from the complement dual. See {@link #regressive(Multivector)}. */
    public Multivector vee(Multivector o) {
        return dual().op(o.dual()).dual();
    }

    /** Readable alias for {@link #vee(Multivector)}. */
    public Multivector regressive(Multivector o) {
        return vee(o);
    }

    /** Right contraction A ⌊ B. Enumerates subsets of A: O(3^n) pairs. */
    public Multivector rightContraction(Multivector o) {
        check(o);
        double[] r = new double[c.length];
        int[][] signs = sig.signTable();
        for (int a = 0; a < c.length; a++) {
            double ca = c[a];
            if (ca == 0) continue;
            int[] row = signs != null ? signs[a] : null;
            for (int b = a; ; b = (b - 1) & a) {
                double cb = o.c[b];
                if (cb != 0) {
                    int sign = row != null ? row[b] : sig.computeBladeProductSign(a, b);
                    if (sign != 0) r[a ^ b] += sign * ca * cb;
                }
                if (b == 0) break;
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
            r[m] = c[m] * Signature.reverseSign(sig.grade(m));
        }
        return new Multivector(sig, r);
    }

    /** Clifford conjugate: reversion times grade involution. */
    public Multivector conjugate() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            r[m] = c[m] * Signature.conjugateSign(sig.grade(m));
        }
        return new Multivector(sig, r);
    }

    /** Grade involution: (-1)^k on grade-k blades. */
    public Multivector involution() {
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            r[m] = c[m] * Signature.involutionSign(sig.grade(m));
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
            r[full ^ m] += sig.dualSign(m) * c[m];
        }
        return new Multivector(sig, r);
    }

    /**
     * Inverse of {@link #dual()} under the Poincare complement convention:
     * {@code a.dual().undual() == a} and {@code a.undual().dual() == a}.
     *
     * <p>This is not the C++ versor Hodge undual ({@code A * I}); for a
     * non-degenerate algebra the two differ by a sign convention that depends
     * on the grade and dimension.
     */
    public Multivector undual() {
        int full = sig.fullMask();
        double[] r = new double[c.length];
        for (int m = 0; m < c.length; m++) {
            int comp = full ^ m;
            if (c[comp] == 0) continue;
            int k = Integer.bitCount(m);
            int sign = sig.dualSign(comp) * (((k * (sig.n - k)) & 1) == 0 ? 1 : -1);
            r[m] += sign * c[comp];
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
     * Result of {@link #inverseGeneralResult()}.
     *
     * @param value            the computed inverse
     * @param residual         {@code max |A x - 1|} after iterative refinement
     * @param relativeResidual residual divided by {@code ||A||_inf ||x||_inf + 1}
     * @param condition        infinity-norm estimate {@code ||A||_inf ||A^-1||_inf}
     */
    public record GeneralInverse(Multivector value, double residual,
                                 double relativeResidual, double condition) {}

    /** Default maximum relative backward error accepted by {@link #inverseGeneral()}. */
    public static final double INVERSE_TOLERANCE = 1e-12;

    /** Default condition-number ceiling accepted by {@link #inverseGeneral()}. */
    public static final double INVERSE_CONDITION_LIMIT = 1e12;

    /**
     * Inverse for an arbitrary invertible multivector, including elements
     * whose reverse product is not scalar. Solves the linear system
     * {@code A x = 1} on the coefficient space, costing {@code O(size^3)}
     * scalar operations; {@link #inverse()} is the {@code O(size^2)} fast path
     * for versors.
     *
     * <p>The coefficient matrix is equilibrated before elimination and the
     * result is checked with a backward-error bound. An
     * {@link ArithmeticException} is thrown when the relative residual exceeds
     * {@link #INVERSE_TOLERANCE}; use {@link #inverseGeneral(double)} to choose
     * another tolerance, or {@link #inverseGeneralResult()} to inspect an
     * ill-conditioned result instead of rejecting it. The no-argument form also
     * rejects condition estimates above {@link #INVERSE_CONDITION_LIMIT}.
     */
    public Multivector inverseGeneral() {
        GeneralInverse result = inverseGeneral(INVERSE_TOLERANCE);
        if (!(result.condition() <= INVERSE_CONDITION_LIMIT)) {
            throw new ArithmeticException("ill-conditioned inverse: condition "
                    + result.condition() + " > " + INVERSE_CONDITION_LIMIT);
        }
        return result.value();
    }

    /** Checked inversion with an explicit relative backward-error tolerance. */
    public GeneralInverse inverseGeneral(double tolerance) {
        GeneralInverse result = inverseGeneralResult();
        if (!(result.relativeResidual() <= tolerance)) {
            throw new ArithmeticException("ill-conditioned inverse: relative residual "
                    + result.relativeResidual() + " > " + tolerance
                    + ", condition ~ " + result.condition());
        }
        return result;
    }

    /**
     * Computes the general inverse and its error estimate without rejecting a
     * large residual. Exact singularities still throw {@link ArithmeticException}.
     */
    public GeneralInverse inverseGeneralResult() {
        int n = c.length;
        double[][] a = new double[n][n];
        for (int col = 0; col < n; col++) {
            Multivector product = gp(new Multivector(sig, unitMask(col)));
            for (int row = 0; row < n; row++) {
                a[row][col] = product.get(row);
            }
        }

        double[] rowScale = new double[n];
        for (int row = 0; row < n; row++) {
            double max = 0;
            for (int col = 0; col < n; col++) max = Math.max(max, Math.abs(a[row][col]));
            if (max == 0) throw new ArithmeticException("element is not invertible");
            rowScale[row] = max;
        }
        double[] colScale = new double[n];
        for (int col = 0; col < n; col++) {
            double max = 0;
            for (int row = 0; row < n; row++) {
                max = Math.max(max, Math.abs(a[row][col]) / rowScale[row]);
            }
            if (max == 0) throw new ArithmeticException("element is not invertible");
            colScale[col] = max;
        }

        double[][] s = new double[n][n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                s[row][col] = a[row][col] / (rowScale[row] * colScale[col]);
            }
        }
        double[][] sInv = invertEquilibrated(s);
        double[][] inv = new double[n][n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                inv[row][col] = sInv[row][col] / (colScale[row] * rowScale[col]);
            }
        }

        double[] x = new double[n];
        for (int row = 0; row < n; row++) x[row] = inv[row][0];
        double residual = inverseResidual(a, x);
        for (int step = 0; step < 2 && residual > 1e-14; step++) {
            double[] correction = matVec(inv, inverseErrorVector(a, x));
            double[] candidate = x.clone();
            for (int row = 0; row < n; row++) candidate[row] -= correction[row];
            double candidateResidual = inverseResidual(a, candidate);
            if (!(candidateResidual < residual)) break;
            x = candidate;
            residual = candidateResidual;
        }

        double normA = maxRowSum(a);
        double normInv = maxRowSum(inv);
        double condition = normA * normInv;
        double scale = normA * maxAbs(x) + 1;
        double relativeResidual = residual / scale;
        if (!Double.isFinite(relativeResidual) || !Double.isFinite(condition)) {
            throw new ArithmeticException("inverse is not numerically finite");
        }
        return new GeneralInverse(new Multivector(sig, x), residual, relativeResidual, condition);
    }

    /** Gauss-Jordan inverse of an already equilibrated coefficient matrix. */
    private static double[][] invertEquilibrated(double[][] a) {
        int n = a.length;
        double[][] m = new double[n][2 * n];
        for (int row = 0; row < n; row++) {
            System.arraycopy(a[row], 0, m[row], 0, n);
            m[row][n + row] = 1;
        }
        double pivotTolerance = Math.ulp(1.0) * n;
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(m[row][col]) > Math.abs(m[pivot][col])) pivot = row;
            }
            if (!(Math.abs(m[pivot][col]) > pivotTolerance)) {
                throw new ArithmeticException("element is not invertible");
            }
            double[] swap = m[col];
            m[col] = m[pivot];
            m[pivot] = swap;
            double p = m[col][col];
            for (int j = col; j < 2 * n; j++) m[col][j] /= p;
            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double f = m[row][col];
                if (f == 0) continue;
                for (int j = col; j < 2 * n; j++) m[row][j] -= f * m[col][j];
            }
        }
        double[][] inv = new double[n][n];
        for (int row = 0; row < n; row++) {
            System.arraycopy(m[row], n, inv[row], 0, n);
        }
        return inv;
    }

    private double inverseResidual(double[][] a, double[] x) {
        return maxAbs(inverseErrorVector(a, x));
    }

    private double[] inverseErrorVector(double[][] a, double[] x) {
        Multivector ax = gp(new Multivector(sig, x.clone()));
        double[] error = new double[a.length];
        for (int row = 0; row < a.length; row++) error[row] = ax.get(row);
        error[0] -= 1;
        return error;
    }

    private static double[] matVec(double[][] matrix, double[] vector) {
        double[] out = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++) {
            double sum = 0;
            for (int col = 0; col < matrix.length; col++) sum += matrix[row][col] * vector[col];
            out[row] = sum;
        }
        return out;
    }

    private static double maxRowSum(double[][] matrix) {
        double max = 0;
        for (double[] row : matrix) {
            double sum = 0;
            for (double value : row) sum += Math.abs(value);
            max = Math.max(max, sum);
        }
        return max;
    }

    private static double maxAbs(double[] values) {
        double max = 0;
        for (double value : values) max = Math.max(max, Math.abs(value));
        return max;
    }

    private double[] unitMask(int mask) {
        double[] v = new double[sig.size];
        v[mask] = 1;
        return v;
    }

    /** Sandwich product M x ~M (a transformation by the versor M). See {@link #apply(Multivector)}. */
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

    /** Readable alias for {@link #sandwich(Multivector)}. */
    public Multivector apply(Multivector x) {
        return sandwich(x);
    }

    /**
     * Reflection of this by an odd unit versor: {@code b x^ b~}, where
     * {@code x^} is the grade involution. This is the C++ versor
     * {@code x.reflect(b)} convention; {@code b} must be unit.
     */
    public Multivector reflect(Multivector reflector) {
        check(reflector);
        return reflector.gp(involution()).gp(reflector.reverse());
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

    /** True when every coefficient is finite. */
    public boolean isFinite() {
        for (double v : c) {
            if (!Double.isFinite(v)) return false;
        }
        return true;
    }

    /** True when only the scalar blade is non-zero within a relative tolerance. */
    public boolean isScalar(double eps) {
        if (!isFinite()) return false;
        double limit = eps * maxAbs();
        for (int m = 1; m < c.length; m++) {
            if (Math.abs(c[m]) > limit) return false;
        }
        return true;
    }

    /** True when all non-zero coefficients have the same grade. */
    public boolean isHomogeneous(double eps) {
        if (!isFinite()) return false;
        double limit = eps * maxAbs();
        int grade = -1;
        for (int m = 0; m < c.length; m++) {
            if (Math.abs(c[m]) <= limit) continue;
            int k = Integer.bitCount(m);
            if (grade < 0) grade = k;
            else if (grade != k) return false;
        }
        return true;
    }

    /**
     * True when this is a simple blade: homogeneous, and the map
     * {@code v -> v ⌋ A} has rank equal to the grade. Vectors and scalars are
     * blades. This covers grades above two, where {@code A ^ A = 0} alone is
     * not sufficient.
     */
    public boolean isBlade(double eps) {
        if (!isFinite() || !isHomogeneous(eps)) return false;
        int grade = -1;
        double limit = eps * maxAbs();
        for (int m = 0; m < c.length; m++) {
            if (Math.abs(c[m]) > limit) {
                grade = Integer.bitCount(m);
                break;
            }
        }
        if (grade <= 1) return true;
        if (grade > sig.n) return false;

        double[][] image = new double[sig.n][sig.size];
        double scale = 0;
        for (int i = 0; i < sig.n; i++) {
            Multivector contracted = Multivector.blade(sig, 1 << i, 1).ip(this);
            image[i] = contracted.coefficients();
            for (double v : image[i]) scale = Math.max(scale, Math.abs(v));
        }
        return rank(image, eps * scale) == grade;
    }

    private static int rank(double[][] matrix, double tolerance) {
        int rows = matrix.length;
        int cols = rows == 0 ? 0 : matrix[0].length;
        int rank = 0;
        int pivotRow = 0;
        for (int col = 0; col < cols && pivotRow < rows; col++) {
            int best = pivotRow;
            double bestAbs = Math.abs(matrix[best][col]);
            for (int r = pivotRow + 1; r < rows; r++) {
                double v = Math.abs(matrix[r][col]);
                if (v > bestAbs) {
                    best = r;
                    bestAbs = v;
                }
            }
            if (bestAbs <= tolerance) continue;
            if (best != pivotRow) {
                double[] tmp = matrix[pivotRow];
                matrix[pivotRow] = matrix[best];
                matrix[best] = tmp;
            }
            double pivot = matrix[pivotRow][col];
            for (int r = pivotRow + 1; r < rows; r++) {
                double factor = matrix[r][col] / pivot;
                if (factor == 0) continue;
                for (int k = col; k < cols; k++) {
                    matrix[r][k] -= factor * matrix[pivotRow][k];
                }
            }
            pivotRow++;
            rank++;
        }
        return rank;
    }

    /**
     * True when {@code A ~A} is a non-zero scalar, i.e. this is an invertible
     * versor. This covers vectors, rotors, translators and motors.
     */
    public boolean isVersor(double eps) {
        if (!isFinite()) return false;
        Multivector square = gp(reverse());
        if (!square.isScalar(eps)) return false;
        return Math.abs(square.scalarPart()) > eps * square.maxAbs();
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
