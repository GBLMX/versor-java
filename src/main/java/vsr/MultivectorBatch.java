package vsr;

/**
 * Structure-of-arrays batch of multivectors sharing one signature.
 *
 * <p>Coefficients are stored per blade ({@code blades[mask][i]}), so the inner
 * loops run over contiguous runs of one coefficient and avoid the per-element
 * object and array allocation of {@code Multivector[]}. Batch arithmetic is
 * element-wise: element {@code i} of the result only depends on element
 * {@code i} of the inputs.
 */
public final class MultivectorBatch {
    public final Signature sig;
    private final int count;
    private final double[][] blades;

    public MultivectorBatch(Signature sig, int count) {
        if (count < 0) throw new IllegalArgumentException("negative count");
        this.sig = sig;
        this.count = count;
        this.blades = new double[sig.size][count];
    }

    private MultivectorBatch(Signature sig, int count, double[][] blades) {
        this.sig = sig;
        this.count = count;
        this.blades = blades;
    }

    /** Copies the given multivectors into a batch. */
    public static MultivectorBatch of(Multivector... values) {
        if (values.length == 0) {
            throw new IllegalArgumentException("empty batch");
        }
        Signature sig = values[0].sig;
        MultivectorBatch batch = new MultivectorBatch(sig, values.length);
        for (int i = 0; i < values.length; i++) {
            batch.set(i, values[i]);
        }
        return batch;
    }

    public int size() {
        return count;
    }

    public Multivector get(int i) {
        checkIndex(i);
        double[] c = new double[sig.size];
        for (int m = 0; m < c.length; m++) {
            c[m] = blades[m][i];
        }
        return new Multivector(sig, c);
    }

    public void set(int i, Multivector value) {
        checkIndex(i);
        checkSignature(value);
        for (int m = 0; m < blades.length; m++) {
            blades[m][i] = value.get(m);
        }
    }

    public double get(int mask, int i) {
        checkIndex(i);
        return blades[mask][i];
    }

    public void set(int mask, int i, double value) {
        checkIndex(i);
        blades[mask][i] = value;
    }

    /** Copies elements {@code [from, to)} into a new batch. */
    public MultivectorBatch slice(int from, int to) {
        if (from < 0 || to > count || from > to) {
            throw new IndexOutOfBoundsException("slice " + from + ".." + to
                    + " outside 0.." + count);
        }
        MultivectorBatch out = new MultivectorBatch(sig, to - from);
        for (int m = 0; m < blades.length; m++) {
            System.arraycopy(blades[m], from, out.blades[m], 0, to - from);
        }
        return out;
    }

    public MultivectorBatch copy() {
        double[][] copy = new double[blades.length][count];
        for (int m = 0; m < blades.length; m++) {
            System.arraycopy(blades[m], 0, copy[m], 0, count);
        }
        return new MultivectorBatch(sig, count, copy);
    }

    /** {@code left * element} for every element in the batch. */
    public MultivectorBatch gp(Multivector left) {
        checkSignature(left);
        double[] la = left.coefficients();
        int size = sig.size;
        boolean[] active = new boolean[size];
        for (int b = 0; b < size; b++) active[b] = hasNonzero(blades[b]);
        double[][] r = new double[size][count];
        for (int a = 0; a < size; a++) {
            double ca = la[a];
            if (ca == 0) continue;
            for (int b = 0; b < size; b++) {
                if (!active[b]) continue;
                int sign = sig.bladeProductSign(a, b);
                if (sign == 0) continue;
                double factor = ca * sign;
                double[] src = blades[b];
                double[] dst = r[a ^ b];
                for (int i = 0; i < count; i++) {
                    dst[i] += factor * src[i];
                }
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    /** {@code element * right} for every element in the batch. */
    public MultivectorBatch gpRight(Multivector right) {
        checkSignature(right);
        double[] ra = right.coefficients();
        int size = sig.size;
        double[][] r = new double[size][count];
        for (int a = 0; a < size; a++) {
            double[] src = blades[a];
            if (!hasNonzero(src)) continue;
            for (int b = 0; b < size; b++) {
                double cb = ra[b];
                if (cb == 0) continue;
                int sign = sig.bladeProductSign(a, b);
                if (sign == 0) continue;
                double factor = cb * sign;
                double[] dst = r[a ^ b];
                for (int i = 0; i < count; i++) {
                    dst[i] += factor * src[i];
                }
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    /**
     * {@code versor * element * ~versor} for every element in the batch.
     *
     * <p>E3 vectors and PGA3 trivectors take the same unrolled formulas as
     * {@link Multivector#sandwich(Multivector)}; the coefficients are read and
     * written in place, so the batch path allocates only its result arrays.
     * Keep these formulas in sync with the scalar fast paths.
     */
    public MultivectorBatch sandwich(Multivector versor) {
        checkSignature(versor);
        if (isE3() && evenOnly(versor) && onlyBlades(1, 2, 4)) {
            return sandwichE3(versor, new double[sig.size][count]);
        }
        if (isPga3() && evenOnly(versor) && onlyBlades(7, 11, 13, 14)) {
            return sandwichPga(versor, new double[sig.size][count]);
        }
        return gp(versor).gpRight(versor.reverse());
    }

    /**
     * Writes {@code versor * element * ~versor} into {@code out} and returns it,
     * reusing its coefficient arrays. The caller owns {@code out}; it must have
     * the same signature and size and must not alias this batch.
     */
    public MultivectorBatch sandwich(Multivector versor, MultivectorBatch out) {
        checkSignature(versor);
        if (out.sig != sig || out.count != count) {
            throw new IllegalArgumentException("batch shape mismatch");
        }
        out.clear();
        if (isE3() && evenOnly(versor) && onlyBlades(1, 2, 4)) {
            return sandwichE3(versor, out.blades);
        }
        if (isPga3() && evenOnly(versor) && onlyBlades(7, 11, 13, 14)) {
            return sandwichPga(versor, out.blades);
        }
        double[] m = versor.coefficients();
        double[] mr = versor.reverse().coefficients();
        double[][] r = out.blades;
        for (int a = 0; a < m.length; a++) {
            if (m[a] == 0) continue;
            for (int c = 0; c < mr.length; c++) {
                if (mr[c] == 0) continue;
                for (int b = 0; b < blades.length; b++) {
                    if (!hasNonzero(blades[b])) continue;
                    int mid = a ^ b;
                    int sign = sig.bladeProductSign(a, b);
                    if (sign == 0) continue;
                    int sign2 = sig.bladeProductSign(mid, c);
                    if (sign2 == 0) continue;
                    double factor = m[a] * mr[c] * sign * sign2;
                    double[] src = blades[b];
                    double[] dst = r[mid ^ c];
                    for (int i = 0; i < count; i++) dst[i] += factor * src[i];
                }
            }
        }
        return out;
    }

    /** {@code left ^ element} for every element in the batch. */
    public MultivectorBatch op(Multivector left) {
        return product(left, true, true);
    }

    /** {@code element ^ right} for every element in the batch. */
    public MultivectorBatch opRight(Multivector right) {
        return product(right, true, false);
    }

    /** {@code left ⌋ element} for every element in the batch. */
    public MultivectorBatch ip(Multivector left) {
        return product(left, false, true);
    }

    /** {@code element ⌋ right} for every element in the batch. */
    public MultivectorBatch ipRight(Multivector right) {
        return product(right, false, false);
    }

    /** Poincare dual of every element. */
    public MultivectorBatch dual() {
        int full = sig.fullMask();
        double[][] r = new double[sig.size][count];
        for (int m = 0; m < sig.size; m++) {
            if (!hasNonzero(blades[m])) continue;
            int sign = sig.dualSign(m);
            double[] src = blades[m];
            double[] dst = r[full ^ m];
            for (int i = 0; i < count; i++) dst[i] += sign * src[i];
        }
        return new MultivectorBatch(sig, count, r);
    }

    /** Inverse of {@link #dual()} for every element. */
    public MultivectorBatch undual() {
        int full = sig.fullMask();
        double[][] r = new double[sig.size][count];
        for (int m = 0; m < sig.size; m++) {
            int comp = full ^ m;
            if (!hasNonzero(blades[comp])) continue;
            int k = Integer.bitCount(m);
            int sign = sig.dualSign(comp) * (((k * (sig.n - k)) & 1) == 0 ? 1 : -1);
            double[] src = blades[comp];
            double[] dst = r[m];
            for (int i = 0; i < count; i++) dst[i] += sign * src[i];
        }
        return new MultivectorBatch(sig, count, r);
    }

    private MultivectorBatch product(Multivector fixed, boolean outer, boolean fixedLeft) {
        checkSignature(fixed);
        double[] f = fixed.coefficients();
        double[][] r = new double[sig.size][count];
        int size = sig.size;
        for (int a = 0; a < size; a++) {
            if (fixedLeft) {
                if (f[a] == 0) continue;
                for (int b = 0; b < size; b++) {
                    if (!hasNonzero(blades[b])) continue;
                    if (outer ? (a & b) != 0 : (a & b) != a) continue;
                    int sign = sig.bladeProductSign(a, b);
                    if (sign == 0) continue;
                    double factor = f[a] * sign;
                    addScaled(r[a ^ b], blades[b], factor);
                }
            } else {
                if (!hasNonzero(blades[a])) continue;
                for (int b = 0; b < size; b++) {
                    if (f[b] == 0) continue;
                    if (outer ? (a & b) != 0 : (a & b) != a) continue;
                    int sign = sig.bladeProductSign(a, b);
                    if (sign == 0) continue;
                    double factor = f[b] * sign;
                    addScaled(r[a ^ b], blades[a], factor);
                }
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    private void addScaled(double[] dst, double[] src, double factor) {
        for (int i = 0; i < count; i++) dst[i] += factor * src[i];
    }

    /** Zeroes every coefficient of every element. */
    public void clear() {
        for (double[] row : blades) java.util.Arrays.fill(row, 0);
    }

    private MultivectorBatch sandwichE3(Multivector versor, double[][] r) {
        double s = versor.get(0);
        double b1 = versor.get(3);
        double b2 = versor.get(5);
        double b3 = versor.get(6);
        double[] x = blades[1];
        double[] y = blades[2];
        double[] z = blades[4];
        for (int i = 0; i < count; i++) {
            double vx = x[i];
            double vy = y[i];
            double vz = z[i];
            r[1][i] = (s * s - b1 * b1 - b2 * b2 + b3 * b3) * vx
                    + 2 * (b1 * (s * vy + b3 * vz) + b2 * (s * vz - b3 * vy));
            r[2][i] = (s * s - b1 * b1 + b2 * b2 - b3 * b3) * vy
                    + 2 * (-vx * (b1 * s + b2 * b3) + vz * (b3 * s - b1 * b2));
            r[4][i] = (s * s + b1 * b1 - b2 * b2 - b3 * b3) * vz
                    + 2 * (vx * (b1 * b3 - b2 * s) - vy * (b1 * b2 + b3 * s));
        }
        return new MultivectorBatch(sig, count, r);
    }

    private MultivectorBatch sandwichPga(Multivector versor, double[][] r) {
        double s = versor.get(0);
        double b01 = versor.get(3);
        double b02 = versor.get(5);
        double b12 = versor.get(6);
        double b03 = versor.get(9);
        double b13 = versor.get(10);
        double b23 = versor.get(12);
        double ps = versor.get(15);
        double[] p7 = blades[7];
        double[] p11 = blades[11];
        double[] p13 = blades[13];
        double[] p14 = blades[14];
        for (int i = 0; i < count; i++) {
            double a7 = p7[i];
            double a11 = p11[i];
            double a13 = p13[i];
            double a14 = p14[i];
            r[7][i] = (s * s + b12 * b12 - b13 * b13 - b23 * b23) * a7
                    + 2 * (b12 * (b13 * a11 + b23 * a13) + s * (b23 * a11 - b13 * a13)
                            + a14 * (ps * b12 - b01 * b13 - b02 * b23 + b03 * s));
            r[11][i] = (s * s - b12 * b12 + b13 * b13 - b23 * b23) * a11
                    + 2 * (b12 * (b13 * a7 + s * a13) + b13 * b23 * a13 - s * b23 * a7
                            + a14 * (ps * b13 + b01 * b12 - b02 * s - b03 * b23));
            r[13][i] = (s * s - b12 * b12 - b13 * b13 + b23 * b23) * a13
                    + 2 * (b12 * (b23 * a7 - s * a11) + b13 * (b23 * a11 + s * a7)
                            + a14 * (ps * b23 + b01 * s + b02 * b12 + b03 * b13));
            r[14][i] = (s * s + b12 * b12 + b13 * b13 + b23 * b23) * a14;
        }
        return new MultivectorBatch(sig, count, r);
    }

    private boolean isE3() {
        return sig.n == 3 && sig.metric[0] == 1 && sig.metric[1] == 1 && sig.metric[2] == 1;
    }

    private boolean isPga3() {
        return sig.n == 4 && sig.metric[0] == 0 && sig.metric[1] == 1
                && sig.metric[2] == 1 && sig.metric[3] == 1;
    }

    private static boolean evenOnly(Multivector m) {
        for (int mask = 1; mask < m.sig.size; mask++) {
            if ((Integer.bitCount(mask) & 1) != 0 && m.get(mask) != 0) return false;
        }
        return true;
    }

    private boolean onlyBlades(int m0, int m1, int m2) {
        for (int mask = 0; mask < sig.size; mask++) {
            if (mask != m0 && mask != m1 && mask != m2 && !allZero(mask)) return false;
        }
        return true;
    }

    private boolean onlyBlades(int m0, int m1, int m2, int m3) {
        for (int mask = 0; mask < sig.size; mask++) {
            if (mask != m0 && mask != m1 && mask != m2 && mask != m3 && !allZero(mask)) return false;
        }
        return true;
    }

    private boolean allZero(int mask) {
        for (double v : blades[mask]) {
            if (v != 0) return false;
        }
        return true;
    }

    public MultivectorBatch add(MultivectorBatch o) {
        checkSame(o);
        double[][] r = new double[blades.length][count];
        for (int m = 0; m < blades.length; m++) {
            for (int i = 0; i < count; i++) {
                r[m][i] = blades[m][i] + o.blades[m][i];
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    public MultivectorBatch scale(double s) {
        double[][] r = new double[blades.length][count];
        for (int m = 0; m < blades.length; m++) {
            for (int i = 0; i < count; i++) {
                r[m][i] = blades[m][i] * s;
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    public MultivectorBatch reverse() {
        double[][] r = new double[blades.length][count];
        for (int m = 0; m < blades.length; m++) {
            double s = Signature.reverseSign(sig.grade(m));
            for (int i = 0; i < count; i++) {
                r[m][i] = blades[m][i] * s;
            }
        }
        return new MultivectorBatch(sig, count, r);
    }

    public Multivector[] toMultivectors() {
        Multivector[] out = new Multivector[count];
        for (int i = 0; i < count; i++) {
            out[i] = get(i);
        }
        return out;
    }

    private static boolean hasNonzero(double[] values) {
        for (double v : values) {
            if (v != 0) return true;
        }
        return false;
    }

    private void checkIndex(int i) {
        if (i < 0 || i >= count) {
            throw new IndexOutOfBoundsException("index " + i + " outside 0.." + (count - 1));
        }
    }

    private void checkSignature(Multivector value) {
        if (value.sig != sig) {
            throw new IllegalArgumentException("different algebras");
        }
    }

    private void checkSame(MultivectorBatch o) {
        if (o.sig != sig || o.count != count) {
            throw new IllegalArgumentException("batch shape mismatch");
        }
    }
}
