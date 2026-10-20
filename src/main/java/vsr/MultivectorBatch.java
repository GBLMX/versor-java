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

    /** {@code versor * element * ~versor} for every element in the batch. */
    public MultivectorBatch sandwich(Multivector versor) {
        checkSignature(versor);
        return gp(versor).gpRight(versor.reverse());
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
            int k = sig.grade(m);
            double s = ((k * (k - 1) / 2) & 1) == 0 ? 1 : -1;
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
