package vsr;

import java.util.Arrays;
import java.util.TreeMap;

/**
 * Sparse multivector (COO: sorted masks + values) for high-dimensional and
 * low-density algebras. Products cost O(k1*k2) with an open-addressing
 * accumulator, or O(k1*k2 + 2^n) for dense-ish inputs, instead of O(4^n).
 */
public final class SparseMultivector {
    public final Signature sig;
    private final int[] masks;
    private final double[] values;

    private SparseMultivector(Signature sig, int[] masks, double[] values) {
        this.sig = sig;
        this.masks = masks;
        this.values = values;
    }

    public static SparseMultivector zero(Signature sig) {
        return new SparseMultivector(sig, new int[0], new double[0]);
    }

    public static SparseMultivector fromDense(Multivector m) {
        int nz = 0;
        for (int mask = 0; mask < m.sig.size; mask++) {
            if (m.get(mask) != 0) nz++;
        }
        int[] masks = new int[nz];
        double[] values = new double[nz];
        int p = 0;
        for (int mask = 0; mask < m.sig.size; mask++) {
            double v = m.get(mask);
            if (v != 0) {
                masks[p] = mask;
                values[p] = v;
                p++;
            }
        }
        return new SparseMultivector(m.sig, masks, values);
    }

    /** Builds from arbitrary mask/value pairs; duplicates are summed. */
    public static SparseMultivector of(Signature sig, int[] masks, double[] values) {
        if (masks.length != values.length) {
            throw new IllegalArgumentException("mask/value length mismatch");
        }
        TreeMap<Integer, Double> merged = new TreeMap<>();
        for (int i = 0; i < masks.length; i++) {
            if (masks[i] < 0 || masks[i] >= sig.size) {
                throw new IllegalArgumentException("mask outside signature: " + masks[i]);
            }
            if (values[i] == 0) continue;
            merged.merge(masks[i], values[i], Double::sum);
        }
        int[] m = new int[merged.size()];
        double[] v = new double[merged.size()];
        int p = 0;
        for (var e : merged.entrySet()) {
            if (e.getValue() == 0) continue;
            m[p] = e.getKey();
            v[p] = e.getValue();
            p++;
        }
        return new SparseMultivector(sig, Arrays.copyOf(m, p), Arrays.copyOf(v, p));
    }

    public int nnz() {
        return masks.length;
    }

    public int mask(int i) {
        return masks[i];
    }

    public double value(int i) {
        return values[i];
    }

    public Multivector toDense() {
        double[] c = new double[sig.size];
        for (int i = 0; i < masks.length; i++) c[masks[i]] = values[i];
        return new Multivector(sig, c);
    }

    /**
     * Geometric product. Low pair counts use an open-addressing accumulator
     * (no {@code 2^n} scan); dense-ish inputs switch to a dense accumulator.
     */
    public SparseMultivector gp(SparseMultivector o) {
        check(o);
        long pairs = (long) masks.length * o.masks.length;
        if (pairs * 4 < sig.size) {
            return gpHashed(o);
        }
        return gpDense(o);
    }

    private SparseMultivector gpHashed(SparseMultivector o) {
        int capacity = 16;
        while (capacity < Math.max(16, masks.length * o.masks.length * 2)) capacity <<= 1;
        int[] keys = new int[capacity];
        Arrays.fill(keys, -1);
        double[] acc = new double[capacity];
        int maskBits = capacity - 1;
        int[][] table = sig.signTable();
        for (int i = 0; i < masks.length; i++) {
            int ma = masks[i];
            double va = values[i];
            int[] row = table != null ? table[ma] : null;
            for (int j = 0; j < o.masks.length; j++) {
                int mb = o.masks[j];
                int sign = row != null ? row[mb] : sig.computeBladeProductSign(ma, mb);
                if (sign == 0) continue;
                int dest = ma ^ mb;
                int slot = mix(dest) & maskBits;
                while (keys[slot] != -1 && keys[slot] != dest) slot = (slot + 1) & maskBits;
                keys[slot] = dest;
                acc[slot] += sign * va * o.values[j];
            }
        }
        int count = 0;
        for (int slot = 0; slot < capacity; slot++) {
            if (keys[slot] != -1 && acc[slot] != 0) count++;
        }
        long[] order = new long[count];
        int p = 0;
        for (int slot = 0; slot < capacity; slot++) {
            if (keys[slot] != -1 && acc[slot] != 0) {
                order[p++] = ((long) keys[slot] << 32) | slot;
            }
        }
        Arrays.sort(order);
        int[] m = new int[count];
        double[] v = new double[count];
        for (int i = 0; i < count; i++) {
            m[i] = (int) (order[i] >>> 32);
            v[i] = acc[(int) order[i]];
        }
        return new SparseMultivector(sig, m, v);
    }

    private int mix(int x) {
        x ^= x >>> 16;
        x *= 0x7feb352d;
        x ^= x >>> 15;
        x *= 0x846ca68b;
        return x ^ (x >>> 16);
    }

    /** Geometric product with a dense accumulator over the output support. */
    private SparseMultivector gpDense(SparseMultivector o) {
        int size = sig.size;
        double[] acc = new double[size];
        boolean[] seen = new boolean[size];
        int capacity = (int) Math.min((long) size, (long) masks.length * o.masks.length);
        int[] touched = new int[Math.max(0, capacity)];
        int count = 0;
        int[][] table = sig.signTable();
        for (int i = 0; i < masks.length; i++) {
            int ma = masks[i];
            double va = values[i];
            int[] row = table != null ? table[ma] : null;
            for (int j = 0; j < o.masks.length; j++) {
                int mb = o.masks[j];
                int sign = row != null ? row[mb] : sig.computeBladeProductSign(ma, mb);
                if (sign == 0) continue;
                int dest = ma ^ mb;
                if (!seen[dest]) {
                    seen[dest] = true;
                    touched[count++] = dest;
                }
                acc[dest] += sign * va * o.values[j];
            }
        }
        int[] sorted = Arrays.copyOf(touched, count);
        Arrays.sort(sorted);
        int nz = 0;
        for (int mask : sorted) {
            if (acc[mask] != 0) nz++;
        }
        int[] m = new int[nz];
        double[] v = new double[nz];
        int p = 0;
        for (int mask : sorted) {
            if (acc[mask] != 0) {
                m[p] = mask;
                v[p] = acc[mask];
                p++;
            }
        }
        return new SparseMultivector(sig, m, v);
    }

    /** {@code this * x * ~this}. */
    public SparseMultivector sandwich(SparseMultivector x) {
        return gp(x).gp(reverse());
    }

    public SparseMultivector reverse() {
        double[] v = new double[values.length];
        for (int i = 0; i < masks.length; i++) {
            int k = Integer.bitCount(masks[i]);
            v[i] = (((k * (k - 1) / 2) & 1) == 0 ? 1 : -1) * values[i];
        }
        return new SparseMultivector(sig, masks.clone(), v);
    }

    public SparseMultivector add(SparseMultivector o) {
        check(o);
        TreeMap<Integer, Double> merged = new TreeMap<>();
        for (int i = 0; i < masks.length; i++) merged.merge(masks[i], values[i], Double::sum);
        for (int i = 0; i < o.masks.length; i++) merged.merge(o.masks[i], o.values[i], Double::sum);
        int[] m = new int[merged.size()];
        double[] v = new double[merged.size()];
        int p = 0;
        for (var e : merged.entrySet()) {
            if (e.getValue() == 0) continue;
            m[p] = e.getKey();
            v[p] = e.getValue();
            p++;
        }
        return new SparseMultivector(sig, Arrays.copyOf(m, p), Arrays.copyOf(v, p));
    }

    public SparseMultivector scale(double s) {
        double[] v = new double[values.length];
        for (int i = 0; i < values.length; i++) v[i] = values[i] * s;
        return new SparseMultivector(sig, masks.clone(), v);
    }

    public double norm2() {
        double sum = 0;
        for (int i = 0; i < masks.length; i++) {
            int mask = masks[i];
            int metric = bladesSquareSign(mask);
            sum += metric * values[i] * values[i];
        }
        return sum;
    }

    /** Product of metric entries for the generators in the blade. */
    private int bladesSquareSign(int mask) {
        int sign = 1;
        for (int j = 0; j < sig.n; j++) {
            if ((mask & (1 << j)) == 0) continue;
            int m = sig.metric[j];
            if (m == 0) return 0;
            sign *= m;
        }
        return sign;
    }

    public double maxAbs() {
        double max = 0;
        for (double v : values) max = Math.max(max, Math.abs(v));
        return max;
    }

    public double distanceTo(SparseMultivector o) {
        check(o);
        return toDense().distanceTo(o.toDense());
    }

    private void check(SparseMultivector o) {
        if (o.sig != sig) throw new IllegalArgumentException("different algebras");
    }
}
