package vsr;

import java.util.Arrays;

/**
 * A Clifford algebra Cl(p,q,r): generators square to +1, -1 or 0.
 *
 * <p>Generators are indexed 0..n-1. A basis blade is encoded as a bit mask
 * (bit i set means generator i is present). The full set of 2^n blades is
 * ordered by mask value.
 */
public final class Signature {
    public final int n;
    public final int[] metric;
    public final int size;
    private final int full;
    private final int nameOffset;

    private Signature(int nameOffset, int[] metric) {
        if (metric.length > 10) {
            throw new IllegalArgumentException("Clifford algebras above 10 generators are not supported");
        }
        this.n = metric.length;
        this.nameOffset = nameOffset;
        this.metric = metric.clone();
        this.size = 1 << n;
        this.full = size - 1;
        for (int m : metric) {
            if (m < -1 || m > 1) {
                throw new IllegalArgumentException("metric entries must be -1, 0 or 1");
            }
        }
    }

    /** Cl(metric...) with generator names e0, e1, ... */
    public static Signature of(int... metric) {
        return new Signature(0, metric);
    }

    public static Signature euclidean(int n) {
        int[] m = new int[n];
        Arrays.fill(m, 1);
        return new Signature(1, m);
    }

    /** 3D projective GA: e1^2=e2^2=e3^2=1, e0^2=0 (e0 first). */
    public static Signature pga3() {
        return new Signature(0, new int[]{0, 1, 1, 1});
    }

    /** 3D Euclidean GA (quaternions). */
    public static Signature e3() {
        return new Signature(1, new int[]{1, 1, 1});
    }

    /** Spacetime algebra Cl(1,3); e0 is timelike. */
    public static Signature sta() {
        return new Signature(0, new int[]{1, -1, -1, -1});
    }

    /** 3D conformal model: e1..e4 spacelike, e5 timelike. */
    public static Signature cga3() {
        return new Signature(1, new int[]{1, 1, 1, 1, -1});
    }

    public int fullMask() {
        return full;
    }

    public int grade(int mask) {
        return Integer.bitCount(mask);
    }

    public boolean isScalar(int mask) {
        return mask == 0;
    }

    /** Name of a basis blade, e.g. "e0e1". */
    public String bladeName(int mask) {
        if (mask == 0) return "1";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            if ((mask & (1 << i)) != 0) sb.append('e').append(i + nameOffset);
        }
        return sb.toString();
    }

    /**
     * Sign and metric factor of the geometric product of two basis blades.
     * Returns 0 when a null generator would be squared.
     */
    public int bladeProductSign(int a, int b) {
        int sign = 1;
        for (int j = 0; j < n; j++) {
            if ((b & (1 << j)) == 0) continue;
            // move e_j left through the factors of a that have a higher index
            if ((Integer.bitCount(a >>> (j + 1)) & 1) != 0) sign = -sign;
            if ((a & (1 << j)) != 0) {
                int m = metric[j];
                if (m == 0) return 0;
                sign *= m;
            }
        }
        return sign;
    }

    @Override
    public String toString() {
        return "Cl(" + Arrays.toString(metric) + ")";
    }
}
