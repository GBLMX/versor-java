package vsr;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * Optional SIMD versions of the E3 and PGA3 batch sandwich kernels.
 *
 * <p>Requires the incubator module:
 * {@code javac --release 25 --add-modules jdk.incubator.vector}.
 * This class is a prototype and is not referenced by the normal build; the
 * scalar {@link MultivectorBatch} paths remain the production implementation.
 */
public final class SimdSandwich {
    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;
    private static final DoubleVector TWO = DoubleVector.broadcast(SPECIES, 2.0);

    private SimdSandwich() {}

    /** E3 even-versor sandwich over three contiguous coordinate arrays. */
    public static void e3(double s, double b1, double b2, double b3,
                          double[] x, double[] y, double[] z,
                          double[] rx, double[] ry, double[] rz) {
        checkLength(x, y, z, rx, ry, rz);
        DoubleVector vs = DoubleVector.broadcast(SPECIES, s);
        DoubleVector v1 = DoubleVector.broadcast(SPECIES, b1);
        DoubleVector v2 = DoubleVector.broadcast(SPECIES, b2);
        DoubleVector v3 = DoubleVector.broadcast(SPECIES, b3);
        DoubleVector s2 = vs.mul(vs);
        DoubleVector b1s = v1.mul(v1);
        DoubleVector b2s = v2.mul(v2);
        DoubleVector b3s = v3.mul(v3);
        DoubleVector x11 = s2.sub(b1s).sub(b2s).add(b3s);
        DoubleVector x22 = s2.sub(b1s).add(b2s).sub(b3s);
        DoubleVector x33 = s2.add(b1s).sub(b2s).sub(b3s);
        int i = 0;
        int bound = SPECIES.loopBound(x.length);
        for (; i < bound; i += SPECIES.length()) {
            DoubleVector vx = DoubleVector.fromArray(SPECIES, x, i);
            DoubleVector vy = DoubleVector.fromArray(SPECIES, y, i);
            DoubleVector vz = DoubleVector.fromArray(SPECIES, z, i);
            DoubleVector outX = x11.mul(vx).add(TWO.mul(
                    v1.mul(vs.mul(vy).add(v3.mul(vz)))
                            .add(v2.mul(vs.mul(vz).sub(v3.mul(vy))))));
            DoubleVector outY = x22.mul(vy).add(TWO.mul(
                    vx.neg().mul(v1.mul(vs).add(v2.mul(v3)))
                            .add(vz.mul(v3.mul(vs).sub(v1.mul(v2))))));
            DoubleVector outZ = x33.mul(vz).add(TWO.mul(
                    vx.mul(v1.mul(v3).sub(v2.mul(vs)))
                            .sub(vy.mul(v1.mul(v2).add(v3.mul(vs))))));
            outX.intoArray(rx, i);
            outY.intoArray(ry, i);
            outZ.intoArray(rz, i);
        }
        for (; i < x.length; i++) {
            double vx = x[i];
            double vy = y[i];
            double vz = z[i];
            rx[i] = (s * s - b1 * b1 - b2 * b2 + b3 * b3) * vx
                    + 2 * (b1 * (s * vy + b3 * vz) + b2 * (s * vz - b3 * vy));
            ry[i] = (s * s - b1 * b1 + b2 * b2 - b3 * b3) * vy
                    + 2 * (-vx * (b1 * s + b2 * b3) + vz * (b3 * s - b1 * b2));
            rz[i] = (s * s + b1 * b1 - b2 * b2 - b3 * b3) * vz
                    + 2 * (vx * (b1 * b3 - b2 * s) - vy * (b1 * b2 + b3 * s));
        }
    }

    /** PGA3 point sandwich over the four trivector coefficient arrays. */
    public static void pga(double s, double b01, double b02, double b12,
                           double b03, double b13, double b23, double ps,
                           double[] p7, double[] p11, double[] p13, double[] p14,
                           double[] r7, double[] r11, double[] r13, double[] r14) {
        checkLength(p7, p11, p13, p14);
        checkLength(r7, r11, r13, r14);
        DoubleVector vs = DoubleVector.broadcast(SPECIES, s);
        DoubleVector v01 = DoubleVector.broadcast(SPECIES, b01);
        DoubleVector v02 = DoubleVector.broadcast(SPECIES, b02);
        DoubleVector v12 = DoubleVector.broadcast(SPECIES, b12);
        DoubleVector v03 = DoubleVector.broadcast(SPECIES, b03);
        DoubleVector v13 = DoubleVector.broadcast(SPECIES, b13);
        DoubleVector v23 = DoubleVector.broadcast(SPECIES, b23);
        DoubleVector vps = DoubleVector.broadcast(SPECIES, ps);
        DoubleVector s2 = vs.mul(vs);
        DoubleVector b12s = v12.mul(v12);
        DoubleVector b13s = v13.mul(v13);
        DoubleVector b23s = v23.mul(v23);
        DoubleVector k7 = s2.add(b12s).sub(b13s).sub(b23s);
        DoubleVector k11 = s2.sub(b12s).add(b13s).sub(b23s);
        DoubleVector k13 = s2.sub(b12s).sub(b13s).add(b23s);
        DoubleVector k14 = s2.add(b12s).add(b13s).add(b23s);
        DoubleVector t7ps = vps.mul(v12).sub(v01.mul(v13)).sub(v02.mul(v23)).add(v03.mul(vs));
        DoubleVector t11ps = vps.mul(v13).add(v01.mul(v12)).sub(v02.mul(vs)).sub(v03.mul(v23));
        DoubleVector t13ps = vps.mul(v23).add(v01.mul(vs)).add(v02.mul(v12)).add(v03.mul(v13));
        int i = 0;
        int bound = SPECIES.loopBound(p7.length);
        for (; i < bound; i += SPECIES.length()) {
            DoubleVector a7 = DoubleVector.fromArray(SPECIES, p7, i);
            DoubleVector a11 = DoubleVector.fromArray(SPECIES, p11, i);
            DoubleVector a13 = DoubleVector.fromArray(SPECIES, p13, i);
            DoubleVector a14 = DoubleVector.fromArray(SPECIES, p14, i);
            DoubleVector o7 = k7.mul(a7).add(TWO.mul(
                    v12.mul(v13.mul(a11).add(v23.mul(a13)))
                            .add(vs.mul(v23.mul(a11).sub(v13.mul(a13))))
                            .add(a14.mul(t7ps))));
            DoubleVector o11 = k11.mul(a11).add(TWO.mul(
                    v12.mul(v13.mul(a7).add(vs.mul(a13)))
                            .add(v13.mul(v23).mul(a13))
                            .sub(vs.mul(v23).mul(a7))
                            .add(a14.mul(t11ps))));
            DoubleVector o13 = k13.mul(a13).add(TWO.mul(
                    v12.mul(v23.mul(a7).sub(vs.mul(a11)))
                            .add(v13.mul(v23.mul(a11).add(vs.mul(a7))))
                            .add(a14.mul(t13ps))));
            o7.intoArray(r7, i);
            o11.intoArray(r11, i);
            o13.intoArray(r13, i);
            k14.mul(a14).intoArray(r14, i);
        }
        for (; i < p7.length; i++) {
            double a7 = p7[i];
            double a11 = p11[i];
            double a13 = p13[i];
            double a14 = p14[i];
            r7[i] = (s * s + b12 * b12 - b13 * b13 - b23 * b23) * a7
                    + 2 * (b12 * (b13 * a11 + b23 * a13) + s * (b23 * a11 - b13 * a13)
                    + a14 * (ps * b12 - b01 * b13 - b02 * b23 + b03 * s));
            r11[i] = (s * s - b12 * b12 + b13 * b13 - b23 * b23) * a11
                    + 2 * (b12 * (b13 * a7 + s * a13) + b13 * b23 * a13 - s * b23 * a7
                    + a14 * (ps * b13 + b01 * b12 - b02 * s - b03 * b23));
            r13[i] = (s * s - b12 * b12 - b13 * b13 + b23 * b23) * a13
                    + 2 * (b12 * (b23 * a7 - s * a11) + b13 * (b23 * a11 + s * a7)
                    + a14 * (ps * b23 + b01 * s + b02 * b12 + b03 * b13));
            r14[i] = (s * s + b12 * b12 + b13 * b13 + b23 * b23) * a14;
        }
    }

    private static void checkLength(double[]... arrays) {
        int length = arrays[0].length;
        for (double[] array : arrays) {
            if (array.length != length) throw new IllegalArgumentException("array length mismatch");
        }
    }
}
