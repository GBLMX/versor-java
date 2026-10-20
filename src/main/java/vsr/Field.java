package vsr;

/** A multivector-valued field on R^3: F(x, y, z). */
@FunctionalInterface
public interface Field {
    Multivector at(double x, double y, double z);

    default Multivector at(double[] p) {
        return at(p[0], p[1], p[2]);
    }
}
