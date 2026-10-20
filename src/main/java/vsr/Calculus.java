package vsr;

/**
 * Numerical geometric calculus on R^3 (Cl(3,0)).
 *
 * <p>The vector derivative of a multivector field is
 * <pre>  nabla F = e1 dF/dx + e2 dF/dy + e3 dF/dz</pre>
 * Its scalar part is the divergence and its bivector part is (minus the dual
 * of) the curl, so all three classical operators come from one object.
 */
public final class Calculus {
    /** Central-difference step. */
    public static final double STEP = 1e-4;

    private static final Multivector E1 = Multivector.blade(E3.SIG, 0b001, 1);
    private static final Multivector E2 = Multivector.blade(E3.SIG, 0b010, 1);
    private static final Multivector E3B = Multivector.blade(E3.SIG, 0b100, 1);
    private static final Multivector I = Multivector.blade(E3.SIG, 0b111, 1);

    private Calculus() {}

    /** Partial derivative along axis 0, 1 or 2. */
    public static Multivector partial(Field f, double x, double y, double z, int axis) {
        return partial(f, x, y, z, axis, STEP);
    }

    public static Multivector partial(Field f, double x, double y, double z, int axis, double h) {
        if (axis < 0 || axis > 2) throw new IllegalArgumentException("axis must be 0, 1 or 2");
        double[] a = {x, y, z};
        double[] p = a.clone();
        p[axis] += h;
        Multivector plus = f.at(p);
        p[axis] = a[axis] - h;
        Multivector minus = f.at(p);
        return plus.sub(minus).scale(0.5 / h);
    }

    /** Directional derivative (v . nabla) F; v is normalized internally. */
    public static Multivector directional(Field f, double x, double y, double z,
                                          double vx, double vy, double vz) {
        double len = Math.sqrt(vx * vx + vy * vy + vz * vz);
        if (len == 0) throw new IllegalArgumentException("zero direction");
        double ux = vx / len, uy = vy / len, uz = vz / len;
        Multivector plus = f.at(x + STEP * ux, y + STEP * uy, z + STEP * uz);
        Multivector minus = f.at(x - STEP * ux, y - STEP * uy, z - STEP * uz);
        return plus.sub(minus).scale(0.5 / STEP);
    }

    /** Vector derivative nabla F. */
    public static Multivector vectorDerivative(Field f, double x, double y, double z) {
        Multivector dx = partial(f, x, y, z, 0);
        Multivector dy = partial(f, x, y, z, 1);
        Multivector dz = partial(f, x, y, z, 2);
        return E1.gp(dx).add(E2.gp(dy)).add(E3B.gp(dz));
    }

    /** Gradient of a scalar field (grade-1 part of nabla F). */
    public static Multivector gradient(Field f, double x, double y, double z) {
        return vectorDerivative(f, x, y, z).grade(1);
    }

    /** Divergence of a vector field (scalar part of nabla F). */
    public static double divergence(Field f, double x, double y, double z) {
        return vectorDerivative(f, x, y, z).scalarPart();
    }

    /** Curl of a vector field: the vector dual of the bivector part, -I (nabla ^ F). */
    public static Multivector curl(Field f, double x, double y, double z) {
        Multivector bivector = vectorDerivative(f, x, y, z).grade(2);
        return I.gp(bivector).neg();
    }
}
