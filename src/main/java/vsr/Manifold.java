package vsr;

/**
 * An embedded submanifold of R^3 given by a parameterization.
 *
 * <p>A curve has dimension 1 and a surface dimension 2. Derivatives of the
 * parameterization give the tangent blade; in 3D the dual of a surface's
 * tangent bivector is its normal vector.
 */
public final class Manifold {
    public interface Param {
        double[] at(double... u);
    }

    private static final double H = 1e-5;
    private static final Multivector I = Multivector.blade(E3.SIG, 0b111, 1);

    private final int dimension;
    private final Param param;

    private Manifold(int dimension, Param param) {
        this.dimension = dimension;
        this.param = param;
    }

    public static Manifold curve(Param param) {
        return new Manifold(1, param);
    }

    public static Manifold surface(Param param) {
        return new Manifold(2, param);
    }

    public int dimension() {
        return dimension;
    }

    public double[] point(double... u) {
        return param.at(u);
    }

    /** Partial derivative of the parameterization, as a vector. */
    public Multivector partial(double[] u, int axis) {
        if (axis < 0 || axis >= dimension) throw new IllegalArgumentException("bad axis");
        double[] plus = u.clone();
        double[] minus = u.clone();
        plus[axis] += H;
        minus[axis] -= H;
        double[] a = param.at(plus);
        double[] b = param.at(minus);
        return E3.vec((a[0] - b[0]) / (2 * H), (a[1] - b[1]) / (2 * H), (a[2] - b[2]) / (2 * H));
    }

    /** Tangent blade: a vector for a curve, a bivector for a surface. */
    public Multivector tangentBlade(double... u) {
        Multivector t = partial(u, 0);
        if (dimension == 1) return t;
        return t.op(partial(u, 1));
    }

    /** Jacobian: speed for a curve, area element for a surface. */
    public double jacobian(double... u) {
        return dimension == 1 ? partial(u, 0).norm() : tangentBlade(u).norm();
    }

    /** Unit normal of a surface (dual of the tangent bivector). */
    public Multivector normal(double... u) {
        if (dimension != 2) throw new IllegalStateException("normal is defined for surfaces");
        Multivector blade = I.gp(tangentBlade(u)).neg().grade(1);
        return blade.unit();
    }

    /** Unit tangent of a curve. */
    public Multivector tangent(double... u) {
        if (dimension != 1) throw new IllegalStateException("tangent is defined for curves");
        return partial(u, 0).unit();
    }
}
