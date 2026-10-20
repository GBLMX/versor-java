package vsr;

/**
 * Conformal model of 3D Euclidean space in {@link Signature#cga3()}:
 * Cl(4,1) with e1..e4 squaring to +1 and e5 squaring to -1.
 *
 * <p>Embeddings follow the standard null-basis convention
 * <pre>
 *   n0 = (e5 - e4)/2,  ninf = e4 + e5,  n0.ninf = -1
 *   P(x) = x + |x|^2/2 ninf + n0
 * </pre>
 * so {@code P.P = 0} and {@code P.Q = -|p - q|^2 / 2}. A sphere is
 * {@code S = P(c) - r^2/2 ninf}, a plane is {@code pi = n - d ninf} for
 * {@code n.x + d = 0}, and a circle is the wedge of its sphere and plane.
 *
 * <p>Plane constants match {@link P3#plane}: {@code n.x + d = 0}.
 */
public final class Conformal {
    public static final Signature SIG = Signature.cga3();

    private static final int E1 = 1;
    private static final int E2 = 2;
    private static final int E3 = 4;
    private static final int E4 = 8;
    private static final int E5 = 16;
    private static final Multivector NINF =
            Multivector.blade(SIG, E4, 1.0).add(Multivector.blade(SIG, E5, 1.0));
    private static final Multivector N0 =
            Multivector.blade(SIG, E5, 0.5).sub(Multivector.blade(SIG, E4, 0.5));

    private Conformal() {}

    public static Multivector point(double x, double y, double z) {
        double half = 0.5 * (x * x + y * y + z * z);
        return Multivector.vector(SIG, x, y, z, 0.0, 0.0)
                .add(NINF.scale(half))
                .add(N0);
    }

    /** Inverse of {@link #point(double, double, double)}. */
    public static double[] pointCoords(Multivector p) {
        double w = p.get(E5) - p.get(E4);
        if (Math.abs(w) < 1e-12) {
            throw new IllegalArgumentException("point at infinity");
        }
        return new double[]{p.get(E1) / w, p.get(E2) / w, p.get(E3) / w};
    }

    /** Sphere with the given centre and radius. */
    public static Multivector sphere(double cx, double cy, double cz, double radius) {
        return point(cx, cy, cz).sub(NINF.scale(0.5 * radius * radius));
    }

    /** Plane {@code nx*x + ny*y + nz*z + d = 0}, same convention as P3. */
    public static Multivector plane(double nx, double ny, double nz, double d) {
        return Multivector.vector(SIG, nx, ny, nz, 0.0, 0.0)
                .add(NINF.scale(-d));
    }

    /**
     * Circle around {@code c} in the plane with normal {@code n}, as the OPNS
     * trivector through three points of the circle. A point {@code P} lies on
     * it when {@code P ^ circle == 0}.
     */
    public static Multivector circle(double cx, double cy, double cz,
            double nx, double ny, double nz, double radius) {
        double length = Math.sqrt(nx * nx + ny * ny + nz * nz);
        if (length < 1e-12) {
            throw new IllegalArgumentException("zero normal");
        }
        double ux = nx / length, uy = ny / length, uz = nz / length;
        // A helper axis that is never parallel to the normal.
        double hx = Math.abs(ux) < 0.9 ? 1.0 : 0.0;
        double hy = Math.abs(ux) < 0.9 ? 0.0 : 1.0;
        double hz = 0.0;
        double ax = hy * uz - hz * uy;
        double ay = hz * ux - hx * uz;
        double az = hx * uy - hy * ux;
        double aLength = Math.sqrt(ax * ax + ay * ay + az * az);
        ax /= aLength;
        ay /= aLength;
        az /= aLength;
        double bx = uy * az - uz * ay;
        double by = uz * ax - ux * az;
        double bz = ux * ay - uy * ax;

        Multivector p1 = point(cx + radius * ax, cy + radius * ay, cz + radius * az);
        Multivector p2 = point(cx + radius * bx, cy + radius * by, cz + radius * bz);
        Multivector p3 = point(cx - radius * ax, cy - radius * ay, cz - radius * az);
        return p1.op(p2).op(p3);
    }

    /** Inner product of two conformal vectors. */
    public static double dot(Multivector a, Multivector b) {
        return a.gp(b).scalarPart();
    }

    /** Squared Euclidean distance between two (possibly non-normalised) points. */
    public static double distanceSq(Multivector a, Multivector b) {
        return -2.0 * dot(a, b);
    }
}
