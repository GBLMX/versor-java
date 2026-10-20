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

    /** Centre (x, y, z) and radius of an OPNS sphere. */
    public static double[] sphereCoords(Multivector sphere) {
        check(sphere);
        double w = sphere.get(E5) - sphere.get(E4);
        if (Math.abs(w) < 1e-12) {
            throw new IllegalArgumentException("sphere at infinity");
        }
        double[] c = {sphere.get(E1) / w, sphere.get(E2) / w, sphere.get(E3) / w};
        double r2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2]
                - (sphere.get(E4) + sphere.get(E5)) / w;
        if (r2 < -1e-9) {
            throw new IllegalArgumentException("sphere has negative radius squared");
        }
        return new double[]{c[0], c[1], c[2], Math.sqrt(Math.max(0.0, r2))};
    }

    /** Plane coefficients {@code (nx, ny, nz, d)} of an encoded plane. */
    public static double[] planeCoords(Multivector plane) {
        check(plane);
        return new double[]{
            plane.get(E1), plane.get(E2), plane.get(E3),
            -0.5 * (plane.get(E4) + plane.get(E5))
        };
    }

    /** Signed distance from a point to a plane; positive on the normal side. */
    public static double signedDistancePointPlane(Multivector plane, Multivector query) {
        double[] pi = planeCoords(plane);
        double[] p = pointCoords(query);
        double n2 = pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2];
        if (n2 < 1e-24) {
            throw new IllegalArgumentException("degenerate plane");
        }
        return (pi[0] * p[0] + pi[1] * p[1] + pi[2] * p[2] + pi[3]) / Math.sqrt(n2);
    }

    /** Euclidean distance from a point to a plane. */
    public static double distancePointPlane(Multivector plane, Multivector query) {
        return Math.abs(signedDistancePointPlane(plane, query));
    }

    /** Orthogonal projection of a point onto a plane. */
    public static Multivector closestPointOnPlane(Multivector plane, Multivector query) {
        double[] pi = planeCoords(plane);
        double[] p = pointCoords(query);
        double n2 = pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2];
        if (n2 < 1e-24) {
            throw new IllegalArgumentException("degenerate plane");
        }
        double t = (pi[0] * p[0] + pi[1] * p[1] + pi[2] * p[2] + pi[3]) / n2;
        return point(p[0] - t * pi[0], p[1] - t * pi[1], p[2] - t * pi[2]);
    }

    /** Euclidean distance from a point to the surface of a sphere. */
    public static double distancePointSphere(Multivector sphere, Multivector query) {
        double[] s = sphereCoords(sphere);
        double[] p = pointCoords(query);
        double dx = p[0] - s[0], dy = p[1] - s[1], dz = p[2] - s[2];
        return Math.abs(Math.sqrt(dx * dx + dy * dy + dz * dz) - s[3]);
    }

    /** Closest point on a sphere's surface; a sphere centre maps to (c+r, c_y, c_z). */
    public static Multivector closestPointOnSphere(Multivector sphere, Multivector query) {
        double[] s = sphereCoords(sphere);
        double[] p = pointCoords(query);
        double dx = p[0] - s[0], dy = p[1] - s[1], dz = p[2] - s[2];
        double length = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (length < 1e-12) {
            return point(s[0] + s[3], s[1], s[2]);
        }
        double scale = s[3] / length;
        return point(s[0] + dx * scale, s[1] + dy * scale, s[2] + dz * scale);
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

    private static void check(Multivector m) {
        if (m.sig != SIG) {
            throw new IllegalArgumentException("expected conformal Cl(4,1) multivector");
        }
    }
}
