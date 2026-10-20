package vsr;

/**
 * 3D projective GA for graphics: points are trivectors, planes are vectors,
 * lines are bivectors, rigid motions are motors.
 *
 * <p>Basis order is e0 (null), e1, e2, e3.
 * Point: P(x,y,z) = e123 - z e012 + y e013 - x e023.
 * Plane: pi = nx e1 + ny e2 + nz e3 + d e0, with nx*x + ny*y + nz*z + d = 0.
 * Translator: T = 1 - (tx e01 + ty e02 + tz e03)/2 (T p ~T moves p by +t).
 */
public final class P3 {
    public static final Signature SIG = Signature.pga3();
    public static final Multivector I = Multivector.blade(SIG, 0b1111, 1);

    private P3() {}

    // blade masks
    private static final int E01 = 0b0011;
    private static final int E02 = 0b0101;
    private static final int E03 = 0b1001;
    private static final int E12 = 0b0110;
    private static final int E13 = 0b1010;
    private static final int E23 = 0b1100;
    private static final int E123 = 0b1110;
    private static final int E012 = 0b0111;
    private static final int E013 = 0b1011;
    private static final int E023 = 0b1101;

    public static Multivector point(double x, double y, double z) {
        return Multivector.blade(SIG, E123, 1)
                .add(Multivector.blade(SIG, E012, -z))
                .add(Multivector.blade(SIG, E013, y))
                .add(Multivector.blade(SIG, E023, -x));
    }

    /** Euclidean coordinates of a point trivector. */
    public static double[] pointCoords(Multivector p) {
        double w = p.get(E123);
        if (Math.abs(w) < 1e-12) throw new IllegalArgumentException("point at infinity");
        return new double[]{-p.get(E023) / w, p.get(E013) / w, -p.get(E012) / w};
    }

    /** Plane nx*x + ny*y + nz*z + d = 0. */
    public static Multivector plane(double nx, double ny, double nz, double d) {
        return Multivector.vector(SIG, d, nx, ny, nz);
    }

    public static double[] planeCoords(Multivector plane) {
        return new double[]{
                plane.get(1 << 1), plane.get(1 << 2), plane.get(1 << 3), plane.get(1)
        };
    }

    /** T = 1 - (t ^ e0)/2 = 1 - (tx e01 + ty e02 + tz e03)/2. */
    public static Multivector translator(double tx, double ty, double tz) {
        return Multivector.scalar(SIG, 1)
                .sub(Multivector.blade(SIG, E01, 0.5 * tx))
                .sub(Multivector.blade(SIG, E02, 0.5 * ty))
                .sub(Multivector.blade(SIG, E03, 0.5 * tz));
    }

    /** Rotation about the unit axis (ax,ay,az), same convention as {@link E3}. */
    public static Multivector rotor(double ax, double ay, double az, double angle) {
        double len = Math.sqrt(ax * ax + ay * ay + az * az);
        if (len == 0) throw new IllegalArgumentException("zero axis");
        Multivector b = Multivector.blade(SIG, E23, ax / len)
                .add(Multivector.blade(SIG, E13, -ay / len))
                .add(Multivector.blade(SIG, E12, az / len));
        return Multivector.scalar(SIG, Math.cos(angle / 2))
                .sub(b.scale(Math.sin(angle / 2)));
    }

    /** Rigid motion M = T R: rotate first, then translate. */
    public static Multivector motor(Multivector rotor, Multivector translator) {
        return translator.gp(rotor);
    }

    /** M x ~M */
    public static Multivector apply(Multivector motor, Multivector x) {
        return motor.sandwich(x);
    }

    /** Join of two point trivectors: the line bivector through them. */
    public static Multivector line(Multivector p1, Multivector p2) {
        return p1.vee(p2);
    }

    /** Meet of two planes: the line bivector where they intersect. */
    public static Multivector meetPlanes(Multivector pi1, Multivector pi2) {
        return pi1.op(pi2);
    }

    /** Meet of a plane and a line: a point trivector. */
    public static Multivector meetPlaneLine(Multivector plane, Multivector line) {
        return plane.op(line);
    }

    /** Point at infinity in direction (x,y,z): P_inf = -z e012 + y e013 - x e023. */
    public static Multivector pointAtInfinity(double x, double y, double z) {
        return Multivector.blade(SIG, E012, -z)
                .add(Multivector.blade(SIG, E013, y))
                .add(Multivector.blade(SIG, E023, -x));
    }

    public static double[] pointAtInfinityCoords(Multivector p) {
        return new double[]{-p.get(E023), p.get(E013), -p.get(E012)};
    }

    /** Translation part of a motor: where it sends the origin. */
    public static double[] translationOf(Multivector motor) {
        return pointCoords(apply(motor, point(0, 0, 0)));
    }

    /** Rotation matrix of a motor; columns are the images of e1,e2,e3. */
    public static double[][] rotationMatrix(Multivector motor) {
        double[][] m = new double[3][3];
        Multivector[] dirs = {
                pointAtInfinity(1, 0, 0),
                pointAtInfinity(0, 1, 0),
                pointAtInfinity(0, 0, 1)
        };
        for (int col = 0; col < 3; col++) {
            double[] c = pointAtInfinityCoords(apply(motor, dirs[col]));
            for (int row = 0; row < 3; row++) m[row][col] = c[row];
        }
        return m;
    }

    /** Rotor part of a motor, recovered from its action on points at infinity. */
    public static Multivector rotationOf(Multivector motor) {
        return rotorFromMatrix(rotationMatrix(motor));
    }

    /** Shepperd's quaternion extraction; returns w - (x e23 - y e13 + z e12). */
    public static Multivector rotorFromMatrix(double[][] m) {
        double tr = m[0][0] + m[1][1] + m[2][2];
        double w, x, y, z;
        if (tr > 0) {
            double s = Math.sqrt(tr + 1.0) * 2;
            w = 0.25 * s;
            x = (m[2][1] - m[1][2]) / s;
            y = (m[0][2] - m[2][0]) / s;
            z = (m[1][0] - m[0][1]) / s;
        } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
            double s = Math.sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2;
            w = (m[2][1] - m[1][2]) / s;
            x = 0.25 * s;
            y = (m[0][1] + m[1][0]) / s;
            z = (m[0][2] + m[2][0]) / s;
        } else if (m[1][1] > m[2][2]) {
            double s = Math.sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2;
            w = (m[0][2] - m[2][0]) / s;
            x = (m[0][1] + m[1][0]) / s;
            y = 0.25 * s;
            z = (m[1][2] + m[2][1]) / s;
        } else {
            double s = Math.sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2;
            w = (m[1][0] - m[0][1]) / s;
            x = (m[0][2] + m[2][0]) / s;
            y = (m[1][2] + m[2][1]) / s;
            z = 0.25 * s;
        }
        return Multivector.scalar(SIG, w)
                .sub(Multivector.blade(SIG, E23, x))
                .add(Multivector.blade(SIG, E13, y))
                .sub(Multivector.blade(SIG, E12, z));
    }

    /** Interpolate rigid motions: rotation slerp + linear translation. */
    public static Multivector slerp(Multivector a, Multivector b, double t) {
        Multivector r = Gen.slerp(rotationOf(a), rotationOf(b), t);
        double[] ta = translationOf(a);
        double[] tb = translationOf(b);
        return motor(r, translator(
                ta[0] + (tb[0] - ta[0]) * t,
                ta[1] + (tb[1] - ta[1]) * t,
                ta[2] + (tb[2] - ta[2]) * t));
    }
}
