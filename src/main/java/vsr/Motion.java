package vsr;

/**
 * Rigid-body kinematics in PGA3: M(t) = exp(t * twist), where the twist is
 * a bivector generator (rotation and/or translation).
 *
 * <p>For a point P(t) = M p ~M the derivatives are
 * <pre>
 *   P'  = W P - P W
 *   P'' = W^2 P + P W^2 - 2 W P W
 *   P'''= W^3 P - 3 W^2 P W + 3 W P W^2 - P W^3
 * </pre>
 * with W the constant twist, so position, velocity, acceleration and jerk
 * are exact rather than finite differences.
 */
public final class Motion {
    private final Multivector twist;
    private final Multivector twist2;
    private final Multivector twist3;

    public Motion(Multivector twist) {
        if (twist.grade(2).distanceTo(twist) > 1e-12) {
            throw new IllegalArgumentException("twist must be a bivector");
        }
        this.twist = twist;
        this.twist2 = twist.gp(twist);
        this.twist3 = twist2.gp(twist);
    }

    /** Rotation about the origin with angular velocity omega. */
    public static Motion rotation(double ax, double ay, double az, double omega) {
        // log of a unit rotor is -0.5 * (unit bivector), so scale by omega
        Multivector unit = Gen.log(P3.rotor(ax, ay, az, 1.0));
        return new Motion(unit.scale(omega));
    }

    /** Constant linear velocity (a pure translation generator). */
    public static Motion translation(double vx, double vy, double vz) {
        return new Motion(P3.translator(vx, vy, vz)
                .sub(Multivector.scalar(P3.SIG, 1)));
    }

    /** Rotation about the origin plus constant translation. */
    public static Motion screw(double ax, double ay, double az, double omega,
                               double vx, double vy, double vz) {
        Multivector rot = Gen.log(P3.rotor(ax, ay, az, 1.0)).scale(omega);
        Multivector trans = P3.translator(vx, vy, vz)
                .sub(Multivector.scalar(P3.SIG, 1));
        return new Motion(rot.add(trans));
    }

    public Multivector twist() {
        return twist;
    }

    /** Pose at time t. */
    public Multivector motor(double t) {
        return twist.scale(t).exp();
    }

    /** Transformed point at time t (a trivector). */
    public Multivector pointAt(double t, double x, double y, double z) {
        return P3.apply(motor(t), P3.point(x, y, z));
    }

    public double[] position(double t, double x, double y, double z) {
        return P3.pointCoords(pointAt(t, x, y, z));
    }

    public double[] velocity(double t, double x, double y, double z) {
        Multivector p = pointAt(t, x, y, z);
        return P3.pointAtInfinityCoords(commutator(twist, p));
    }

    public double[] acceleration(double t, double x, double y, double z) {
        Multivector p = pointAt(t, x, y, z);
        Multivector a = twist2.gp(p).add(p.gp(twist2))
                .sub(twist.gp(p).gp(twist).scale(2));
        return P3.pointAtInfinityCoords(a);
    }

    public double[] jerk(double t, double x, double y, double z) {
        Multivector p = pointAt(t, x, y, z);
        Multivector j = twist3.gp(p)
                .sub(twist2.gp(p).gp(twist).scale(3))
                .add(twist.gp(p).gp(twist2).scale(3))
                .sub(p.gp(twist3));
        return P3.pointAtInfinityCoords(j);
    }

    private static Multivector commutator(Multivector a, Multivector b) {
        return a.gp(b).sub(b.gp(a));
    }
}
