package vsr;

/** 3D Euclidean GA (quaternion / rotor model). */
public final class E3 {
    public static final Signature SIG = Signature.e3();

    private E3() {}

    public static Multivector vec(double x, double y, double z) {
        return Multivector.vector(SIG, x, y, z);
    }

    /** Bivector for an axis: dual(a) = a1 e23 + a2 e31 + a3 e12. */
    public static Multivector bivector(double ax, double ay, double az) {
        return Multivector.blade(SIG, 0b110, ax)   // e2e3
                .add(Multivector.blade(SIG, 0b101, -ay)) // e3e1 = -e1e3
                .add(Multivector.blade(SIG, 0b011, az)); // e1e2
    }

    /** R = exp(-theta/2 * B) for the unit axis (ax,ay,az). */
    public static Multivector rotor(double ax, double ay, double az, double angle) {
        double len = Math.sqrt(ax * ax + ay * ay + az * az);
        if (len == 0) throw new IllegalArgumentException("zero axis");
        Multivector b = bivector(ax / len, ay / len, az / len);
        return Multivector.scalar(SIG, Math.cos(angle / 2))
                .sub(b.scale(Math.sin(angle / 2)));
    }

    /** Rotate a vector (or any multivector) with R x ~R. */
    public static Multivector rotate(Multivector x, Multivector rotor) {
        return rotor.sandwich(x);
    }

    public static double[] coords(Multivector v) {
        return new double[]{v.get(1), v.get(2), v.get(4)};
    }
}
