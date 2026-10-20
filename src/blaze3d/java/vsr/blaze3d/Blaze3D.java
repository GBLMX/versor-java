package vsr.blaze3d;

import org.joml.Matrix4f;
import org.joml.Quaternionf;
import org.joml.Vector3f;

import vsr.E3;
import vsr.Multivector;
import vsr.P3;
import vsr.Signature;

/**
 * Adapter between the GA core and Mojang's Blaze3D, which consumes JOML
 * types (Matrix4f, Quaternionf, Vector3f) through PoseStack.
 *
 * <p>The core stays dependency-free: this class lives in a separate source
 * root and is only compiled when JOML is on the classpath.
 *
 * <p>Rotor convention matches {@link E3}/{@link P3}:
 * R = cos(theta/2) - sin(theta/2)(ax e23 + ay e31 + az e12), so the
 * quaternion is (x,y,z,w) = (-R_e23, R_e13, -R_e12, R_1).
 */
public final class Blaze3D {
    private Blaze3D() {}

    /** JOML quaternion for an E3 or P3 rotor. */
    public static Quaternionf quaternion(Multivector rotor) {
        Signature sig = rotor.sig;
        int m12 = blade(sig, "e1e2");
        int m13 = blade(sig, "e1e3");
        int m23 = blade(sig, "e2e3");
        double w = rotor.scalarPart();
        double x = -rotor.get(m23);
        double y = rotor.get(m13);
        double z = -rotor.get(m12);
        double n = Math.sqrt(w * w + x * x + y * y + z * z);
        if (n < 1e-12) throw new ArithmeticException("degenerate rotor");
        return new Quaternionf(
                (float) (x / n), (float) (y / n), (float) (z / n), (float) (w / n));
    }

    /** E3 rotor for a JOML quaternion. */
    public static Multivector rotor(Quaternionf q) {
        return rotor(E3.SIG, q);
    }

    /** Rotor in the given 3D algebra (E3 or P3) for a JOML quaternion. */
    public static Multivector rotor(Signature sig, Quaternionf q) {
        double n = Math.sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
        if (n < 1e-12) throw new ArithmeticException("degenerate quaternion");
        double x = q.x / n, y = q.y / n, z = q.z / n, w = q.w / n;
        return Multivector.scalar(sig, w)
                .sub(Multivector.blade(sig, blade(sig, "e2e3"), x))
                .add(Multivector.blade(sig, blade(sig, "e1e3"), y))
                .sub(Multivector.blade(sig, blade(sig, "e1e2"), z));
    }

    /** P3 motor (rotation + translation) for a JOML quaternion and offset. */
    public static Multivector motor(Quaternionf q, Vector3f translation) {
        return P3.motor(rotor(P3.SIG, q),
                P3.translator(translation.x, translation.y, translation.z));
    }

    /** 4x4 rigid transform (rotation + translation), column-major like OpenGL. */
    public static Matrix4f matrix(Multivector motor) {
        double[][] r = P3.rotationMatrix(motor);
        double[] t = P3.translationOf(motor);
        float[] a = {
                (float) r[0][0], (float) r[1][0], (float) r[2][0], 0f,
                (float) r[0][1], (float) r[1][1], (float) r[2][1], 0f,
                (float) r[0][2], (float) r[1][2], (float) r[2][2], 0f,
                (float) t[0], (float) t[1], (float) t[2], 1f
        };
        return new Matrix4f().set(a);
    }

    /** P3 motor for a rigid transform matrix (rotation must be orthonormal). */
    public static Multivector motor(Matrix4f m) {
        double[][] r = {
                {m.m00(), m.m01(), m.m02()},
                {m.m10(), m.m11(), m.m12()},
                {m.m20(), m.m21(), m.m22()}
        };
        Vector3f t = m.getTranslation(new Vector3f());
        return P3.motor(P3.rotorFromMatrix(r), P3.translator(t.x, t.y, t.z));
    }

    /** Translation part of a motor as a JOML vector. */
    public static Vector3f translation(Multivector motor) {
        double[] t = P3.translationOf(motor);
        return new Vector3f((float) t[0], (float) t[1], (float) t[2]);
    }

    /** Transform a vector by a rigid motor (rotation then translation). */
    public static Vector3f transform(Multivector motor, Vector3f v) {
        return transform(motor, v, new Vector3f());
    }

    public static Vector3f transform(Multivector motor, Vector3f v, Vector3f dest) {
        double[][] r = P3.rotationMatrix(motor);
        double[] t = P3.translationOf(motor);
        dest.x = (float) (r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z + t[0]);
        dest.y = (float) (r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z + t[1]);
        dest.z = (float) (r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z + t[2]);
        return dest;
    }

    private static int blade(Signature sig, String name) {
        for (int m = 0; m < sig.size; m++) {
            if (sig.bladeName(m).equals(name)) return m;
        }
        throw new IllegalArgumentException("missing blade " + name + " in " + sig);
    }
}
