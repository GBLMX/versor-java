package vsr.blaze3d;

import org.joml.Matrix4f;
import org.joml.Quaternionf;
import org.joml.Vector3f;

import vsr.Multivector;
import vsr.P3;

/** Verifies the JOML adapter against the GA core and Blaze3D conventions. */
public final class Blaze3DTest {
    private static int passed = 0;

    public static void main(String[] args) {
        Multivector motor = P3.motor(
                P3.rotor(0, 0, 1, Math.PI / 2),
                P3.translator(1, 2, 3));

        Matrix4f m = Blaze3D.matrix(motor);
        Vector3f v = Blaze3D.transform(motor, new Vector3f(1, 0, 0));
        Vector3f byJoml = m.transformPosition(new Vector3f(1, 0, 0));
        check("transform matches matrix", v.distance(byJoml) < 1e-6f);
        check("transform value", Math.abs(v.x - 1) < 1e-6
                && Math.abs(v.y - 3) < 1e-6 && Math.abs(v.z - 3) < 1e-6);

        double[] core = P3.pointCoords(P3.apply(motor, P3.point(1, 0, 0)));
        check("matches core", Math.abs(core[0] - v.x) < 1e-6
                && Math.abs(core[1] - v.y) < 1e-6 && Math.abs(core[2] - v.z) < 1e-6);

        // JOML names elements m<column><row>; col 0 is the image of e1 = (0,1,0)
        check("column-major layout", Math.abs(m.m00()) < 1e-6
                && Math.abs(m.m01() - 1) < 1e-6
                && Math.abs(m.m10() + 1) < 1e-6
                && Math.abs(m.m30() - 1) < 1e-6
                && Math.abs(m.m31() - 2) < 1e-6
                && Math.abs(m.m32() - 3) < 1e-6);

        Quaternionf q = Blaze3D.quaternion(P3.rotationOf(motor));
        check("quaternion components", Math.abs(q.x) < 1e-6
                && Math.abs(q.y) < 1e-6
                && Math.abs(q.z - Math.sqrt(0.5)) < 1e-6
                && Math.abs(q.w - Math.sqrt(0.5)) < 1e-6);

        Multivector fromQ = Blaze3D.rotor(P3.SIG, q);
        double[] rotated = P3.pointCoords(P3.apply(fromQ, P3.point(1, 0, 0)));
        check("quaternion round-trip", Math.abs(rotated[0]) < 1e-6
                && Math.abs(rotated[1] - 1) < 1e-6);

        Multivector fromM = Blaze3D.motor(m);
        double[] p = P3.pointCoords(P3.apply(fromM, P3.point(0, 0, 0)));
        check("matrix round-trip", Math.abs(p[0] - 1) < 1e-6
                && Math.abs(p[1] - 2) < 1e-6 && Math.abs(p[2] - 3) < 1e-6);

        Multivector fromQt = Blaze3D.motor(q, new Vector3f(1, 2, 3));
        double[] qt = P3.pointCoords(P3.apply(fromQt, P3.point(1, 0, 0)));
        check("quaternion+offset motor", Math.abs(qt[0] - 1) < 1e-6
                && Math.abs(qt[1] - 3) < 1e-6 && Math.abs(qt[2] - 3) < 1e-6);

        Vector3f t = Blaze3D.translation(motor);
        check("translation", Math.abs(t.x - 1) < 1e-6
                && Math.abs(t.y - 2) < 1e-6 && Math.abs(t.z - 3) < 1e-6);

        Multivector pureT = Blaze3D.motor(new Matrix4f().translation(4, 5, 6));
        double[] pt = P3.pointCoords(P3.apply(pureT, P3.point(0, 0, 0)));
        check("pure translation matrix", Math.abs(pt[0] - 4) < 1e-6
                && Math.abs(pt[1] - 5) < 1e-6 && Math.abs(pt[2] - 6) < 1e-6);

        Quaternionf qi = Blaze3D.quaternion(Multivector.scalar(P3.SIG, 1));
        check("identity quaternion", Math.abs(qi.x) < 1e-6 && Math.abs(qi.y) < 1e-6
                && Math.abs(qi.z) < 1e-6 && Math.abs(qi.w - 1) < 1e-6);

        System.out.println("OK: " + passed + " Blaze3D checks passed");
    }

    private static void check(String name, boolean ok) {
        if (!ok) throw new AssertionError("FAILED: " + name);
        passed++;
        System.out.println("ok  " + name);
    }
}
