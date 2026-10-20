package vsr;

import org.joml.Matrix4f;
import org.joml.Quaternionf;
import org.joml.Vector3f;

import vsr.blaze3d.Blaze3D;

/**
 * Head-to-head comparison with JOML (float quaternions and 4x4 matrices).
 * JOML is the adapter target used by the Fabric client, so these numbers show
 * what the geometric-algebra path costs relative to the standard library.
 */
public final class JomlCompare {
    public static void main(String[] args) {
        Multivector rotorA = E3.rotor(1, 2, 3, 0.7);
        Multivector rotorB = E3.rotor(-2, 1, 0.5, -0.8);
        Quaternionf qa = Blaze3D.quaternion(rotorA);
        Quaternionf qb = Blaze3D.quaternion(rotorB);
        Vector3f v = new Vector3f(1.5f, -2.0f, 0.25f);
        Multivector e3v = E3.vec(1.5, -2.0, 0.25);

        System.out.println("## E3 vector rotation");
        System.out.printf("E3,rotate,%.2f%n", GaBench.ns(() -> rotorA.sandwich(e3v).get(0)));
        System.out.printf("JOML,rotateInPlace,%.2f%n", GaBench.ns(() -> qa.transform(v).x));
        System.out.printf("JOML,rotateAlloc,%.2f%n", GaBench.ns(() -> qa.transform(new Vector3f(v)).x));
        System.out.printf("E3,toQuaternion,%.2f%n", GaBench.ns(() -> Blaze3D.quaternion(rotorA).w));
        System.out.printf("JOML,toE3Rotor,%.2f%n", GaBench.ns(() -> Blaze3D.rotor(qa).get(0)));

        System.out.println("## rotor slerp");
        System.out.printf("E3,slerp,%.2f%n", GaBench.ns(() -> Gen.slerp(rotorA, rotorB, 0.3).get(0)));
        System.out.printf("JOML,slerpAlloc,%.2f%n", GaBench.ns(() -> new Quaternionf(qa).slerp(qb, 0.3f).w));
        Quaternionf qd = new Quaternionf();
        System.out.printf("JOML,slerpInto,%.2f%n", GaBench.ns(() -> qd.set(qa).slerp(qb, 0.3f).w));

        System.out.println("## P3 rigid transform of a point");
        Multivector motor = P3.motor(P3.rotor(0, 0, 1, 0.6), P3.translator(1, 2, 3));
        Multivector point = P3.point(1.5f, 2.5f, -0.5f);
        Matrix4f matrix = Blaze3D.matrix(motor);
        Vector3f p = new Vector3f(1.5f, 2.5f, -0.5f);
        System.out.printf("P3,applyPoint,%.2f%n", GaBench.ns(() -> {
            Multivector img = motor.sandwich(point);
            return img.get(14);
        }));
        System.out.printf("P3,applyAndProject,%.2f%n",
                GaBench.ns(() -> P3.pointCoords(motor.sandwich(point))[0]));
        System.out.printf("JOML,matrixTransformInPlace,%.2f%n",
                GaBench.ns(() -> matrix.transformPosition(p).x));
        System.out.printf("JOML,matrixTransformAlloc,%.2f%n",
                GaBench.ns(() -> matrix.transformPosition(new Vector3f(p)).x));
        System.out.printf("JOML,toMatrix4f,%.2f%n", GaBench.ns(() -> Blaze3D.matrix(motor).m00()));
        System.out.printf("JOML,fromMatrix4f,%.2f%n", GaBench.ns(() -> Blaze3D.motor(matrix).get(0)));
        System.out.printf("JOML,adapterTransform,%.2f%n",
                GaBench.ns(() -> Blaze3D.transform(motor, p, new Vector3f()).x));

        System.out.println("## accuracy (float JOML vs double GA)");
        Vector3f viaJoml = matrix.transformPosition(new Vector3f(1.5f, 2.5f, -0.5f));
        double[] viaGa = P3.pointCoords(motor.sandwich(point));
        System.out.printf("ACCURACY,maxAbs=%.3e%n",
                Math.max(Math.abs(viaJoml.x - viaGa[0]),
                        Math.max(Math.abs(viaJoml.y - viaGa[1]),
                                Math.abs(viaJoml.z - viaGa[2]))));
    }
}
