package vsr;

import java.util.Random;

/** Sparse-support fast path in the dense Multivector representation. */
public final class GpSupportBench {
    public static void main(String[] args) {
        Random rnd = new Random(8);
        Multivector motor = P3.motor(P3.rotor(0.3, -0.7, 0.2, 0.9), P3.translator(1, 2, 3));
        Multivector point = P3.point(1.5, -2.0, 0.75);
        Multivector rotor = E3.rotor(0.3, -0.7, 0.2, 0.9);
        Multivector vec = E3.vec(1, 2, 3);
        double p3 = GaBench.ns(() -> motor.gp(point).get(14));
        double e3 = GaBench.ns(() -> rotor.gp(vec).get(1));
        System.out.printf("SUPPORT,P3,motor_gp_point,%.1f%n", p3);
        System.out.printf("SUPPORT,E3,rotor_gp_vec,%.1f%n", e3);

        Signature sig = Signature.euclidean(10);
        for (int k : new int[]{2, 4, 10, 32, 128, 512, 1024}) {
            Multivector a = sparse(sig, k, rnd);
            Multivector b = sparse(sig, k, rnd);
            double ns = GaBench.ns(() -> a.gp(b).get(0));
            System.out.printf("SUPPORT,10,%d,%.1f%n", k, ns);
        }
    }

    private static Multivector sparse(Signature sig, int k, Random rnd) {
        double[] c = new double[sig.size];
        for (int i = 0; i < k; i++) c[rnd.nextInt(sig.size)] = rnd.nextDouble() * 2 - 1;
        return new Multivector(sig, c);
    }
}
