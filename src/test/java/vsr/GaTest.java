package vsr;

import java.util.Random;

/**
 * Self-contained test harness (no external test framework available on this
 * machine). Exits with status 1 on the first failed assertion.
 */
public final class GaTest {
    private static int passed = 0;

    public static void main(String[] args) {
        testSignatureFactories();
        testE3Basis();
        testProducts();
        testUnaryOps();
        testCoords();
        testReverseAndDual();
        testRotor();
        testSlerpAndLog();
        testStaBoostInterpolation();
        testPgaPoint();
        testPgaTranslator();
        testPgaRotor();
        testPgaMotorComposition();
        testPgaIntersections();
        testKinematics();
        testFieldsAndManifolds();
        testWolframGolden();
        testMotorInterpolation();
        testConformal();
        testSandwichFastPaths();
        System.out.println("OK: " + passed + " checks passed");
    }

    private static void testSignatureFactories() {
        Signature sta = Signature.sta();
        check("STA metric", sta.bladeProductSign(1, 1) == 1
                && sta.bladeProductSign(2, 2) == -1);
        check("STA names", sta.bladeName(1).equals("e0"));
        Signature cga = Signature.cga3();
        check("CGA metric", cga.bladeProductSign(1 << 4, 1 << 4) == -1
                && cga.bladeProductSign(1, 1) == 1);
        check("CGA names", cga.bladeName(1 << 4).equals("e5"));
        check("custom signature", Signature.of(1, 1, 1).bladeProductSign(1, 1) == 1);
        check("Euclidean names", Signature.euclidean(2).bladeName(2).equals("e2"));
    }

    private static void testE3Basis() {
        Multivector e1 = Multivector.blade(E3.SIG, 1, 1);
        Multivector e2 = Multivector.blade(E3.SIG, 2, 1);
        check("e1*e1 = 1", e1.gp(e1).distanceTo(Multivector.scalar(E3.SIG, 1)) < 1e-12);
        check("e1*e2 = -e2*e1", e1.gp(e2).distanceTo(e2.gp(e1).neg()) < 1e-12);
        check("e1*e2*e2 = e1", e1.gp(e2).gp(e2).distanceTo(e1) < 1e-12);
    }

    private static void testProducts() {
        Random rnd = new Random(42);
        for (int trial = 0; trial < 20; trial++) {
            Multivector a = random(E3.SIG, rnd);
            Multivector b = random(E3.SIG, rnd);
            Multivector c = random(E3.SIG, rnd);
            check("associativity", a.gp(b).gp(c).distanceTo(a.gp(b.gp(c))) < 1e-10);
            check("(ab)~ = b~ a~",
                    a.gp(b).reverse().distanceTo(b.reverse().gp(a.reverse())) < 1e-10);
            Multivector e1 = Multivector.blade(E3.SIG, 1, 1);
            Multivector e2 = Multivector.blade(E3.SIG, 2, 1);
            check("e1∧e1 = 0", e1.op(e1).isZero(1e-12));
            check("e1∧e2 = -e2∧e1", e1.op(e2).distanceTo(e2.op(e1).neg()) < 1e-12);
        }
    }

    private static void testUnaryOps() {
        Random rnd = new Random(7);
        Multivector a = random(E3.SIG, rnd);
        check("even + odd = a", a.even().add(a.odd()).distanceTo(a) < 1e-12);
        Multivector sum = Multivector.zero(E3.SIG);
        for (int k = 0; k <= 3; k++) sum = sum.add(a.grade(k));
        check("grade projections sum to a", sum.distanceTo(a) < 1e-12);
        check("conjugate = reverse then involution",
                a.conjugate().distanceTo(a.reverse().involution()) < 1e-12);
        check("unary operations commute",
                a.reverse().involution().distanceTo(a.involution().reverse()) < 1e-12);
    }

    private static void testCoords() {
        double[] p = P3.planeCoords(P3.plane(1, 2, 3, -4));
        check("plane coords", Math.abs(p[0] - 1) < 1e-12 && Math.abs(p[1] - 2) < 1e-12
                && Math.abs(p[2] - 3) < 1e-12 && Math.abs(p[3] + 4) < 1e-12);
        double[] v = E3.coords(E3.vec(4, 5, 6));
        check("E3 vector coords", Math.abs(v[0] - 4) < 1e-12
                && Math.abs(v[1] - 5) < 1e-12 && Math.abs(v[2] - 6) < 1e-12);
    }

    private static void testReverseAndDual() {
        Multivector e1 = Multivector.blade(E3.SIG, 1, 1);
        Multivector e2 = Multivector.blade(E3.SIG, 2, 1);
        Multivector I = Multivector.blade(E3.SIG, 7, 1);
        check("e1 dual", e1.dual().distanceTo(e2.gp(Multivector.blade(E3.SIG, 4, 1))) < 1e-12
                || e1.dual().distanceTo(Multivector.blade(E3.SIG, 0b110, 1)) < 1e-12);
        check("e1 ∧ dual(e1) = I", e1.op(e1.dual()).distanceTo(I) < 1e-12
                || e1.op(e1.dual()).distanceTo(I.neg()) < 1e-12);
    }

    private static void testRotor() {
        Multivector e1 = E3.vec(1, 0, 0);
        Multivector e2 = E3.vec(0, 1, 0);
        Multivector r = E3.rotor(0, 0, 1, Math.PI / 2);
        check("R e1 R~ = e2", E3.rotate(e1, r).distanceTo(e2) < 1e-12);
        check("R e2 R~ = -e1", E3.rotate(e2, r).distanceTo(e1.neg()) < 1e-12);
        check("R ~R = 1", r.gp(r.reverse()).distanceTo(Multivector.scalar(E3.SIG, 1)) < 1e-12);
        check("R norm 1", Math.abs(r.norm() - 1) < 1e-12);
    }

    private static void testSlerpAndLog() {
        Multivector a = E3.rotor(0, 0, 1, 0);
        Multivector b = E3.rotor(0, 0, 1, Math.PI / 2);
        Multivector mid = Gen.slerp(a, b, 0.5);
        Multivector expect = E3.rotor(0, 0, 1, Math.PI / 4);
        check("slerp midpoint", mid.distanceTo(expect) < 1e-12);
        check("exp(log R) = R",
                Gen.expBivector(Gen.log(b)).distanceTo(b) < 1e-12);
        Multivector minus = Multivector.scalar(E3.SIG, -1);
        Multivector halfTurn = Gen.slerp(Multivector.scalar(E3.SIG, 1), minus, 0.5);
        check("antipodal slerp", halfTurn.distanceTo(Multivector.blade(E3.SIG, 3, 1)) < 1e-12);
        Multivector gen = E3.bivector(1, 2, -1).scale(0.3);
        check("series exp = closed form",
                Gen.expBivector(gen).distanceTo(gen.exp()) < 1e-12);
        check("exp(0) = 1",
                Multivector.zero(E3.SIG).exp()
                        .distanceTo(Multivector.scalar(E3.SIG, 1)) < 1e-12);
    }

    private static void testStaBoostInterpolation() {
        Signature sta = Signature.sta();
        Multivector identity = Multivector.scalar(sta, 1);
        Multivector boost = Multivector.blade(sta, 0b0011, 0.7); // 0.7 e0e1
        Multivector rotor = Gen.expBivector(boost);
        check("STA boost log", Gen.log(rotor).distanceTo(boost) < 1e-10);
        check("STA boost slerp", Gen.slerp(identity, rotor, 0.5)
                .distanceTo(Gen.expBivector(boost.scale(0.5))) < 1e-10);
    }

    private static void testPgaPoint() {
        Multivector p = P3.point(1, 2, 3);
        double[] xyz = P3.pointCoords(p);
        check("point round-trip", Math.abs(xyz[0] - 1) < 1e-12
                && Math.abs(xyz[1] - 2) < 1e-12 && Math.abs(xyz[2] - 3) < 1e-12);
        check("point is trivector", p.grade(3).distanceTo(p) < 1e-12);
        Multivector d = p.dual();
        check("dual point is plane (vector)", d.grade(1).distanceTo(d) < 1e-12);
    }

    private static void testPgaTranslator() {
        Multivector origin = P3.point(0, 0, 0);
        Multivector t = P3.translator(1, 2, 3);
        Multivector moved = P3.apply(t, origin);
        double[] xyz = P3.pointCoords(moved);
        check("translator moves origin", Math.abs(xyz[0] - 1) < 1e-12
                && Math.abs(xyz[1] - 2) < 1e-12 && Math.abs(xyz[2] - 3) < 1e-12);
        check("translator = exp(bivector)",
                t.sub(1).exp().distanceTo(t) < 1e-12);
        check("translator is unit versor",
                t.gp(t.reverse()).distanceTo(Multivector.scalar(P3.SIG, 1)) < 1e-12);
        Multivector p = P3.point(4, 5, 6);
        double[] q = P3.pointCoords(P3.apply(t, p));
        check("translator moves point", Math.abs(q[0] - 5) < 1e-12
                && Math.abs(q[1] - 7) < 1e-12 && Math.abs(q[2] - 9) < 1e-12);
    }

    private static void testPgaRotor() {
        Multivector p = P3.point(1, 0, 0);
        Multivector r = P3.rotor(0, 0, 1, Math.PI / 2);
        double[] q = P3.pointCoords(P3.apply(r, p));
        check("PGA rotor z 90deg", Math.abs(q[0]) < 1e-12
                && Math.abs(q[1] - 1) < 1e-12 && Math.abs(q[2]) < 1e-12);
        check("PGA rotor unit", r.gp(r.reverse())
                .distanceTo(Multivector.scalar(P3.SIG, 1)) < 1e-12);
    }

    private static void testPgaMotorComposition() {
        Multivector p = P3.point(1, 0, 0);
        Multivector r = P3.rotor(0, 0, 1, Math.PI / 2);
        Multivector t = P3.translator(0, 2, 0);
        Multivector m = P3.motor(r, t);
        double[] q = P3.pointCoords(P3.apply(m, p));
        // rotate (1,0,0) -> (0,1,0), then translate by (0,2,0) -> (0,3,0)
        check("motor = T R", Math.abs(q[0]) < 1e-12
                && Math.abs(q[1] - 3) < 1e-12 && Math.abs(q[2]) < 1e-12);
        Multivector back = P3.apply(m.inverse(), P3.apply(m, p));
        double[] b = P3.pointCoords(back);
        check("motor inverse", Math.abs(b[0] - 1) < 1e-12
                && Math.abs(b[1]) < 1e-12 && Math.abs(b[2]) < 1e-12);
    }

    private static void testPgaIntersections() {
        Multivector px = P3.plane(1, 0, 0, 0);   // x = 0
        Multivector py = P3.plane(0, 1, 0, 0);   // y = 0
        Multivector line = P3.meetPlanes(px, py); // z axis
        check("meet of planes is a bivector", line.grade(2).distanceTo(line) < 1e-12);

        Multivector pz = P3.plane(0, 0, 1, -2);  // z = 2
        double[] hit = P3.pointCoords(P3.meetPlaneLine(pz, line));
        check("plane/line meet", Math.abs(hit[0]) < 1e-12
                && Math.abs(hit[1]) < 1e-12 && Math.abs(hit[2] - 2) < 1e-12);

        Multivector a = P3.point(1, 2, 3);
        Multivector b = P3.point(4, 5, 6);
        Multivector l = P3.line(a, b);
        check("point join is a bivector", l.grade(2).distanceTo(l) < 1e-12);

        Multivector z3 = P3.plane(0, 0, 1, -3);
        double[] p1 = P3.pointCoords(P3.meetPlaneLine(z3, l));
        check("line through points meets plane", Math.abs(p1[0] - 1) < 1e-12
                && Math.abs(p1[1] - 2) < 1e-12 && Math.abs(p1[2] - 3) < 1e-12);
    }

    private static void testKinematics() {
        Motion spin = Motion.rotation(0, 0, 1, 1.0);
        for (double t : new double[]{0.0, 0.7, 1.2}) {
            double c = Math.cos(t), sn = Math.sin(t);
            double[] p = spin.position(t, 1, 0, 0);
            double[] v = spin.velocity(t, 1, 0, 0);
            double[] a = spin.acceleration(t, 1, 0, 0);
            double[] j = spin.jerk(t, 1, 0, 0);
            check("rotation position", near(p, new double[]{c, sn, 0}));
            check("rotation velocity", near(v, new double[]{-sn, c, 0}));
            check("rotation acceleration", near(a, new double[]{-c, -sn, 0}));
            check("rotation jerk", near(j, new double[]{sn, -c, 0}));
        }

        Motion drift = Motion.translation(2, 1, 0);
        check("translation velocity", near(drift.velocity(3.0, 1, 1, 1), new double[]{2, 1, 0}));
        check("translation acceleration",
                near(drift.acceleration(3.0, 1, 1, 1), new double[]{0, 0, 0}));
        check("translation position",
                near(drift.position(3.0, 1, 1, 1), new double[]{7, 4, 1}));

        Motion screw = Motion.screw(0, 0, 1, 1.0, 0, 0, 1);
        double t = 0.8;
        double c = Math.cos(t), sn = Math.sin(t);
        check("screw position", near(screw.position(t, 1, 0, 0), new double[]{c, sn, t}));
        check("screw velocity", near(screw.velocity(t, 1, 0, 0),
                new double[]{-sn, c, 1}));
        check("screw acceleration", near(screw.acceleration(t, 1, 0, 0),
                new double[]{-c, -sn, 0}));

        boolean threw = false;
        try {
            Motion.rotation(0, 0, 0, 1);
        } catch (IllegalArgumentException e) {
            threw = true;
        }
        check("zero axis rejected", threw);

        threw = false;
        try {
            new Motion(Multivector.scalar(P3.SIG, 1));
        } catch (IllegalArgumentException e) {
            threw = true;
        }
        check("non-bivector twist rejected", threw);

        check("zero twist is identity",
                near(Motion.rotation(0, 0, 1, 0).velocity(1.0, 1, 2, 3), new double[]{0, 0, 0}));
        check("zero twist keeps position",
                near(Motion.rotation(0, 0, 1, 0).position(1.0, 1, 2, 3), new double[]{1, 2, 3}));

        // finite-difference cross-check on a general screw
        Motion gen = Motion.screw(1, 2, 3, 0.8, 0.3, -0.2, 0.5);
        double h = 1e-4;
        double x = 1, y = -1, z = 0.5;
        double[] px = gen.position(0.6 + h, x, y, z);
        double[] pm = gen.position(0.6 - h, x, y, z);
        double[] fd = {(px[0] - pm[0]) / (2 * h),
                (px[1] - pm[1]) / (2 * h), (px[2] - pm[2]) / (2 * h)};
        check("velocity matches finite difference",
                near(gen.velocity(0.6, x, y, z), fd, 1e-5));

        double[] p0 = gen.position(0.6, x, y, z);
        double[] fd2 = {(px[0] - 2 * p0[0] + pm[0]) / (h * h),
                (px[1] - 2 * p0[1] + pm[1]) / (h * h),
                (px[2] - 2 * p0[2] + pm[2]) / (h * h)};
        check("acceleration matches finite difference",
                near(gen.acceleration(0.6, x, y, z), fd2, 1e-4));

        double[] aPlus = gen.acceleration(0.6 + h, x, y, z);
        double[] aMinus = gen.acceleration(0.6 - h, x, y, z);
        double[] fd3 = {(aPlus[0] - aMinus[0]) / (2 * h),
                (aPlus[1] - aMinus[1]) / (2 * h), (aPlus[2] - aMinus[2]) / (2 * h)};
        check("jerk matches finite difference",
                near(gen.jerk(0.6, x, y, z), fd3, 1e-5));
    }

    private static void testFieldsAndManifolds() {
        // F(x, y, z) = (-y, x, 0): divergence 0, curl 2 e3, nabla F = 2 e12
        Field swirl = (x, y, z) -> E3.vec(-y, x, 0);
        check("swirl divergence", Math.abs(Calculus.divergence(swirl, 1, 2, 3)) < 1e-6);
        check("swirl curl", near(E3.coords(Calculus.curl(swirl, 1, 2, 3)),
                new double[]{0, 0, 2}, 1e-6));
        Multivector nablaF = Calculus.vectorDerivative(swirl, 1, 2, 3);
        check("swirl vector derivative", Math.abs(nablaF.get(0b011) - 2) < 1e-6
                && Math.abs(nablaF.scalarPart()) < 1e-6);

        Field scalar = (x, y, z) -> Multivector.scalar(E3.SIG, x * x + 2 * y * y + 3 * z * z);
        check("gradient", near(E3.coords(Calculus.gradient(scalar, 1, 2, 3)),
                new double[]{2, 8, 18}, 1e-5));
        check("directional", Math.abs(Calculus.directional(scalar, 1, 2, 3, 1, 0, 0).scalarPart() - 2) < 1e-5);

        Manifold plane = Manifold.surface(u -> new double[]{u[0], u[1], 0});
        check("plane normal", near(E3.coords(plane.normal(1, 2)), new double[]{0, 0, 1}, 1e-6));
        check("plane tangent blade", Math.abs(plane.tangentBlade(1, 2).get(0b011) - 1) < 1e-9);

        Manifold sphere = Manifold.surface(u -> new double[]{
                2 * Math.cos(u[0]) * Math.cos(u[1]),
                2 * Math.sin(u[0]) * Math.cos(u[1]),
                2 * Math.sin(u[1])});
        check("sphere normal", near(E3.coords(sphere.normal(0, 0)),
                new double[]{1, 0, 0}, 1e-5));

        check("plane jacobian", Math.abs(plane.jacobian(1, 2) - 1) < 1e-9);
        check("sphere jacobian", Math.abs(sphere.jacobian(0, 0) - 4) < 1e-4);

        Field radial = (x, y, z) -> E3.vec(x, y, z);
        check("sphere flux", Math.abs(Calculus.flux(radial, sphere,
                0, 2 * Math.PI, -Math.PI / 2, Math.PI / 2, 128, 64) - 32 * Math.PI) < 0.05);
        Field up = (x, y, z) -> E3.vec(0, 0, 1);
        check("plane flux", Math.abs(Calculus.flux(up, plane, 0, 1, 0, 1, 16, 16) - 1) < 1e-9);

        Manifold helix = Manifold.curve(u -> new double[]{Math.cos(u[0]), Math.sin(u[0]), u[0]});
        double s = 1 / Math.sqrt(2);
        check("helix tangent", near(E3.coords(helix.tangent(0)), new double[]{0, s, s}, 1e-5));
    }

    private static boolean near(double[] a, double[] b) {
        return near(a, b, 1e-9);
    }

    private static boolean near(double[] a, double[] b, double tol) {
        for (int i = 0; i < a.length; i++) {
            if (Math.abs(a[i] - b[i]) > tol) return false;
        }
        return true;
    }

    /** Values cross-checked with galgebra (PyPI) and Wolfram CliffordAlgebra. */
    private static void testWolframGolden() {
        double[][] golden = {
                {1, 2, -3, 1.2, 0.5, 0, 0.7, 0},
                {4, -1, 2, 0, 0.25, -0.8, 0, 0},
                {-3.875, 9.8, -8.625, 5.24, -0.75, 0.2, 2.01, -2.8}
        };
        Multivector a = fromMasks(E3.SIG, golden[0]);
        Multivector b = fromMasks(E3.SIG, golden[1]);
        Multivector ab = fromMasks(E3.SIG, golden[2]);
        check("Wolfram product", a.gp(b).distanceTo(ab) < 1e-12);

        Multivector r = E3.rotor(1, 2, 3, 0.7);
        double[] rGolden = {0.9393727128473789, 0, 0, -0.2749298816087739,
                0, 0.1832865877391826, -0.0916432938695913, 0};
        check("Wolfram rotor", r.distanceTo(fromMasks(E3.SIG, rGolden)) < 1e-12);

        double[] e1img = {0, 0.781639173907025, 0.5501172307043584, 0, -0.2939578784385805, 0, 0, 0};
        double[] e2img = {0, -0.4829292842142122, 0.8320301337746345, 0, 0.27295633888831433, 0, 0, 0};
        check("Wolfram R e1 R~",
                E3.rotate(E3.vec(1, 0, 0), r).distanceTo(fromMasks(E3.SIG, e1img)) < 1e-12);
        check("Wolfram R e2 R~",
                E3.rotate(E3.vec(0, 1, 0), r).distanceTo(fromMasks(E3.SIG, e2img)) < 1e-12);
    }

    private static void testMotorInterpolation() {
        Multivector identity = Multivector.scalar(P3.SIG, 1);
        Multivector m1 = P3.motor(
                P3.rotor(0, 0, 1, Math.PI / 2),
                P3.translator(1, 2, 3));
        Multivector mid = P3.slerp(identity, m1, 0.5);
        double[] q = P3.pointCoords(P3.apply(mid, P3.point(1, 0, 0)));
        double s = Math.sqrt(0.5);
        check("motor slerp midpoint", Math.abs(q[0] - (s + 0.5)) < 1e-12
                && Math.abs(q[1] - (s + 1.0)) < 1e-12
                && Math.abs(q[2] - 1.5) < 1e-12);
        check("translationOf", Math.abs(P3.translationOf(m1)[0] - 1) < 1e-12
                && Math.abs(P3.translationOf(m1)[1] - 2) < 1e-12
                && Math.abs(P3.translationOf(m1)[2] - 3) < 1e-12);
        double[] rot = P3.pointCoords(P3.apply(P3.rotationOf(m1), P3.point(1, 0, 0)));
        check("rotationOf", Math.abs(rot[0]) < 1e-12 && Math.abs(rot[1] - 1) < 1e-12);
        Multivector pureTranslation = P3.slerp(
                P3.translator(0, 0, 0), P3.translator(2, 4, 6), 0.5);
        double[] halfT = P3.translationOf(pureTranslation);
        check("translation slerp", Math.abs(halfT[0] - 1) < 1e-12
                && Math.abs(halfT[1] - 2) < 1e-12 && Math.abs(halfT[2] - 3) < 1e-12);
    }

    private static Multivector fromMasks(Signature sig, double[] coeffs) {
        Multivector m = Multivector.zero(sig);
        for (int i = 0; i < coeffs.length; i++) {
            m = m.add(Multivector.blade(sig, i, coeffs[i]));
        }
        return m;
    }

    private static Multivector random(Signature sig, Random rnd) {
        Multivector m = Multivector.zero(sig);
        for (int i = 0; i < sig.size; i++) {
            if (rnd.nextBoolean()) m = m.add(Multivector.blade(sig, i, rnd.nextGaussian()));
        }
        return m;
    }

    private static void testConformal() {
        Multivector p = Conformal.point(1, 2, 3);
        check("CGA point is null", Math.abs(Conformal.dot(p, p)) < 1e-9);
        double[] coords = Conformal.pointCoords(p);
        check("CGA point roundtrip",
                Math.abs(coords[0] - 1) < 1e-9
                        && Math.abs(coords[1] - 2) < 1e-9
                        && Math.abs(coords[2] - 3) < 1e-9);
        check("CGA distance squared",
                Math.abs(Conformal.distanceSq(
                        Conformal.point(1, 2, 3), Conformal.point(4, 6, 3)) - 25.0) < 1e-9);

        Multivector sphere = Conformal.sphere(1, 1, 1, 2);
        check("CGA sphere radius", Math.abs(Conformal.dot(sphere, sphere) - 4.0) < 1e-9);
        check("CGA sphere surface point",
                Math.abs(Conformal.dot(Conformal.point(3, 1, 1), sphere)) < 1e-9);
        check("CGA sphere interior point",
                Math.abs(Conformal.dot(Conformal.point(1, 1, 1), sphere) - 2.0) < 1e-9);

        Multivector plane = Conformal.plane(0, 0, 1, -2);
        check("CGA plane incidence",
                Math.abs(Conformal.dot(Conformal.point(1, 1, 2), plane)) < 1e-9);
        check("CGA plane offset",
                Math.abs(Conformal.dot(Conformal.point(1, 1, 3), plane) - 1.0) < 1e-9);

        Multivector circle = Conformal.circle(0, 0, 2, 0, 0, 1, 1);
        check("CGA circle grade", circle.grade(3).distanceTo(circle) < 1e-12);
        check("CGA circle surface point",
                Conformal.point(1, 0, 2).op(circle).maxAbs() < 1e-9);
        check("CGA circle off point",
                Conformal.point(2, 0, 2).op(circle).maxAbs() > 0.5);
        check("CGA point at infinity rejected", throwsIllegal(() ->
                Conformal.pointCoords(Conformal.point(1, 1, 1).add(
                        Multivector.blade(Conformal.SIG, 8, 1)))));
    }

    private static void testSandwichFastPaths() {
        Random rnd = new Random(20251021);
        double e3Error = 0;
        double pgaError = 0;
        for (int i = 0; i < 50; i++) {
            Multivector rotor = Multivector.scalar(E3.SIG, rnd.nextDouble() * 2 - 1)
                    .add(E3.bivector(rnd.nextDouble() * 2 - 1,
                            rnd.nextDouble() * 2 - 1, rnd.nextDouble() * 2 - 1));
            Multivector vector = E3.vec(rnd.nextDouble() * 2 - 1,
                    rnd.nextDouble() * 2 - 1, rnd.nextDouble() * 2 - 1);
            e3Error = Math.max(e3Error,
                    rotor.sandwich(vector).distanceTo(rotor.gp(vector).gp(rotor.reverse())));

            Multivector motor = Multivector.scalar(P3.SIG, rnd.nextDouble() * 2 - 1)
                    .add(Multivector.blade(P3.SIG, 3, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 5, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 6, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 9, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 10, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 12, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 15, rnd.nextDouble() * 2 - 1));
            Multivector point = Multivector.blade(P3.SIG, 7, rnd.nextDouble() * 2 - 1)
                    .add(Multivector.blade(P3.SIG, 11, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 13, rnd.nextDouble() * 2 - 1))
                    .add(Multivector.blade(P3.SIG, 14,
                            i % 5 == 0 ? 0 : rnd.nextDouble() * 2 - 1));
            pgaError = Math.max(pgaError,
                    motor.sandwich(point).distanceTo(motor.gp(point).gp(motor.reverse())));
        }
        check("E3 vector sandwich fast path", e3Error < 1e-12);
        check("PGA point sandwich fast path", pgaError < 1e-12);
    }

    private static boolean throwsIllegal(Runnable runnable) {
        try {
            runnable.run();
            return false;
        } catch (IllegalArgumentException e) {
            return true;
        }
    }

    private static void check(String name, boolean ok) {
        if (!ok) {
            throw new AssertionError("FAILED: " + name);
        }
        passed++;
        System.out.println("ok  " + name);
    }
}
