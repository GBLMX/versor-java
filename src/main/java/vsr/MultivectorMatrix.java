package vsr;

/**
 * Dense matrix whose entries are multivectors of one signature.
 *
 * <p>Multiplication is the noncommutative geometric product. {@link #inverse()}
 * uses Gauss-Jordan elimination with left row operations and partial pivoting
 * by {@link Multivector#norm()}; it supports matrices whose pivots have a
 * nonzero scalar norm (the same restriction as {@link Multivector#inverse()}).
 */
public final class MultivectorMatrix {
    public final Signature sig;
    public final int rows;
    public final int cols;
    private final Multivector[][] entries;

    public MultivectorMatrix(Signature sig, int rows, int cols) {
        if (rows < 0 || cols < 0) throw new IllegalArgumentException("negative dimension");
        this.sig = sig;
        this.rows = rows;
        this.cols = cols;
        this.entries = new Multivector[rows][cols];
        Multivector zero = Multivector.zero(sig);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                entries[r][c] = zero;
            }
        }
    }

    public static MultivectorMatrix identity(Signature sig, int n) {
        MultivectorMatrix m = new MultivectorMatrix(sig, n, n);
        Multivector one = Multivector.scalar(sig, 1);
        for (int i = 0; i < n; i++) m.entries[i][i] = one;
        return m;
    }

    public Multivector get(int r, int c) {
        return entries[r][c];
    }

    public void set(int r, int c, Multivector value) {
        checkSignature(value);
        entries[r][c] = value;
    }

    public MultivectorMatrix copy() {
        MultivectorMatrix m = new MultivectorMatrix(sig, rows, cols);
        for (int r = 0; r < rows; r++) {
            System.arraycopy(entries[r], 0, m.entries[r], 0, cols);
        }
        return m;
    }

    public MultivectorMatrix add(MultivectorMatrix o) {
        checkShape(o);
        MultivectorMatrix m = new MultivectorMatrix(sig, rows, cols);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                m.entries[r][c] = entries[r][c].add(o.entries[r][c]);
            }
        }
        return m;
    }

    public MultivectorMatrix sub(MultivectorMatrix o) {
        checkShape(o);
        MultivectorMatrix m = new MultivectorMatrix(sig, rows, cols);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                m.entries[r][c] = entries[r][c].sub(o.entries[r][c]);
            }
        }
        return m;
    }

    public MultivectorMatrix scale(double s) {
        MultivectorMatrix m = new MultivectorMatrix(sig, rows, cols);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                m.entries[r][c] = entries[r][c].scale(s);
            }
        }
        return m;
    }

    /** {@code this * o}; the inner dimension must agree. */
    public MultivectorMatrix multiply(MultivectorMatrix o) {
        if (o.sig != sig || cols != o.rows) {
            throw new IllegalArgumentException("matrix shape mismatch");
        }
        MultivectorMatrix m = new MultivectorMatrix(sig, rows, o.cols);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < o.cols; c++) {
                Multivector sum = Multivector.zero(sig);
                for (int k = 0; k < cols; k++) {
                    sum = sum.add(entries[r][k].gp(o.entries[k][c]));
                }
                m.entries[r][c] = sum;
            }
        }
        return m;
    }

    public MultivectorMatrix transpose() {
        MultivectorMatrix m = new MultivectorMatrix(sig, cols, rows);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                m.entries[c][r] = entries[r][c];
            }
        }
        return m;
    }

    /**
     * Matrix inverse by Gauss-Jordan elimination. Row operations act on the
     * left, so this computes the left inverse; for a square matrix over a
     * division algebra that is also the right inverse. Throws when no pivot
     * with a nonzero scalar norm remains.
     */
    public MultivectorMatrix inverse() {
        if (rows != cols) throw new IllegalStateException("matrix is not square");
        int n = rows;
        Multivector[][] aug = new Multivector[n][2 * n];
        Multivector zero = Multivector.zero(sig);
        Multivector one = Multivector.scalar(sig, 1);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) aug[r][c] = entries[r][c];
            for (int c = n; c < 2 * n; c++) aug[r][c] = zero;
            aug[r][n + r] = one;
        }
        for (int col = 0; col < n; col++) {
            int pivot = -1;
            Multivector inv = null;
            boolean[] tried = new boolean[n];
            for (int attempt = col; attempt < n && pivot < 0; attempt++) {
                int candidate = -1;
                double best = -1;
                for (int r = col; r < n; r++) {
                    if (tried[r]) continue;
                    double score = aug[r][col].maxAbs();
                    if (score > best) {
                        best = score;
                        candidate = r;
                    }
                }
                if (candidate < 0 || best == 0) break;
                tried[candidate] = true;
                try {
                    inv = aug[candidate][col].inverseGeneral();
                    pivot = candidate;
                } catch (ArithmeticException ignored) {
                    // non-invertible candidate; try the next row
                }
            }
            if (pivot < 0) {
                throw new ArithmeticException("matrix is singular");
            }
            Multivector[] swap = aug[col];
            aug[col] = aug[pivot];
            aug[pivot] = swap;
            for (int c = 0; c < 2 * n; c++) {
                aug[col][c] = inv.gp(aug[col][c]);
            }
            for (int r = 0; r < n; r++) {
                if (r == col) continue;
                Multivector factor = aug[r][col];
                if (factor.isZero(0)) continue;
                for (int c = 0; c < 2 * n; c++) {
                    aug[r][c] = aug[r][c].sub(factor.gp(aug[col][c]));
                }
            }
        }
        MultivectorMatrix m = new MultivectorMatrix(sig, n, n);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) m.entries[r][c] = aug[r][n + c];
        }
        return m;
    }

    /** Largest absolute coefficient difference between two matrices. */
    public double distanceTo(MultivectorMatrix o) {
        checkShape(o);
        double max = 0;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                max = Math.max(max, entries[r][c].distanceTo(o.entries[r][c]));
            }
        }
        return max;
    }

    private void checkShape(MultivectorMatrix o) {
        if (o.sig != sig || o.rows != rows || o.cols != cols) {
            throw new IllegalArgumentException("matrix shape mismatch");
        }
    }

    private void checkSignature(Multivector value) {
        if (value.sig != sig) {
            throw new IllegalArgumentException("different algebras");
        }
    }
}
