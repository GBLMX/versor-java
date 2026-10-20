/**
 * Geometric algebra for graphics and kinematics.
 *
 * <h2>Layers</h2>
 * <ul>
 *   <li>{@link vsr.Signature}: metric, blade encoding (bit mask), grade and
 *       product-sign rules. All operands must share one instance.</li>
 *   <li>{@link vsr.Multivector}: dense coefficient array indexed by mask,
 *       with an adaptive sparse-support path for the geometric product.</li>
 *   <li>{@link vsr.SparseMultivector}: COO (sorted masks + values) for high
 *       dimensions or low density; picks a hash or dense accumulator per
 *       operation.</li>
 *   <li>{@link vsr.MultivectorBatch}: structure-of-arrays batch with optional
 *       output reuse and slicing for parallel blocks.</li>
 *   <li>{@link vsr.MultivectorMatrix}: matrix of multivectors with
 *       Gauss-Jordan inversion.</li>
 *   <li>{@link vsr.E3}, {@link vsr.P3}, {@link vsr.Conformal},
 *       {@link vsr.Motion}, {@link vsr.Field}, {@link vsr.Calculus},
 *       {@link vsr.Manifold}: geometry and kinematics on top of the core.</li>
 * </ul>
 *
 * <p>Blades are ordered by their integer mask; a product of basis blades is
 * either zero (null metric) or a signed basis blade. The sign table is built
 * lazily for n &le; 10; larger algebras compute signs on demand.
 */
package vsr;
