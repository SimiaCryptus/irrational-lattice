// Quadratic field arithmetic over Q(sqrt(D)).
// Supports both floating-point embeddings for display and exact BigInt
// arithmetic in Z[sqrt(D)] / Q(sqrt(D)) for exact fractional parts,
// comparisons, and nearest-lattice snapping.

/**
 * Exact floor(sqrt(n)) for a non-negative BigInt n via Newton's method.
 */
export function isqrtBig(n) {
  n = BigInt(n);
  if (n < 0n) throw new RangeError('isqrtBig of negative number');
  if (n < 2n) return n;
  // Seed from bit length.
  const bits = n.toString(2).length;
  let x = 1n << BigInt((bits + 1) >> 1);
  while (true) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
}

/**
 * Exact floor(n * sqrt(D)) for integer n (BigInt or Number) and square-free
 * positive integer D (BigInt or Number).
 * Since D is square-free > 1, n*sqrt(D) is an integer iff n == 0.
 */
export function floorSqrtDBig(n, D) {
  const nb = BigInt(n);
  const Db = BigInt(D);
  if (nb === 0n) return 0n;
  const mag = isqrtBig(nb * nb * Db);
  return nb > 0n ? mag : -mag - 1n;
}

/**
 * Exact sign of (numA + numB * sqrt(D)) / den, where numA, numB, D, den are
 * BigInts and den > 0. Returns -1, 0, or 1.
 *
 * Uses the algebraic identity:
 *   sign(a' + b*sqrt(D)) = sign(a') * sign(a'^2 - D * b^2)
 * when sign(a') != sign(b).
 */
export function signQuadBig(numA, numB, D) {
  const a = BigInt(numA);
  const b = BigInt(numB);
  const Db = BigInt(D);
  if (a === 0n && b === 0n) return 0;
  if (b === 0n) return a > 0n ? 1 : -1;
  if (a === 0n) return b > 0n ? 1 : -1;
  const sa = a > 0n ? 1 : -1;
  const sb = b > 0n ? 1 : -1;
  if (sa === sb) return sa;
  const diff = a * a - Db * b * b;
  if (diff === 0n) return 0;
  return diff > 0n ? sa : -sa;
}

/**
 * Exact nearest integer (round half to +infinity) of (numA + numB*sqrt(D)) / den
 * for BigInt numA, numB, den > 0, and square-free D > 1.
 * Returns a BigInt m such that (numA + numB*sqrt(D))/den - m in [-1/2, 1/2).
 */
export function roundQuadBig(numA, numB, den, D) {
  const a = BigInt(numA);
  const b = BigInt(numB);
  const q = BigInt(den);
  const Db = BigInt(D);
  // Floor of b*sqrt(D):
  const f = floorSqrtDBig(b, Db);
  // Estimate integer part of (a + f) / q (Euclidean floor division).
  const sum = a + f;
  let m = sum >= 0n ? sum / q : (sum - q + 1n) / q;
  // Refine m using exact comparison against (m + 1/2) and (m - 1/2):
  // We want 2*(a + b*sqrt(D)) - (2*m + 1)*q < 0  and  2*(a + b*sqrt(D)) - (2*m - 1)*q >= 0.
  while (signQuadBig(2n * a - (2n * m + 1n) * q, 2n * b, Db) >= 0) {
    m += 1n;
  }
  while (signQuadBig(2n * a - (2n * m - 1n) * q, 2n * b, Db) < 0) {
    m -= 1n;
  }
  return m;
}

export class QuadField {
  constructor(D) {
    this.D = D;
    this.Db = BigInt(D);
    this.sqrtD = Math.sqrt(D);
  }

  // Construct element from rational/irrational parts.
  el(a, b) {
    return [a, b];
  }

  add(x, y) {
    return [x[0] + y[0], x[1] + y[1]];
  }
  sub(x, y) {
    return [x[0] - y[0], x[1] - y[1]];
  }
  scale(x, s) {
    return [x[0] * s, x[1] * s];
  }

  // Exact multiplication in Q(sqrt(D)): (a1 + b1*sqrt(D))*(a2 + b2*sqrt(D)).
  mul(x, y) {
    return [x[0] * y[0] + this.D * x[1] * y[1], x[0] * y[1] + x[1] * y[0]];
  }

  // Algebraic norm: a^2 - D*b^2.
  norm(x) {
    return x[0] * x[0] - this.D * x[1] * x[1];
  }

  // Real-valued embedding: a + b*sqrt(D).
  toReal(x) {
    return x[0] + x[1] * this.sqrtD;
  }

  // Galois conjugate embedding: a - b*sqrt(D) (internal-space coordinate).
  conj(x) {
    return x[0] - x[1] * this.sqrtD;
  }
}