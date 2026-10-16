/**
 * fractals.js
 *
 * Six domain-specific fractal generators built around Q(sqrt(D)),
 * its Galois conjugate embedding (a - b*sqrt(D)), exact BigInt quadratic
 * integer arithmetic, continued-fraction Gauss dynamics, and iterated
 * lattice snapping.
 *
 * Each mode returns { values, values2, width, height, stats } matching
 * computeField(opts) so renderField(canvas, result) works unmodified.
 */

import { makeFrequencies, etaAt } from './field.js';

const PHI = (1 + Math.sqrt(5)) / 2;
const PLASTIC = 1.324717957244746;

function fractCentered(t) {
  return (t - Math.floor(t)) - 0.5;
}

/**
 * Exact integer floor(sqrt(N)) for BigInt N >= 0 via Newton's method (Part II §6).
 */
function isqrtBig(N) {
  if (N < 0n) throw new RangeError('Negative BigInt in isqrtBig');
  if (N < 2n) return N;
  let x0 = 1n << BigInt((N.toString(2).length + 1) >> 1);
  let x1 = (x0 + N / x0) >> 1n;
  while (x1 < x0) {
    x0 = x1;
    x1 = (x0 + N / x0) >> 1n;
  }
  return x0;
}

/**
 * Exact floor(n * sqrt(D)) for signed integer n and positive square-free D.
 * Uses floor(n * sqrt(D)) = sgn(n)*isqrt(n^2 * D) - [n < 0] (Part II §6).
 */
function exactFloorNSqrtD(nInt, DInt) {
  const n = BigInt(Math.trunc(nInt));
  const D = BigInt(Math.trunc(DInt));
  if (n === 0n) return 0n;
  const mag = isqrtBig(n * n * D);
  return n > 0n ? mag : -mag - 1n;
}

/**
 * Exact sign test in Q(sqrt(D)) for (a + b*sqrt(D)) where a, b are BigInt numerators
 * over a common denominator (Part II §4).
 * Returns -1, 0, or +1 with zero floating-point error.
 */
function exactSignQuad(a, b, D) {
  const sA = a > 0n ? 1 : a < 0n ? -1 : 0;
  const sB = b > 0n ? 1 : b < 0n ? -1 : 0;
  if (sB === 0) return sA;
  if (sA === 0) return sB;
  if (sA === sB) return sA;
  const diff = a * a - D * b * b;
  const sDiff = diff > 0n ? 1 : diff < 0n ? -1 : 0;
  return sA * sDiff;
}

/**
 * Evaluate both the real embedding r = a + b*sqrt(D), the Galois conjugate
 * r_bar = a - b*sqrt(D), and the quadratic field norm N(xi) = a^2 - D*b^2.
 */
function evalAlgebraicJet(x1, x2, freqs, sqrtD, D) {
   const e = etaAt(x1, x2, freqs);
   const ax = e[0][0];
   const bx = e[0][1];
   const ay = e[1][0];
   const by = e[1][1];
   const rx = ax + bx * sqrtD;
   const ry = ay + by * sqrtD;
  // Galois conjugate embedding sigma_-(a + b*sqrt(D)) = a - b*sqrt(D)
  const rxConj = ax - bx * sqrtD;
  const ryConj = ay - by * sqrtD;
  // Field norm in K = Q(sqrt(D)): N(a + b*sqrt(D)) = a^2 - D*b^2
  const normX = ax * ax - D * bx * bx;
  const normY = ay * ay - D * by * by;
  return { ax, bx, ay, by, rx, ry, rxConj, ryConj, normX, normY };
}

/**
 * 1. Iterated Snap Basin Boundaries ("Lattice Julia Set")
 *
 * Iterates z_{n+1} = z_n + eps * eta(z_n) - gamma * (z_n - round(z_n)).
 * Primary channel: smoothed locking time + log residual.
 * Secondary channel: Golden-ratio spatial hash of the terminal lattice basin (m, n) in Z^2.
 */
function pixelSnapJulia(x1, x2, freqs, sqrtD, D, opts) {
  const maxIter = opts.fractalIter || 18;
   const eps = opts.epsilon ?? opts.eps ?? 0.35;
  const gamma = opts.snapGravity ?? 0.42;
  const lockTol2 = 1e-5;

  let zx = x1;
  let zy = x2;
  let smoothIter = maxIter;
  let normAccum = 0;

  for (let i = 0; i < maxIter; i++) {
    const jet = evalAlgebraicJet(zx, zy, freqs, sqrtD, D);
    const nx = zx + eps * jet.rx;
    const ny = zy + eps * jet.ry;

    const sx = Math.round(nx);
    const sy = Math.round(ny);
    const dx = nx - sx;
    const dy = ny - sy;
    const dist2 = dx * dx + dy * dy;

    normAccum += Math.abs(jet.normX + jet.normY);

    // Pull toward nearest lattice site while perturbed by incommensurate field
    zx = nx - gamma * dx;
    zy = ny - dy * gamma;

    if (dist2 < lockTol2) {
      const logDist = Math.log(Math.max(dist2, 1e-16)) * 0.5;
      smoothIter = i + 1 - Math.log(Math.abs(logDist)) / Math.LN2;
      break;
    }
  }

  // Final lattice site index (m, n) in Z^2
  const m = Math.round(zx);
  const n = Math.round(zy);
  // Quasi-random algebraic hash of the terminal lattice basin
  const basinHash = fractCentered(m * sqrtD + n * PHI + 0.5) + 0.5;
  const seamResidual = Math.hypot(zx - m, zy - n);

  const v1 = smoothIter + 0.15 * Math.log1p(normAccum) - 2.0 * seamResidual;
  const v2 = basinHash;
  return [v1, v2];
}

/**
 * 2. Galois Conjugate Hofstadter Carpet (Minkowski Internal-Space IFS)
 *
 * Uses the Minkowski embedding (xi, xi_bar) -> R^2.
 * Folds coordinates by the real embedding while shifting by the Galois conjugate
 * and coloring by the hyperbolic norm form N(a + b*sqrt(D)) = a^2 - D*b^2.
 */
function pixelGaloisCarpet(x1, x2, freqs, sqrtD, D, opts) {
  const maxIter = opts.fractalIter || 14;
   const eps = opts.epsilon ?? opts.eps ?? 0.4;
  const windowRadius = opts.cutWindow || 1.25;

  let zx = x1;
  let zy = x2;
  let minNormHyp = Infinity;
  let windowHits = 0;
  let orbitTrap = 0;

  // Fundamental inflation factor tied to sqrt(D)
  const scale = Math.min(sqrtD, 2.6180339887);

  for (let i = 0; i < maxIter; i++) {
    const jet = evalAlgebraicJet(zx, zy, freqs, sqrtD, D);

    // Hyperbolic norm distance to the light-cone a^2 - D*b^2 = 0
    const hyper = Math.abs(jet.normX * jet.normY);
    if (hyper < minNormHyp) minNormHyp = hyper;

    // Cut-and-project acceptance test in internal (conjugate) space
    const conjMag = Math.hypot(jet.rxConj, jet.ryConj);
    if (conjMag < windowRadius) {
      windowHits += 1;
    }

    // Distance to half-integer snap seams in physical space
    const fx = fractCentered(zx + eps * jet.rx);
    const fy = fractCentered(zy + eps * jet.ry);
    const seamDist = Math.min(Math.abs(fx), Math.abs(fy));
    orbitTrap += Math.exp(-18.0 * seamDist) / (i + 1);

    // Minkowski fold: inflate fractional part by sqrt(D) and twist by conjugate
    zx = scale * fx + 0.35 * jet.rxConj;
    zy = scale * fy + 0.35 * jet.ryConj;
  }

  const v1 = Math.log1p(orbitTrap * 4.0) - 0.25 * Math.log(minNormHyp + 1e-6);
  const v2 = windowHits / maxIter;
  return [v1, v2];
}

/**
 * 3. Metallic Weierstrass-Beatty Cascade (Exact Z[sqrt(D)] Sturmian Fractal)
 *
 * Replaces trigonometric functions with a multi-scale Weierstrass cascade of
 * exact algebraic fractional parts {n * sqrt(D)} computed via BigInt isqrt.
 */
function pixelBeattyCascade(x1, x2, freqs, sqrtD, D, opts) {
  const octaves = Math.min(opts.fractalIter || 10, 12);
  const hurst = opts.hurst ?? 0.55;
  const DInt = Math.max(2, Math.round(D));

  // Precompute golden rotation cos/sin
  const cosA = Math.cos(2 * Math.PI / PHI);
  const sinA = Math.sin(2 * Math.PI / PHI);

  let ux = x1 * (opts.alphaScale || 0.08) * 16.0;
  let uy = x2 * (opts.alphaScale || 0.08) * 16.0;
  let amp = 1.0;
  const lacunarity = sqrtD;
  const gain = Math.pow(lacunarity, -hurst);

  let sumReal = 0;
  let sumConj = 0;

  for (let m = 0; m < octaves; m++) {
    // Discretize scaled coordinate onto integer lattice for exact Beatty step
    const n1 = Math.round(ux * 32);
    const n2 = Math.round(uy * 32);

    // Exact floor(n * sqrt(D)) via BigInt Newton isqrt:
    // {n * sqrt(D)} = n * sqrt(D) - floor(n * sqrt(D)) in Z + Z*sqrt(D)
    const floor1 = Number(exactFloorNSqrtD(n1, DInt));
    const floor2 = Number(exactFloorNSqrtD(n2, DInt));

    // Exact algebraic integer pair (a, b) = (-floor(n*sqrt(D)), n)
    const frac1 = (n1 * sqrtD - floor1) - 0.5;
    const frac2 = (n2 * sqrtD - floor2) - 0.5;

    // Galois conjugate of the exact fractional part: a - b*sqrt(D) normalized
    const conjWave = Math.sin(2 * Math.PI * (floor1 * PHI - floor2 / PLASTIC));

    // Sturmian interference ridge
    const ridge = Math.abs(Math.abs(frac1) - Math.abs(frac2));
    sumReal += amp * (0.5 - ridge);
    sumConj += amp * conjWave * ridge;

    // Rotate by golden angle and inflate by sqrt(D)
    const rx = (ux * cosA - uy * sinA) * lacunarity;
    const ry = (ux * sinA + uy * cosA) * lacunarity;
    ux = rx;
    uy = ry;
    amp *= gain;
  }

  return [sumReal, 0.5 + 0.5 * Math.tanh(sumConj)];
}

/**
 * 4. Diophantine Resonance Flame (Continued-Fraction Gauss Map & Arnold Tongues)
 *
 * Applies the continued-fraction Gauss map t_{n+1} = {1 / t_n} to the local
 * ratio of rational-to-irrational field components, measuring how badly
 * approximable the field is at each point (Part II §10d).
 */
function pixelDiophantineFlame(x1, x2, freqs, sqrtD, D, opts) {
  const maxSteps = opts.fractalIter || 14;
  const jet = evalAlgebraicJet(x1, x2, freqs, sqrtD, D);

  // Ratio of rational to irrational coordinates
  let tx = Math.abs((jet.ax + 0.01 * x1) / (jet.bx * sqrtD + 1e-7));
  let ty = Math.abs((jet.ay + 0.01 * x2) / (jet.by * sqrtD + 1e-7));
  tx = tx - Math.floor(tx);
  ty = ty - Math.floor(ty);

  // Convergent denominators k_n = a_n * k_{n-1} + k_{n-2}
  let kx0 = 0, kx1 = 1;
  let ky0 = 0, ky1 = 1;
  let brjunoSum = 0;
  let nobleScore = 0;
  let weight = 0.5;

  for (let n = 0; n < maxSteps; n++) {
    if (tx < 1e-9 || ty < 1e-9) break;

    const invX = 1.0 / tx;
    const invY = 1.0 / ty;
    const ax = Math.floor(invX);
    const ay = Math.floor(invY);

    tx = invX - ax;
    ty = invY - ay;

    const kx2 = ax * kx1 + kx0;
    const ky2 = ay * ky1 + ky0;
    kx0 = kx1; kx1 = Math.min(kx2, 1e12);
    ky0 = ky1; ky1 = Math.min(ky2, 1e12);

    // Brjuno-like diophantine sum: explodes near small-denominator rationals (Arnold tongues)
    brjunoSum += weight * Math.log1p(ax + ay);

    // Noble / quadratic irrational detector: bounded small partial quotients (a_n <= D)
    if (ax <= D && ay <= D) {
      nobleScore += 1;
    }
    weight *= 0.72;
  }

  const logDenom = Math.log1p(Math.hypot(kx1, ky1));
  return [brjunoSum, nobleScore / maxSteps + 0.05 * logDenom];
}

/**
 * 5. Precision-Horizon Catastrophe Map (Float64 vs. Exact Z[sqrt(D)] Snap Divergence)
 *
 * Directly visualizes Part II §3 & §4: sweeps simulated pan horizons 2^b and
 * compares IEEE-754 float64 rounding against the exact BigInt sign test in Q(sqrt(D)).
 */
function pixelHorizonCatastrophe(x1, x2, freqs, sqrtD, D, opts) {
   const eps = opts.epsilon ?? opts.eps ?? 0.45;
  const DInt = BigInt(Math.max(2, Math.round(D)));
  const jet = evalAlgebraicJet(x1, x2, freqs, sqrtD, D);

  // Quantize (ax, bx) onto high-denominator rationals to represent an exact element of Q(sqrt(D))
  const Q = 1048576n; // 2^20 denominator
  const QNum = 1048576.0;

  const axBig = BigInt(Math.round((x1 + eps * jet.ax) * QNum));
  const bxBig = BigInt(Math.round(eps * jet.bx * QNum));
  const ayBig = BigInt(Math.round((x2 + eps * jet.ay) * QNum));
  const byBig = BigInt(Math.round(eps * jet.by * QNum));

  // Exact integer nearest-lattice candidate via float, then verified by exactSignQuad
  // Test sign of (a + b*sqrt(D)) - (m + 0.5)
  const mCand = Math.floor(x1 + eps * jet.rx);
  const nCand = Math.floor(x2 + eps * jet.ry);

  const halfQ = Q >> 1n;
  const aPrimeX = axBig - BigInt(mCand) * Q - halfQ;
  const aPrimeY = ayBig - BigInt(nCand) * Q - halfQ;

  const exactSnapX = mCand + (exactSignQuad(aPrimeX, bxBig, DInt) >= 0 ? 1 : 0);
  const exactSnapY = nCand + (exactSignQuad(aPrimeY, byBig, DInt) >= 0 ? 1 : 0);

  // Now sweep floating-point pan horizons P = 2^b for b = 1..52 and detect when
  // IEEE-754 loses enough mantissa bits to flip the snap decision
  let firstFlipBit = 53;
  let disagreementCount = 0;

  for (let b = 8; b <= 52; b += 2) {
    const panHorizon = Math.pow(2, b);
    // Simulate evaluating phase & snap at extreme pan (x + panHorizon) - panHorizon
    const shiftedX = (x1 + panHorizon) + eps * Math.sin(2 * Math.PI * ((x1 + panHorizon) * sqrtD));
    const floatSnapX = Math.round(shiftedX - panHorizon + eps * (jet.rx - Math.sin(2 * Math.PI * x1 * sqrtD)));

    const shiftedY = (x2 + panHorizon) + eps * Math.cos(2 * Math.PI * ((x2 + panHorizon) * PHI));
    const floatSnapY = Math.round(shiftedY - panHorizon + eps * (jet.ry - Math.cos(2 * Math.PI * x2 * PHI)));

    if (floatSnapX !== exactSnapX || floatSnapY !== exactSnapY) {
      disagreementCount++;
      if (firstFlipBit === 53) firstFlipBit = b;
    }
  }

  // Distance to the undecidable half-integer boundary a' + b*sqrt(D) = 0
  const seamDistX = Math.abs((x1 + eps * jet.rx) - (mCand + 0.5));
  const seamDistY = Math.abs((x2 + eps * jet.ry) - (nCand + 0.5));
  const seamProximity = -Math.log(Math.min(seamDistX, seamDistY) + 1e-7);

  const v1 = (53 - firstFlipBit) + 0.35 * seamProximity;
  const v2 = disagreementCount / 23.0;
  return [v1, v2];
}

/**
 * 6. Strange Attractor Density of the Quasicrystal Kick Map
 *
 * Whole-viewport histogram renderer: couples physical embedding r = a + b*sqrt(D)
 * with conjugate embedding r_bar = a - b*sqrt(D) in a area-preserving kick map.
 */
function computeQuasicrystalAttractor(opts, freqs, sqrtD, D) {
   const { width: W, height: H, zoom = 1, panX = 0, panY = 0, offsetX = 0, offsetY = 0 } = opts;
  const values = new Float32Array(W * H);
  const values2 = new Float32Array(W * H);

   const step = zoom;
  const spanX = W * step;
  const spanY = H * step;
  const cx = panX + offsetX;
  const cy = panY + offsetY;

   const eps = opts.epsilon ?? opts.eps ?? 0.35;
  const numSeeds = 64;
  const stepsPerSeed = Math.max(4000, Math.floor((W * H * 6) / numSeeds));

  let maxDensity = 1;
   let irrSumSq = 0;
   let irrCount = 0;

  for (let s = 0; s < numSeeds; s++) {
    // Seed points distributed by plastic/golden low-discrepancy sequence
    let x = cx + (fractCentered(s * PHI) * spanX);
    let y = cy + (fractCentered(s * PLASTIC) * spanY);

    for (let i = 0; i < stepsPerSeed; i++) {
      const jet = evalAlgebraicJet(x, y, freqs, sqrtD, D);
       irrSumSq += jet.bx * jet.bx + jet.by * jet.by;
       irrCount++;

      // Symplectic-style twist coupling real and Galois conjugate embeddings
      const nx = x + eps * jet.rx + 0.15 * eps * jet.ryConj;
      const jetNext = evalAlgebraicJet(nx, y, freqs, sqrtD, D);
      const ny = y + eps * jetNext.ryConj - 0.15 * eps * jetNext.rx;

      // Toroidal wrap back into current viewport
      x = cx + fractCentered((nx - cx) / spanX + 0.5) * spanX;
      y = cy + fractCentered((ny - cy) / spanY + 0.5) * spanY;

      const px = Math.floor((x - cx) / step + W * 0.5);
      const py = Math.floor((y - cy) / step + H * 0.5);

      if (px >= 0 && px < W && py >= 0 && py < H) {
        const idx = py * W + px;
        const c = values[idx] + 1;
        values[idx] = c;
        if (c > maxDensity) maxDensity = c;

        // Second channel tracks local algebraic norm sign (hyperbolic sector)
        const normSign = Math.atan(jetNext.normX + jetNext.normY) / Math.PI + 0.5;
        values2[idx] = 0.9 * values2[idx] + 0.1 * normSign;
      }
    }
  }

  // Log-density tone mapping
  const logMax = Math.log1p(maxDensity);
   let minV = Infinity,
     maxV = -Infinity;
   let min2 = Infinity,
     max2 = -Infinity;
  for (let i = 0; i < values.length; i++) {
     const v = Math.log1p(values[i]) / (logMax || 1);
     values[i] = v;
     if (v < minV) minV = v;
     if (v > maxV) maxV = v;
     const v2 = values2[i];
     if (v2 < min2) min2 = v2;
     if (v2 > max2) max2 = v2;
  }

  return {
     data: values,
     chan2: values2,
    values,
    values2,
     min: minV,
     max: maxV,
     rawMin: minV,
     rawMax: maxV,
     min2,
     max2,
    width: W,
    height: H,
     irrRMS: irrCount > 0 ? Math.sqrt(irrSumSq / irrCount) : 0,
    nyquistRatio: (opts.alphaScale || 0.08) * step
  };
}

/**
 * Main dispatcher for all fractal modes.
 *
 * Supported opts.fractalMode values:
 *   - 'snap_julia'          : Iterated Snap Basin Boundaries
 *   - 'galois_carpet'       : Galois Conjugate Hofstadter Carpet
 *   - 'beatty_weierstrass'  : Exact Z[sqrt(D)] Metallic Weierstrass Cascade
 *   - 'diophantine_flame'   : Continued-Fraction Gauss Map & Arnold Tongues
 *   - 'horizon_catastrophe' : Float64 vs Exact BigInt Q(sqrt(D)) Snap Divergence
 *   - 'quasicrystal_attractor' : Strange Attractor Density of Conjugate Kick Map
 */
export function computeFractalField(opts) {
  const {
    D = 5,
    K = 12,
    alphaScale = 0.08,
    seed = 1,
    zoom = 1,
    panX = 0,
    panY = 0,
    offsetX = 0,
    offsetY = 0,
    fractalMode = 'snap_julia'
  } = opts;
   const W = opts.width || opts.size || 256;
   const H = opts.height || opts.size || 256;
   const upsample = opts.upsample && opts.upsample > 1 ? Math.min(8, Math.round(opts.upsample)) : 1;

  const sqrtD = Math.sqrt(D);
   const freqs = makeFrequencies(D, K, alphaScale, seed);

  if (fractalMode === 'quasicrystal_attractor') {
     return computeQuasicrystalAttractor({ ...opts, width: W, height: H }, freqs, sqrtD, D);
  }

  const values = new Float32Array(W * H);
  const values2 = new Float32Array(W * H);

   const step = zoom;
  const halfW = W * 0.5;
  const halfH = H * 0.5;
  const ox = panX + offsetX;
  const oy = panY + offsetY;

  let kernel;
  switch (fractalMode) {
    case 'galois_carpet':
      kernel = pixelGaloisCarpet;
      break;
    case 'beatty_weierstrass':
      kernel = pixelBeattyCascade;
      break;
    case 'diophantine_flame':
      kernel = pixelDiophantineFlame;
      break;
    case 'horizon_catastrophe':
      kernel = pixelHorizonCatastrophe;
      break;
    case 'snap_julia':
    default:
      kernel = pixelSnapJulia;
      break;
  }
   let minV = Infinity,
     maxV = -Infinity;
   let min2 = Infinity,
     max2 = -Infinity;
   const invU = 1 / upsample;
   const invU2 = 1 / (upsample * upsample);


  for (let j = 0; j < H; j++) {
    const rowOffset = j * W;
    for (let i = 0; i < W; i++) {
       let v1Sum = 0;
       let v2Sum = 0;
       for (let sy = 0; sy < upsample; sy++) {
         const subOffY = upsample > 1 ? (sy + 0.5) * invU - 0.5 : 0;
         const x2 = (j + subOffY - halfH) * step + oy;
         for (let sx = 0; sx < upsample; sx++) {
           const subOffX = upsample > 1 ? (sx + 0.5) * invU - 0.5 : 0;
           const x1 = (i + subOffX - halfW) * step + ox;
           const [s1, s2] = kernel(x1, x2, freqs, sqrtD, D, opts);
           v1Sum += isFinite(s1) ? s1 : 0;
           v2Sum += isFinite(s2) ? s2 : 0;
         }
       }
       const v1 = v1Sum * invU2;
       const v2 = v2Sum * invU2;
      values[rowOffset + i] = v1;
      values2[rowOffset + i] = v2;
       if (v1 < minV) minV = v1;
       if (v1 > maxV) maxV = v1;
       if (v2 < min2) min2 = v2;
       if (v2 > max2) max2 = v2;
    }
  }

  return {
     data: values,
     chan2: values2,
    values,
    values2,
     min: minV,
     max: maxV,
     rawMin: minV,
     rawMax: maxV,
     min2,
     max2,
    width: W,
    height: H,
    irrRMS: 0,
     nyquistRatio: (alphaScale * step) / upsample
  };
}