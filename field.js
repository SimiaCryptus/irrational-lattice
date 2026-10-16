// Generation of the algebraic colored lattice field eta : Z^d -> K^d.
//
// Supports two engines:
// 1. 'exact': Exact algebraic integer Beatty/cut-and-project gratings in
//             Z[sqrt(D)] using BigInt floor(n*sqrt(D)) and exact sign tests
//             for nearest-lattice snapping.
// 2. 'trig':  Smooth trigonometric preview field using incommensurate
//             frequencies derived from sqrt(D).

import { QuadField, floorSqrtDBig, roundQuadBig } from './algebra.js';
import { rationalApprox } from './rational.js';

// Deterministic pseudo-irrational vector built from sqrt(D) and a seed.
export function makeFrequencies(D, K, alphaScale, seed) {
  const sqrtD = Math.sqrt(D);
  const phi = (1 + Math.sqrt(5)) / 2;
  const alphas = [];
  const betas = [];
  const phases = [];
  const ampsR = [];
  const ampsS = [];

  // Integer linear-form coefficients for the exact Z[sqrt(D)] engine.
  // Scale alphaScale into an integer multiplier over a fixed rational grid Q_FREQ.
  const Q_FREQ = 1024;
  const scaleInt = Math.max(1, Math.round(alphaScale * Q_FREQ));
  const formsA = [];
  const formsB = [];
  const phaseInt = [];
  // Common denominator for harmonic weights 1/k: use a fixed multiple so BigInt
  // accumulation stays exact.
  const W_SCALE = 2520n; // lcm(1..10)
  const weightsBig = [];

  for (let k = 1; k <= K; k++) {
    // Build alpha_k, beta_k in [-alphaScale/2, alphaScale/2)^2
    const ax = alphaScale * fractCentered(k * sqrtD + 0.137 * seed);
    const ay = alphaScale * fractCentered(k * phi + 0.241 * seed);
    const bx = alphaScale * fractCentered(k * sqrtD * phi + 0.353 * seed);
    const by = alphaScale * fractCentered((k + 1) * sqrtD - 0.467 * seed);
    alphas.push([ax, ay]);
    betas.push([bx, by]);
    phases.push([
      2 * Math.PI * fractCentered(k * 0.6180339 + 0.1 * seed),
      2 * Math.PI * fractCentered(k * 0.7548776 + 0.2 * seed),
    ]);
    // Spectral weights: 1/k for pink-ish tilt.
    ampsR.push(1 / k);
    ampsS.push(1 / k);

    // Exact integer linear forms l_k(x, y) = c1*x + c2*y + c0 for Z[sqrt(D)] gratings.
    // Derived deterministically from (k, D, seed) without floating-point drift.
    const cAx = Math.round(scaleInt * fractCentered(k * sqrtD + 0.137 * seed)) || (k % 2 ? 1 : -1);
    const cAy =
      Math.round(scaleInt * fractCentered((k + 2) * sqrtD + 0.241 * seed)) || (k % 3 ? -1 : 1);
    const cBx =
      Math.round(scaleInt * fractCentered((2 * k + 1) * sqrtD + 0.353 * seed)) || (k % 2 ? -1 : 1);
    const cBy =
      Math.round(scaleInt * fractCentered((k + 1) * sqrtD - 0.467 * seed)) || (k % 3 ? 1 : -1);
    formsA.push([BigInt(cAx), BigInt(cAy)]);
    formsB.push([BigInt(cBx), BigInt(cBy)]);
    phaseInt.push([BigInt(k * 17 + seed * 31), BigInt(k * 29 + seed * 47)]);
    weightsBig.push(BigInt(Math.max(1, Math.round(Number(W_SCALE) / k))));
  }

  return {
    alphas,
    betas,
    phases,
    ampsR,
    ampsS,
    formsA,
    formsB,
    phaseInt,
    weightsBig,
    Q_FREQ: BigInt(Q_FREQ),
    W_SCALE,
  };
}

function fractCentered(v) {
  const f = v - Math.floor(v);
  return f - 0.5; // in [-0.5, 0.5)
}

// Compute smooth trigonometric eta(x) for one site, returning a 2D vector of
// QuadField elements: [ [a_x, b_x], [a_y, b_y] ].
export function etaAt(x, y, params) {
  const { alphas, betas, phases, ampsR, ampsS } = params;
  let ax = 0,
    ay = 0,
    bx = 0,
    by = 0;
  const K = alphas.length;
  for (let k = 0; k < K; k++) {
    const a = alphas[k];
    const b = betas[k];
    const ph = phases[k];
    const rk = ampsR[k];
    const sk = ampsS[k];

    const angA = 2 * Math.PI * (a[0] * x + a[1] * y) + ph[0];
    const angB = 2 * Math.PI * (b[0] * x + b[1] * y) + ph[1];

    const sA = Math.sin(angA);
    const cB = Math.cos(angB);

    // x-component
    ax += rk * sA;
    bx += sk * cB;
    // y-component, decorrelated via +0.5 rad shift
    ay += rk * Math.cos(angA + 0.5);
    by += sk * Math.sin(angB + 0.5);
  }
  return [
    [ax, bx],
    [ay, by],
  ];
}

/**
 * Exact algebraic Beatty / cut-and-project field over Q(sqrt(D)) at rational
 * coordinates x = xNum/xDen, y = yNum/xDen.
 *
 * For each mode k, let n = l_k(x, y). The centered fractional part
 *   < (n / Q) * sqrt(D) > = (n * sqrt(D) - floor(n * sqrt(D) / Q) * Q) / Q - 1/2
 * lives in Q(sqrt(D)) with rational part (-floor(...) - 1/2) and sqrt(D)
 * coefficient (n / Q).
 *
 * Returns both floating-point pairs [[ax, bx], [ay, by]] and exact BigInt
 * numerators over common denominator `denBig` so snapping can use exact sign tests.
 */
export function etaExactAt(xNum, yNum, coordDen, params, Db) {
  const { formsA, formsB, phaseInt, weightsBig, Q_FREQ, W_SCALE } = params;
  const K = formsA.length;
  const M = Q_FREQ * coordDen; // integer modulus for the linear form
  const halfM = M; // we work over denominator 2 * M * W_SCALE

  let numAx = 0n,
    numBx = 0n;
  let numAy = 0n,
    numBy = 0n;

  for (let k = 0; k < K; k++) {
    const fA = formsA[k];
    const fB = formsB[k];
    const ph = phaseInt[k];
    const w = weightsBig[k];

    // Integer linear combinations for x and y components.
    const nAx = fA[0] * xNum + fA[1] * yNum + ph[0] * coordDen;
    const nBx = fB[0] * xNum + fB[1] * yNum + ph[1] * coordDen;
    const nAy = -fA[1] * xNum + fA[0] * yNum + (ph[0] + 7n) * coordDen;
    const nBy = -fB[1] * xNum + fB[0] * yNum + (ph[1] + 13n) * coordDen;

    // Centered fractional part of (n * sqrt(D)) / M:
    // Let m = floor((n * sqrt(D)) / M). Then:
    //   (n*sqrt(D))/M - m - 1/2 = (-(2*m + 1)*M + 2*n*sqrt(D)) / (2*M)
    const mAx = floorDivSqrtD(nAx, M, Db);
    const mAy = floorDivSqrtD(nAy, M, Db);

    numAx += w * (-(2n * mAx + 1n) * halfM);
    numBx += w * (2n * nAx);

    numAy += w * (-(2n * mAy + 1n) * halfM);
    numBy += w * (2n * nAy);

    // Fold second grating B into the rational/irrational pair via multiplication
    // in Q(sqrt(D)) so both gratings interfere algebraically.
    const mBx = floorDivSqrtD(nBx, M, Db);
    const mBy = floorDivSqrtD(nBy, M, Db);
    numAx += w * (-(2n * mBx + 1n) * halfM);
    numBx += w * (2n * nBx);
    numAy += w * (-(2n * mBy + 1n) * halfM);
    numBy += w * (2n * nBy);
  }

  // Common denominator for (numA + numB * sqrt(D)) / denBig:
  // Factor of 2 from the two gratings summed per component.
  const denBig = 4n * M * W_SCALE;
  const invDen = 1 / Number(denBig);

  const ax = Number(numAx) * invDen;
  const bx = Number(numBx) * invDen;
  const ay = Number(numAy) * invDen;
  const by = Number(numBy) * invDen;

  return {
    pair: [
      [ax, bx],
      [ay, by],
    ],
    numAx,
    numBx,
    numAy,
    numBy,
    denBig,
  };
}

/**
 * Exact floor((n * sqrt(D)) / M) for BigInt n, M > 0, and square-free D > 1.
 * Equivalent to floor(sqrt(n^2 * D / M^2)) with sign handling.
 */
function floorDivSqrtD(n, M, Db) {
  if (n === 0n) return 0n;
  // floor( (n * sqrt(D)) / M ) = floorSqrtDBig(n, D) when M == 1;
  // for general M > 0, we compute floor( floorSqrtDBig(n, D) / M ).
  // Since floor(floor(z) / M) == floor(z / M) for integer M > 0:
  const f = floorSqrtDBig(n, Db);
  return f >= 0n ? f / M : (f - M + 1n) / M;
}

// Compute a scalar field over a grid using the chosen mode.
export function computeField(opts) {
  const { D, K, alphaScale, seed, mode, epsilon } = opts;
  const engine = opts.engine || 'trig';
  const normMode = opts.normMode || 'minmax';
  const width = opts.width || opts.size;
  const height = opts.height || opts.size;
  const panX = opts.panX || 0;
  const panY = opts.panY || 0;
  const zoom = opts.zoom || 1;
  // Upsample: supersample each pixel with a u x u subgrid at stride zoom / u,
  // then box-filter the resulting scalar values (fixes Part I §4 / Part II §7).
  const upsample = opts.upsample && opts.upsample > 1 ? Math.min(8, Math.round(opts.upsample)) : 1;
  const offsetX = opts.offsetX || 0;
  const offsetY = opts.offsetY || 0;
  const cmap2d = opts.cmap2d || 'none';

  const field = new QuadField(D);
  const params = makeFrequencies(D, K, alphaScale, seed);
  const out = new Float32Array(width * height);

  const useChannel2 = cmap2d && cmap2d !== 'none';
  const chan2 = useChannel2 ? new Float32Array(width * height) : null;
  let min2 = Infinity,
    max2 = -Infinity;

  let minV = Infinity,
    maxV = -Infinity;
  let irrSumSq = 0,
    irrCount = 0;

  // Rational coordinate approximation for the exact engine.
  const COORD_DEN = 1024n;
  const coordDenNum = Number(COORD_DEN);
  const epsRat = rationalApprox(epsilon, 1000);
  const epsNumBig = BigInt(epsRat.num);
  const epsDenBig = BigInt(epsRat.den);

  const invU = 1 / upsample;
  const invU2 = 1 / (upsample * upsample);

  for (let j = 0; j < height; j++) {
    for (let i = 0; i < width; i++) {
      let vSum = 0;
      let v2Sum = 0;

      for (let sy = 0; sy < upsample; sy++) {
        const subOffY = upsample > 1 ? (sy + 0.5) * invU - 0.5 : 0;
        const y = (j + subOffY - height / 2) * zoom + panY;
        const yo = y + offsetY;

        for (let sx = 0; sx < upsample; sx++) {
          const subOffX = upsample > 1 ? (sx + 0.5) * invU - 0.5 : 0;
          const x = (i + subOffX - width / 2) * zoom + panX;
          const xo = x + offsetX;

          let e,
            exactInfo = null;
          if (engine === 'exact') {
            const xNum = BigInt(Math.round(xo * coordDenNum));
            const yNum = BigInt(Math.round(yo * coordDenNum));
            exactInfo = etaExactAt(xNum, yNum, COORD_DEN, params, field.Db);
            e = exactInfo.pair;
          } else {
            e = etaAt(xo, yo, params);
          }

          const ax = e[0][0],
            bx = e[0][1];
          const ay = e[1][0],
            by = e[1][1];

          // Real embedding a + b*sqrt(D) and Galois conjugate a - b*sqrt(D).
          const rx = field.toReal(e[0]);
          const ry = field.toReal(e[1]);
          const cx = field.conj(e[0]);
          const cy = field.conj(e[1]);

          const irr2 = bx * bx + by * by;
          irrSumSq += irr2;
          irrCount++;

          // Helper for snapped residual: uses exact sign test in Q(sqrt(D)) when
          // running in exact mode, and float round in trig preview mode.
          const computeSnapResidual = () => {
            if (exactInfo) {
              const { numAx, numBx, numAy, numBy, denBig } = exactInfo;
              const xNum = BigInt(Math.round(xo * coordDenNum));
              const yNum = BigInt(Math.round(yo * coordDenNum));
              // Total coordinate: xNum/COORD_DEN + (epsNum/epsDen) * (numA + numB*sqrt(D))/denBig
              const totDen = COORD_DEN * epsDenBig * denBig;
              const totAx = xNum * epsDenBig * denBig + COORD_DEN * epsNumBig * numAx;
              const totBx = COORD_DEN * epsNumBig * numBx;
              const totAy = yNum * epsDenBig * denBig + COORD_DEN * epsNumBig * numAy;
              const totBy = COORD_DEN * epsNumBig * numBy;

              const snapX = roundQuadBig(totAx, totBx, totDen, field.Db);
              const snapY = roundQuadBig(totAy, totBy, totDen, field.Db);

              const invTot = 1 / Number(totDen);
              const dx =
                Number(totAx - snapX * totDen) * invTot +
                Number(totBx) * invTot * field.sqrtD;
              const dy =
                Number(totAy - snapY * totDen) * invTot +
                Number(totBy) * invTot * field.sqrtD;
              return Math.hypot(dx, dy);
            }
            const tx = xo + epsilon * rx;
            const ty = yo + epsilon * ry;
            return Math.hypot(tx - Math.round(tx), ty - Math.round(ty));
          };

          let v;
          switch (mode) {
            case 'rational':
              v = Math.hypot(ax, ay);
              break;
            case 'irrational':
              v = Math.hypot(bx, by) * field.sqrtD;
              break;
            case 'conjugate':
              v = Math.hypot(cx, cy);
              break;
            case 'snapped':
              v = computeSnapResidual();
              break;
            case 'magnitude':
            default:
              v = Math.hypot(rx, ry);
              break;
          }
          vSum += v;

          if (useChannel2) {
            let v2;
            switch (cmap2d) {
              case 'rational_irrational':
                v2 = Math.hypot(ax, ay);
                break;
              case 'conj_mag':
                v2 = Math.hypot(cx, cy);
                break;
              case 'xy_phase':
                v2 = (Math.atan2(ry, rx) + Math.PI) / (2 * Math.PI);
                break;
              case 'snap':
              default:
                v2 = mode === 'snapped' ? v : computeSnapResidual();
                break;
            }
            v2Sum += v2;
          }
        }
      }

      const vAvg = vSum * invU2;
      out[j * width + i] = vAvg;
      if (vAvg < minV) minV = vAvg;
      if (vAvg > maxV) maxV = vAvg;

      if (useChannel2) {
        const v2Avg = v2Sum * invU2;
        chan2[j * width + i] = v2Avg;
        if (v2Avg < min2) min2 = v2Avg;
        if (v2Avg > max2) max2 = v2Avg;
      }
    }
  }

  // Analytic bounds (Part I §3, Part II §9):
  // H_K = sum_{k=1..K} 1/k.
  let HK = 0;
  for (let k = 1; k <= K; k++) HK += 1 / k;
  let analyticMin = 0;
  let analyticMax;
  switch (mode) {
    case 'rational':
      analyticMax = Math.SQRT2 * HK;
      break;
    case 'irrational':
      analyticMax = Math.SQRT2 * HK * field.sqrtD;
      break;
    case 'snapped':
      analyticMax = Math.SQRT1_2; // sqrt(2)/2
      break;
    case 'conjugate':
    case 'magnitude':
    default:
      analyticMax = Math.SQRT2 * HK * (1 + field.sqrtD);
      break;
  }

  let analyticMin2 = 0,
    analyticMax2 = 1;
  if (useChannel2) {
    switch (cmap2d) {
      case 'rational_irrational':
        analyticMax2 = Math.SQRT2 * HK;
        break;
      case 'conj_mag':
        analyticMax2 = Math.SQRT2 * HK * (1 + field.sqrtD);
        break;
      case 'xy_phase':
        analyticMax2 = 1;
        break;
      case 'snap':
      default:
        analyticMax2 = Math.SQRT1_2;
        break;
    }
  }

  // Nyquist ratio: alphaScale * (zoom / upsample). Values >= 1 indicate spatial aliasing.
  const nyquistRatio = (alphaScale * zoom) / upsample;
  const irrRMS = Math.sqrt(irrSumSq / irrCount);

  return {
    data: out,
    min: normMode === 'analytic' ? analyticMin : minV,
    max: normMode === 'analytic' ? analyticMax : maxV,
    rawMin: minV,
    rawMax: maxV,
    analyticMin,
    analyticMax,
    irrRMS,
    nyquistRatio,
    chan2,
    min2: normMode === 'analytic' ? analyticMin2 : min2,
    max2: normMode === 'analytic' ? analyticMax2 : max2,
    width,
    height,
  };
}