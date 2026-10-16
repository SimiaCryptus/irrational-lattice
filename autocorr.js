// Autocorrelation analysis for the lattice field.
//
// We compute the 2D autocorrelation of the (downsampled) field via FFT,
// then locate the strongest off-center peaks. The displacement vectors to
// those peaks (in lattice units) make natural "step" directions for a
// random walk that tends to land the viewport on self-similar features.

import { computeFFT2D } from './fft.js';

// Radix-2 in-place FFT (re/im arrays). Mirrors the one in fft.js but kept
// local so this module is self-contained for inverse transforms.
function fft1d(re, im, inverse) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      const tr = re[i];
      re[i] = re[j];
      re[j] = tr;
      const ti = im[i];
      im[i] = im[j];
      im[j] = ti;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = ((inverse ? 2 : -2) * Math.PI) / len;
    const wpr = Math.cos(ang);
    const wpi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wr = 1,
        wi = 0;
      for (let k = 0; k < len / 2; k++) {
        const a = i + k;
        const b = i + k + len / 2;
        const tr = wr * re[b] - wi * im[b];
        const ti = wr * im[b] + wi * re[b];
        re[b] = re[a] - tr;
        im[b] = im[a] - ti;
        re[a] += tr;
        im[a] += ti;
        const nwr = wr * wpr - wi * wpi;
        wi = wr * wpi + wi * wpr;
        wr = nwr;
      }
    }
  }
}

function largestPow2LE(n) {
  let p = 1;
  while (p * 2 <= n) p *= 2;
  return p;
}

function fft2d(re, im, N, inverse) {
  const row = new Float64Array(N);
  const rowI = new Float64Array(N);
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      row[x] = re[y * N + x];
      rowI[x] = im[y * N + x];
    }
    fft1d(row, rowI, inverse);
    for (let x = 0; x < N; x++) {
      re[y * N + x] = row[x];
      im[y * N + x] = rowI[x];
    }
  }
  const col = new Float64Array(N);
  const colI = new Float64Array(N);
  for (let x = 0; x < N; x++) {
    for (let y = 0; y < N; y++) {
      col[y] = re[y * N + x];
      colI[y] = im[y * N + x];
    }
    fft1d(col, colI, inverse);
    for (let y = 0; y < N; y++) {
      re[y * N + x] = col[y];
      im[y * N + x] = colI[y];
    }
  }
}

// Compute the top-2 autocorrelation displacement vectors (in lattice units).
// `data` is a Float32Array (size*size). `zoom` converts pixel lags into
// lattice units. Returns an array of up to 2 {dx, dy, strength} entries.
export function topAutocorrVectors(data, width, zoom, maxN = 64, height = width) {
  console.group('[autocorr] topAutocorrVectors');
  console.log('inputs:', {
    dataLength: data ? data.length : null,
    width,
    height,
    zoom,
    maxN,
  });
  if (!data || data.length === 0) {
    console.warn('[autocorr] empty/missing data array');
    console.groupEnd();
    return [];
  }
  const minDim = Math.min(width, height);
  const N = largestPow2LE(Math.min(minDim, maxN));
  console.log('[autocorr] downsample grid N =', N, '(expected N*N =', N * N, ')');
  if (N < 2) {
    console.warn('[autocorr] N < 2, grid too small to analyze');
    console.groupEnd();
    return [];
  }
   // Zero-pad to M = 2*N so Wiener-Khinchin computes linear (non-circular)
   // autocorrelation rather than wrapping opposite edges (Part II §8).
   const M = N * 2;
   const samples = new Float64Array(N * N);
   const re = new Float64Array(M * M);
   const im = new Float64Array(M * M);

   // Block-average downsample and remove mean (so DC doesn't dominate).
  const scaleX = width / N;
  const scaleY = height / N;
  let mean = 0;
  let dsMin = Infinity,
    dsMax = -Infinity;
  let nonFinite = 0;
  for (let y = 0; y < N; y++) {
     const y0 = Math.floor(y * scaleY);
     const y1 = Math.max(y0 + 1, Math.min(height, Math.floor((y + 1) * scaleY)));
    for (let x = 0; x < N; x++) {
       const x0 = Math.floor(x * scaleX);
       const x1 = Math.max(x0 + 1, Math.min(width, Math.floor((x + 1) * scaleX)));
       let sum = 0,
         cnt = 0;
       for (let sy = y0; sy < y1; sy++) {
         const rowOff = sy * width;
         for (let sx = x0; sx < x1; sx++) {
           sum += data[rowOff + sx];
           cnt++;
         }
       }
       const v = cnt > 0 ? sum / cnt : 0;
      if (!isFinite(v)) nonFinite++;
      if (v < dsMin) dsMin = v;
      if (v > dsMax) dsMax = v;
       samples[y * N + x] = v;
      mean += v;
    }
  }
  mean /= N * N;
  console.log('[autocorr] downsample stats:', {
    scaleX,
    scaleY,
    mean,
    dsMin,
    dsMax,
    range: dsMax - dsMin,
    nonFinite,
  });
  if (dsMax - dsMin < 1e-12) {
    console.warn(
      '[autocorr] field is (nearly) constant after downsampling — autocorrelation will be flat and no off-center peaks will be found.'
    );
  }
  if (nonFinite > 0) {
    console.warn('[autocorr]', nonFinite, 'non-finite samples in downsampled field');
  }
   // Place zero-mean samples into the top-left N x N quadrant of the M x M buffer.
   for (let y = 0; y < N; y++) {
     for (let x = 0; x < N; x++) {
       re[y * M + x] = samples[y * N + x] - mean;
     }
   }

   // Linear autocorrelation = IFFT( |FFT(zero-padded x)|^2 ).
   fft2d(re, im, M, false);
   for (let i = 0; i < M * M; i++) {
    const p = re[i] * re[i] + im[i] * im[i];
    re[i] = p;
    im[i] = 0;
  }
   fft2d(re, im, M, true);
   // Normalize inverse transform by M*M and sample count N*N.
   const norm = M * M * (N * N);
   for (let i = 0; i < M * M; i++) re[i] /= norm;
  // Diagnostic: report the central (zero-lag) autocorrelation value, which
  // equals the field variance and should be the global maximum.
  {
    let acMin = Infinity,
      acMax = -Infinity,
      acNonFinite = 0;
     for (let i = 0; i < M * M; i++) {
      const v = re[i];
      if (!isFinite(v)) acNonFinite++;
      if (v < acMin) acMin = v;
      if (v > acMax) acMax = v;
    }
    console.log('[autocorr] raw AC (pre-shift) stats:', {
      acMin,
      acMax,
      acNonFinite,
      zeroLagValue: re[0],
    });
  }

   // Extract central N x N lags [-N/2, N/2) from the M x M linear autocorrelation.
  const half = N / 2;
  const ac = new Float64Array(N * N);
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
       const lagX = x - half;
       const lagY = y - half;
       const mx = (lagX + M) % M;
       const my = (lagY + M) % M;
       ac[y * N + x] = re[my * M + mx];
    }
  }

  // Find local maxima away from the central (zero-lag) peak.
  const center = half;
  const minR = 2; // ignore tiny lags near the center
  const candidates = [];
  let belowMinR = 0;
  for (let y = 1; y < N - 1; y++) {
    for (let x = 1; x < N - 1; x++) {
      const dx = x - center;
      const dy = y - center;
      if (Math.hypot(dx, dy) < minR) {
        belowMinR++;
        continue;
      }
      const v = ac[y * N + x];
      // Local maximum test (3x3 neighborhood).
      let isMax = true;
      for (let oy = -1; oy <= 1 && isMax; oy++) {
        for (let ox = -1; ox <= 1; ox++) {
          if (ox === 0 && oy === 0) continue;
          if (ac[(y + oy) * N + (x + ox)] > v) {
            isMax = false;
            break;
          }
        }
      }
      if (isMax && v > 0) candidates.push({ x, y, v, dx, dy });
    }
  }
  candidates.sort((a, b) => b.v - a.v);
  console.log('[autocorr] candidate peaks:', {
    count: candidates.length,
    skippedNearCenter: belowMinR,
    top5: candidates.slice(0, 5).map((c) => ({
      dx: c.dx,
      dy: c.dy,
      v: c.v,
    })),
  });
  if (candidates.length === 0) {
    console.warn(
      '[autocorr] no off-center local maxima with v > 0 found. ' +
        'Field may be too smooth/constant or autocorrelation degenerate.'
    );
  }

  // Pick the top 2 that aren't near-duplicates / mirror images.
  const picked = [];
  // Dedup distance scales with the grid so larger lags aren't all merged.
  const dedupR = Math.max(minR, N * 0.08);
  for (const c of candidates) {
    let dup = false;
    for (const p of picked) {
      if (
        Math.hypot(c.dx - p.rdx, c.dy - p.rdy) < dedupR ||
        Math.hypot(c.dx + p.rdx, c.dy + p.rdy) < dedupR
      ) {
        dup = true;
        break;
      }
    }
    if (dup) continue;
    picked.push({ rdx: c.dx, rdy: c.dy, strength: c.v });
    if (picked.length >= 2) break;
  }
  console.log('[autocorr] picked (after dedup):', {
    dedupR,
    count: picked.length,
    picked: picked.map((p) => ({ rdx: p.rdx, rdy: p.rdy, strength: p.strength })),
  });

  // Convert pixel lags into lattice units (each pixel lag in the downsampled
  // grid corresponds to `scale` source pixels, each `zoom` lattice units).
  // Use the axis-specific scales so non-square fields convert correctly.
  const result = picked.map((p) => ({
    dx: p.rdx * scaleX * zoom,
    dy: p.rdy * scaleY * zoom,
    strength: p.strength,
  }));
  console.log('[autocorr] final vectors (lattice units):', result);
  console.groupEnd();
  return result;
}