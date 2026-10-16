# Algebraic Colored Lattice Fields: Algorithm and Commentary

This document has two parts.

- **Part I** describes the algorithm as implemented in this repository (`field.js`, `render.js`, `fft.js`, `autocorr.js`, `audio.js`, `rational.js`, `main.js`). It records the exact formulas, the coordinate conventions, and the auxiliary analyses (spectrum, autocorrelation, sonification).
- **Part II** discusses the numerical methods behind it and argues that the project is, read carefully, a **critique of the real numbers**. It looks at what the code claims, what a floating-point machine can actually do, and where exact algebraic arithmetic could replace the continuum.

Notation: $\mathbb{Z}^2$ is the integer lattice, $D$ a square-free positive integer, $K = \mathbb{Q}(\sqrt D)$ the real quadratic field, $\varphi = (1+\sqrt5)/2$ the golden ratio, $\{t\}$ the fractional part, and $\langle t\rangle = \{t\} - \tfrac12 \in [-\tfrac12, \tfrac12)$ the centered fractional part (`fractCentered`).

---

# Part I — The Algorithm

## 1. Overview

The pipeline, run once per redraw:

```
UI state ─► readOpts() ─► computeField(opts) ─► renderField(canvas, result)
                                │
                                ├─► computeFFT2D  ─► renderFFT3D        (spectrum window)
                                ├─► topAutocorrVectors ─► random walk   (on demand)
                                └─► FieldAudio.updateField              (if playing)
```

`computeField` evaluates a vector-valued displacement field
$$\eta : \mathbb{R}^2 \to (\mathbb{R}\times\mathbb{R})^2,\qquad \eta(x) = \big(a_x(x) + b_x(x)\sqrt D,\; a_y(x) + b_y(x)\sqrt D\big)$$
at one sample point per output pixel. It then reduces each sample to a scalar according to the selected mode and, optionally, to a second scalar for 2D colormaps.

## 2. Frequency construction (`makeFrequencies`)

For $k = 1, \dots, K$ and integer seed $s$, with $\alpha_{\text{scale}}$ the "frequency scale" slider:

$$
\alpha_k = \alpha_{\text{scale}}\Big(\langle k\sqrt D + 0.137s\rangle,\ \langle k\varphi + 0.241s\rangle\Big)
$$
$$
\beta_k = \alpha_{\text{scale}}\Big(\langle k\sqrt D\,\varphi + 0.353s\rangle,\ \langle (k+1)\sqrt D - 0.467s\rangle\Big)
$$

Phases:
$$
\phi_k = 2\pi\Big(\langle 0.6180339\,k + 0.1s\rangle,\ \langle 0.7548776\,k + 0.2s\rangle\Big)
$$
The constants are $1/\varphi$ and $1/\rho$, where $\rho \approx 1.3247$ is the plastic number. Both are badly approximable, so $\{k c\}$ is spread evenly over $[0,1)$ (Weyl equidistribution, three-gap theorem).

Amplitudes: $r_k = s_k = 1/k$. This gives a $1/k$ spectral tilt, a "pink-ish" weighting in mode index. It is not a pink weighting in spatial frequency, because $|\alpha_k|$ is not monotone in $k$.

Properties:

- Every coordinate of $\alpha_k,\beta_k$ lies in $[-\alpha_{\text{scale}}/2,\ \alpha_{\text{scale}}/2)$.
- In exact arithmetic, the coordinates are irrational for generic $s$. The "incommensurability" comes from the sequences $\{k\sqrt D\}$ and $\{k\varphi\}$.
- The mixed term $k\sqrt D\varphi$ lies in $\mathbb{Q}(\sqrt D, \sqrt 5)$. This is a biquadratic field when $D \ne 5$, so the construction already leaves $\mathbb{Q}(\sqrt D)$ in the frequencies.

## 3. The displacement field (`etaAt`)

For a sample point $x = (x_1, x_2) \in \mathbb{R}^2$, define
$$\theta^A_k = 2\pi\,\alpha_k\cdot x + \phi_{k,1},\qquad \theta^B_k = 2\pi\,\beta_k\cdot x + \phi_{k,2}.$$

The four components are:
$$
\begin{aligned}
a_x &= \textstyle\sum_k r_k \sin\theta^A_k, & b_x &= \textstyle\sum_k s_k \cos\theta^B_k,\\
a_y &= \textstyle\sum_k r_k \cos(\theta^A_k + 0.5), & b_y &= \textstyle\sum_k s_k \sin(\theta^B_k + 0.5).
\end{aligned}
$$

The pair $(a_\bullet, b_\bullet)$ is read as the "rational" and "irrational" coordinates of an element $a + b\sqrt D$. Its real embedding is
$$r_x = a_x + b_x\sqrt D,\qquad r_y = a_y + b_y\sqrt D.$$

The $+0.5$ rad shift decorrelates the $y$-component from the $x$-component. A shift of $\pi/2$ would make $(a_x, a_y)$ a pure rotation of a single phasor; $0.5$ rad deliberately avoids that.

Bounds: $|a_\bullet|, |b_\bullet| \le H_K = \sum_{k\le K} 1/k \approx \ln K + 0.5772$. So $\|r\| \le \sqrt2\,H_K(1+\sqrt D)$.

The `K_field` argument to `etaAt` is unused. The `QuadField` class (`algebra.js`) is only used for its cached `sqrtD`; `conj` (the Galois conjugate embedding $a - b\sqrt D$) is never called. See Part II §9.

## 4. Sampling geometry (`computeField` + `main.js`)

Pixel $(i, j)$ in a render of size $W\times H$ maps to the lattice point
$$
x_1 = (i - W/2)\,z_{\text{eff}}\,u + p_x + o_x,\qquad
x_2 = (j - H/2)\,z_{\text{eff}}\,u + p_y + o_y
$$
where:

- $p = (p_x, p_y)$ is the continuous pan, set by dragging, wheel-zoom anchoring, or the autocorrelation walk.
- $o = (o_x, o_y) \in \mathbb{Z}^2$ is the integer offset (buttons).
- $u$ is the "upsample" slider value.
- $z_{\text{eff}} = z \cdot 512 / \min(W,H)$ is the effective zoom. $z$ is stored in "lattice units per reference pixel" with reference size 512. Changing grid size therefore changes resolution but not the visible region.

The render dimensions follow the display aspect. The grid-size control fixes the smaller of $W,H$, and the other is scaled so pixels stay square (`fitCanvas`).

Note on `upsample`: the comment in `field.js` says it samples at $1/u$-unit steps, but the code multiplies the step by $u$. As written, raising $u$ widens the sampling stride: it zooms out and subsamples, rather than refining. Part II §6 discusses the aliasing consequences.

Wheel zoom multiplies $z$ by the rational granularity $p/q$ (or its inverse) from `rational.js`. It then re-solves for the pan so that the lattice point under the cursor stays fixed:
$$p \leftarrow x_{\text{before}} - (f - c)\,z_{\text{eff}}^{\text{new}}.$$

Drag panning rounds the pixel delta to integers, so the pan moves in whole multiples of $z_{\text{eff}}$.

## 5. Scalar modes

| mode | value $v(x)$ |
|---|---|
| `magnitude` | $\sqrt{r_x^2 + r_y^2}$ |
| `rational` | $\sqrt{a_x^2 + a_y^2}$ |
| `irrational` | $\sqrt D\sqrt{b_x^2 + b_y^2}$ |
| `snapped` | $\operatorname{dist}\big(T_\varepsilon(x), \mathbb{Z}^2\big)$ coordinatewise, see below |

**Snap map.** With $\varepsilon$ the amplitude slider:
$$T_\varepsilon(x) = x + \varepsilon\,(r_x, r_y),\qquad S(x) = \operatorname{round}(T_\varepsilon(x)),\qquad v = \|T_\varepsilon(x) - S(x)\|_2.$$

So the displayed quantity is the residual of the nearest-lattice quantizer, not the snapped point itself. $v \in [0, \sqrt2/2]$. It is zero exactly when the deformed point lands on a lattice site.

At integer sample points $x\in\mathbb{Z}^2$ (zoom 1, integer pan), this reduces to $\|\varepsilon r - \operatorname{round}(\varepsilon r)\|$. This is the distance of the displacement itself to the lattice, which is the "nudge then snap back" picture from the README.

**Diagnostic.** `irrRMS` $= \sqrt{\operatorname{mean}(b_x^2 + b_y^2)}$ over the viewport is reported as a rigidity diagnostic. It is measured, not enforced; no lower bound on $|b|$ is imposed anywhere.

## 6. Second channel and colormaps

If a 2D colormap is selected, a second scalar $v_2$ is computed:

- `snap`: the snap residual above (independent of mode).
- `rational_irrational`: $\|a\|$.
- `xy_phase`: $(\operatorname{atan2}(r_y, r_x) + \pi)/2\pi$.

`render.js` normalizes both channels by min–max over the current viewport:
$$t = \frac{v - \min v}{\max v - \min v}.$$
The 2D maps send $(u,v) = (t_2, t_1)$ to HSV: hue from $u$ plus a cycling phase, value from $v$. The 1D maps are piecewise-linear interpolations of perceptual stops (viridis, magma, plasma, …) or HSV cycles (rainbow, disco, neon).

Color cycling only advances the phase and re-renders the cached field without recomputing it.

## 7. Spectrum (`fft.js`)

1. **Downsample.** Choose $N$ = the largest power of two $\le \min(\min(W,H), 64)$. Fill the $N\times N$ grid by nearest-pixel point sampling with independent strides $W/N$ and $H/N$. The comment says "block average", but the code takes a single nearest sample.
2. **Transform.** Apply a separable 2D FFT: rows then columns, iterative radix-2 Cooley–Tukey with bit-reversal and twiddles updated by recurrence.
3. **Log magnitude.** Compute $\log(1 + |F|)$ and fftshift so that DC is centered.
4. **Display.** Draw as an oblique-projected 3D height surface. Quads are depth-sorted by the projected screen-y of each row (painter's algorithm) and colored by a 5-stop viridis-like ramp.

The mean is not removed, so the DC bin is usually the tallest peak.

## 8. Autocorrelation and the self-similar walk (`autocorr.js`)

1. Downsample as in §7. Subtract the mean.
2. Compute the circular autocorrelation by Wiener–Khinchin:
   $$R = \mathcal{F}^{-1}\big(|\mathcal{F}(f - \bar f)|^2\big)/N^2.$$
   $R(0)$ is $N^2$ times the sample variance (the code divides once by $N^2$ after an unnormalized inverse).
3. fftshift, then find strict-ish local maxima over the 3×3 neighborhood. Exclude lags with $|\tau| < 2$ and keep only $R > 0$.
4. Sort candidates by $R$. Greedily keep up to two, rejecting any within $\max(2, 0.08N)$ of an already-picked $\tau$ or of $-\tau$. Since $R(\tau) = R(-\tau)$, antipodes are duplicates.
5. Convert lags to lattice units: $\Delta = (\tau_x\cdot W/N,\ \tau_y\cdot H/N)\cdot z_{\text{eff}}$.

The **walk** (`acWalkStep`) adds $\pm v_i$ to the pan, each with probability ½ independently. This is a lazy random walk on the 2D lattice generated by the two vectors. Steps run at the chosen rate from the animation loop.

The vectors can also be edited by dragging their endpoint handles. The drag inverts the pixel↔lattice map.

## 9. Sonification (`audio.js`)

**Scan mode.** The bank has 48 sine oscillators at $110\cdot 2^{5i/48}$ Hz, i.e. 5 octaves in equal steps of 1/48 of an octave each. A column cursor advances at $c$ columns/s, wrapping at $W$. For column $j$, oscillator $i$ reads row $\lfloor iH/48\rfloor$. Its gain is set (with a 20 ms smoothing constant) to $t^2/\sqrt{48}$, where $t$ is the min–max-normalized field value. The bank is spectral additive synthesis, so the field column is heard as a spectrum envelope.

**Sample mode.** The field is flattened row-major into a mono buffer, mapped to $[-1,1]$, and looped at the context sample rate. Each row becomes $W$ samples. Row-to-row continuity therefore produces a pitched tone near $f_s/W$, with the field's structure heard as timbre and slow modulation.

## 10. Rational view controls (`rational.js`)

- **Zoom granularity** is an exact rational $p/q$. The slider seeds it through the best continued-fraction convergent with $q \le 100$:
  - expand $x = [a_0; a_1, a_2, \dots]$;
  - build convergents $h_n = a_n h_{n-1} + h_{n-2}$ and $k_n = a_n k_{n-1} + k_{n-2}$;
  - stop before $k_n > q_{\max}$.

  Only full convergents are considered, not semiconvergents. The result is therefore the best approximation of the second kind, not necessarily the closest fraction with $q \le 100$.
- **Grid size** is an integer.

The intent is that the view parameters are exact rationals, even though everything downstream is floating point.

## 11. State persistence

All controls, the pan, the zoom and the offset are serialized into the URL hash. The pan is rounded to 3 decimals and the zoom to 5. A shared link therefore reproduces the view only up to that rounding. This is a small, concrete instance of the theme of Part II: a "point" is only ever shared as a finite rational approximation.

## 12. Complexity

| stage | cost |
|---|---|
| Field | $O(WHK)$ trig evaluations: 4 trig calls per mode per pixel |
| FFT and autocorrelation | $O(N^2\log N)$ with $N \le 64$; negligible |
| Rendering | $O(WH)$ |

At 1024² and $K = 42$ the field costs about $1.8\times10^8$ trig calls, which explains the visible latency at high settings.

---

# Part II — Commentary: Numerical Methods, and a Critique of the Real Numbers

## 1. The thesis

The README promises a perturbation that is "deterministic, algebraically grounded, provably aperiodic … computed in exact integer arithmetic". The code computes it with IEEE-754 doubles, `Math.sin`, and `Math.round`. That gap is not a bug to be patched quietly. It is the most interesting thing about the project.

Every property that makes the construction mathematically special exists in the idealized continuum $\mathbb{R}$. None of those properties survive in what a machine can hold:

- irrationality,
- incommensurability,
- aperiodicity,
- the exact decision "which lattice point is nearest".

The project ends up as a working demonstration of where the real numbers stop paying their way, and of how algebraic number fields can step in.

## 2. A floating-point irrational is a rational

`Math.sqrt(5)` returns a dyadic rational $m/2^{51}$, accurate to about $2^{-53}$ relative error. Every frequency coordinate in `makeFrequencies` is therefore a dyadic rational. A sum of finitely many sinusoids with dyadic frequencies is **periodic**: if every $\alpha_{k,i}$ and $\beta_{k,i}$ is a multiple of $2^{-M}$, then $\eta(x + 2^M e_i) = \eta(x)$ in exact arithmetic on those rational frequencies.

So the computed field is periodic, with a period on the order of $2^{50}$ lattice units or more. The aperiodicity theorem is true of the *idealized* field. What the program shows is a periodic object whose period lies beyond any viewport you can pan to.

This is not pedantry. It shows that "aperiodic" is never verifiable on a finite window. Every finite patch of any field is also a patch of some periodic field. What the mathematics actually guarantees is quantitative:

- **Discrepancy bounds.** How evenly $\{k\sqrt D\}$ is spread, which is controlled by the continued fraction of $\sqrt D$. That expansion is eventually periodic (Lagrange), with bounded partial quotients, so $\sqrt D$ is badly approximable.
- **Near-period penalties.** How far any near-period $\tau$ must be before $\|\eta(\cdot+\tau)-\eta\|$ becomes small.

A finite machine can witness quantitative statements like these. "Aperiodic" in the real-number sense is an infinitary statement it can only gesture at.

## 3. Where precision actually runs out

The phase $\theta = 2\pi\,\alpha\cdot x + \phi$ is evaluated in double precision. Its absolute error is roughly $|\theta|\cdot 2^{-53}$. In addition, `Math.sin` must reduce $\theta$ modulo $2\pi$, and modern engines do this accurately (Payne–Hanek), so the binding error is in forming $\theta$, not in the reduction.

For pans near $10^{12}$ with $\alpha\sim 1$, the phase error is about $6\times10^{12}\cdot 10^{-16}\approx 10^{-3}$ rad, which is still visually invisible. Near $10^{15}$ it reaches order 1, and the field becomes numerical noise wearing an algebraic costume.

Three smaller precision losses:

- `fractCentered(k*sqrtD + 0.137*seed)` loses about $\log_2 k$ bits of the fractional part at mode $k$. This is harmless for $K\le 42$.
- `Float32Array` storage truncates field values to 24 bits. This is irrelevant for display but would matter for any downstream algebra.
- Pan is stored and shared rounded to $10^{-3}$, so the hash-reproduced view differs from the original by up to $5\times10^{-4}$ lattice units. That is enough to flip snap decisions; see §4.

The README's claim of "no floating-point drift" is true in a weaker sense: the computation is **deterministic** (same inputs, same bits, on the same engine). It is not exact. Determinism is a property of the program; exactness would be a property of the numbers.

## 4. The snap is the crux: rounding is not computable on $\mathbb{R}$

The decisive step is $S(x) = \operatorname{round}(x + \varepsilon\eta(x))$. In computable analysis (Turing, Brouwer, Weihrauch), a real number is given by a stream of ever-better rational approximations. A function on such reals is computable only if it is continuous. `round` is discontinuous at every half-integer.

As a result, no algorithm can decide `round(t)` for an arbitrary computable real $t$. If $t$ is exactly $n + \tfrac12$, or merely indistinguishable from it at every finite precision you have looked at so far, the algorithm can never commit.

In floating point the question is answered anyway, because $t$ was silently replaced by a nearby dyadic, and that dyadic's round is decided. The answer is correct for the dyadic. Whether it is correct for the "true" $t$ is unknowable from inside the computation.

Visually, this appears as the bright seams in `snapped` mode: the loci where the residual jumps from $\approx\tfrac{\sqrt2}{2}$ back toward 0. Those seams are exactly the set where the real-number model and the machine model may disagree. They are measure-zero in the continuum, but the pixels sitting on them are real.

**This is where algebraic numbers rescue the construction.** Suppose $t = a + b\sqrt D$ with $a, b \in \mathbb{Q}$. The sign of $t - c$ for rational $c$ is decidable exactly:
$$
\operatorname{sign}(a' + b\sqrt D)\ \text{with}\ a' = a - c:\quad
\begin{cases}
\operatorname{sign}(a') & \text{if } \operatorname{sign}(a') = \operatorname{sign}(b)\ \text{or}\ b=0,\\
\operatorname{sign}(a')\cdot\operatorname{sign}(a'^2 - Db^2) & \text{otherwise.}
\end{cases}
$$
Moreover $a' + b\sqrt D = 0$ only if $a' = b = 0$, because $\sqrt D \notin \mathbb{Q}$. A tie at a half-integer therefore happens only when $b = 0$ and $a$ is a half-integer. That condition is itself decidable, and it is ruled out entirely if $b \ne 0$.

This is the real content of the README's "rigidity" remark: **in $K = \mathbb{Q}(\sqrt D)$, nearest-lattice snapping is a total, exact, decidable function.** In $\mathbb{R}$ it is not.

## 5. The trigonometric leak: why the current $\eta$ is not in $\mathbb{Q}(\sqrt D)$

The pairs $(a_x, b_x)$ produced by `etaAt` are sums of $\sin$ and $\cos$ of real arguments. These values are generally transcendental (Lindemann–Weierstrass), not rational. Writing them as $a + b\sqrt D$ is therefore a labeling, not a membership.

Over $\mathbb{R}$ the decomposition is not even unique: for any real $t$, the map $(a,b)\mapsto(a + t\sqrt D,\ b - t)$ gives the same real value. The "rational part" and "irrational part" are meaningful only because the program *declares* which sum is which. Consequences:

- The `rational` and `irrational` modes display two independent trigonometric fields. Neither is the rational or irrational part of an algebraic number in any intrinsic sense.
- The algebraic structure that is genuinely present lives in the **frequencies**, as integer multiples of $\sqrt D$ and $\varphi$ reduced mod 1. Aperiodicity of the ideal field follows from the frequencies being $\mathbb{Q}$-linearly independent of $\mathbb{Z}^2$ (Kronecker/Weyl). It does not depend on a lower bound on $|b|$.
- No lower bound on $|b|$ is enforced. `irrRMS` reports an average. Pointwise $b_x(x) = b_y(x) = 0$ is entirely possible, since a sum of cosines has zeros.

The README's rigidity argument should therefore be read as a design goal, not a property of `etaAt`.

## 6. An exactly computable alternative

The trig leak can be removed by building $\eta$ from fractional parts of **algebraic integers** instead of sinusoids. Consider $\{m + n\sqrt D\}$ for integers $m, n$:
$$
\lfloor n\sqrt D\rfloor = \operatorname{sgn}(n)\cdot\operatorname{isqrt}(n^2 D) - [\,n<0\,]
$$
Here `isqrt` is the exact integer square root, computed with `BigInt` Newton iteration. The correction term applies for $n<0$, valid because $n\sqrt D$ is never an integer. So
$$
\{n\sqrt D\} = n\sqrt D - \lfloor n\sqrt D\rfloor \in \mathbb{Z} + \mathbb{Z}\sqrt D
$$
is held **exactly** as the integer pair $(-\lfloor n\sqrt D\rfloor,\ n)$.

A deformation built this way stays inside $\mathbb{Z}[\sqrt D]$. For example, take
$$\eta(x) = \sum_k c_k\big(\{\ell_k(x)\sqrt D\} - \tfrac12\big)$$
with integer linear forms $\ell_k$ and rational weights $c_k$. Every snap decision is then made by the exact sign test of §4. The field is the classical **cut-and-project / Beatty** construction. It is genuinely aperiodic (Sturmian in each direction), computed in integer arithmetic, and only converted to floats at the last moment for display.

The README's original promises — integer-arithmetic-clean, provably aperiodic, rigidity from $\sqrt D$ — would then be literally true. The trade-off is spectral smoothness: sawtooth gratings have harmonics. But harmonics can be shaped by the integer weights, which keeps the "spectrally tunable" claim.

The Galois conjugate `conj` in `algebra.js` becomes useful here. For $\xi = a + b\sqrt D$, the pair $(\xi, \bar\xi)$ is the standard embedding of $\mathbb{Z}[\sqrt D]$ as a lattice in $\mathbb{R}^2$. Bounding $|\bar\xi|$ (the "window" in internal space) is exactly how quasicrystal theory proves both uniform discreteness and aperiodicity.

## 7. Sampling, aliasing, and the moiré

A pleasant but uncomfortable observation: some of the "Penrose-like moiré" is **aliasing**. The field contains spatial frequencies up to $|\alpha|_\infty \le \alpha_{\text{scale}}/2$ cycles per lattice unit. Sampled at stride $z_{\text{eff}}\,u$, the Nyquist limit is
$$\alpha_{\text{scale}}\,z_{\text{eff}}\,u < 1 .$$

With the default $\alpha = 0.08$, $z = 1$, a 256 grid and $u = 1$, the product is $0.08\cdot2 = 0.16$, safely sampled. Push $\alpha$ toward 2, zoom out, or raise `upsample` (which, as noted in Part I §4, *coarsens* the stride), and the frequencies fold back into low-frequency beat patterns. These are moiré in the original optical sense: the interference of the field's gratings with the **pixel lattice**, not with each other.

This has a philosophical edge. Aliasing is the discrete lattice $\mathbb{Z}^2$ reasserting itself against the continuum. A real-valued frequency is always seen only through its residue modulo the sampling lattice, which is the same reduction $\langle\cdot\rangle$ that `makeFrequencies` applies on purpose. The viewer and the generator perform the same operation, and the image cannot tell them apart.

Snapped mode is worse. It is discontinuous, so its spectrum decays slowly and it aliases at *every* zoom. That cannot be avoided by point sampling; it needs area sampling (supersampling or analytic coverage of the snap cells).

Two concrete remedies:

1. Make `upsample` mean what its comment says: use stride $z_{\text{eff}}/u$ and box-filter $u\times u$ subsamples per pixel. This turns it into supersampling.
2. Display the Nyquist ratio $\alpha_{\text{scale}} z_{\text{eff}} u$ in the stats panel, so users can tell structure from artifact.

## 8. The spectral and correlation estimators

The README says the ideal field has a **pure-point** spectrum: finitely many Bragg peaks at $\pm\alpha_k$ and $\pm\beta_k$. Nothing a finite FFT computes is pure-point.

- **Leakage.** A rectangular window of side $L$ convolves each delta with a 2D Dirichlet kernel (a periodic sinc). Peaks are smeared over $\sim 1/L$, with sidelobes decaying like $1/f$. A Hann or Kaiser window would trade resolution for much lower sidelobes, and would make the peaks of the ideal spectrum more legible.
- **Pre-FFT aliasing.** The nearest-neighbor downsample to $N\le 64$ adds aliasing of its own: a frequency $f$ appears at $f \bmod (N/L)$. Averaging blocks, as the comment intends, would act as a crude low-pass filter.
- **DC dominance.** The mean is not removed in `computeFFT2D`, so the DC bin is the tallest peak. The log scale hides most of the effect.
- **Circular autocorrelation.** Wiener–Khinchin via FFT computes the *circular* autocorrelation. Lags near $\pm N/2$ mix opposite edges of the window. Zero-padding to $2N$ gives the linear autocorrelation, so the peaks the walk selects would correspond to true near-periods rather than wrap-around coincidences.
- **Lag quantization.** Recovered lags are integers on the $N$ grid, i.e. multiples of $(W/N)z_{\text{eff}}$ lattice units. An irrational near-period gets rounded to that grid. The walk therefore moves along a rational approximation of an irrational near-symmetry: a continued-fraction convergent in disguise.

That last point closes a loop with `rational.js`. The interface represents zoom with continued-fraction convergents on purpose, and the autocorrelation walk represents near-periods with crude rational convergents by accident. In both cases the machine can only hold a rational shadow of an irrational quantity. The interesting question is how good the shadow is.

## 9. Normalization is viewport-relative

Both `render.js` and `audio.js` use min–max normalization over the current window. The color and loudness of a site therefore depend on what else happens to be visible. Two consequences:

- Panning changes colors of *unchanged* sites. This is visible during the autocorrelation walk.
- Comparisons across $D$, $K$, or $\alpha$ are not meaningful pixel-to-pixel.

A fixed normalization by the analytic bound gives an intrinsic scale: $\|r\|\le\sqrt2 H_K(1+\sqrt D)$ for magnitude, or $[0, \sqrt2/2]$ for snapped mode. It is less contrasty, but honest.

## 10. A critique of the real numbers

The pieces above add up to an argument.

**(a) The reals are an idealization that the construction needs but cannot use.** The theorems about this field — aperiodicity, pure-point spectrum, equidistribution of the gratings — are statements about $\mathbb{R}$. The computation never touches $\mathbb{R}$. It uses a finite set of dyadic rationals, plus operations (`sin`, `round`) that the IEEE standard and JS engines define on that set, with rounding rules, ties, and overflow behavior that have no counterpart in the continuum.

**(b) Every interesting operation here is discontinuous, and discontinuous functions are exactly what $\mathbb{R}$ cannot compute.** Snapping, fractional parts, $\operatorname{atan2}$'s branch cut (used by `xy_phase`), and colormap thresholds are all discontinuous. A constructive or computable reading of real analysis says these functions do not exist as total computations. The program uses them freely because it never works with reals. The pictures are made from operations that the real-number model says are uncomputable.

**(c) Equality is undecidable in $\mathbb{R}$ and decidable in $K$.** You cannot, in general, decide whether a computable real is zero. In $\mathbb{Q}(\sqrt D)$ you can: $a + b\sqrt D = 0 \iff a = b = 0$. The README's rigidity principle — "if the irrational part cancels everywhere it must be zero" — is exactly this decidability. It is a property of the number field, not of the real line, and the real line erases it. Once $a + b\sqrt D$ is evaluated to a double, the pair $(a,b)$ cannot be recovered. That is precisely the information the aperiodicity argument needs.

**(d) "Irrational" is relational, not intrinsic, at finite precision.** On a machine, the right notion is not "is $\alpha$ irrational?" but "how well is $\alpha$ approximated by rationals with denominators up to $Q$?" That is a continued-fraction question with a finite answer. Quadratic irrationals are the best possible case: their partial quotients are bounded and eventually periodic, so they stay far from every rational with small denominator. **The reason to prefer $\sqrt D$ is not that it is irrational, but that it is quantitatively, uniformly, and certifiably far from rationals.** That is something a finite computation can actually verify.

**(e) Therefore the honest foundations for this project are algebraic, not real.** The objects with operational meaning are:

- the ring $\mathbb{Z}[\sqrt D]$, which supports exact arithmetic, exact comparison, and exact rounding;
- its Minkowski embedding into $\mathbb{R}^2$ via $(\xi, \bar\xi)$, which turns aperiodicity into a statement about a lattice and a window;
- rational view parameters, as `rational.js` already insists.

The real numbers enter only at the last step, as the target of the display map: a lossy conversion to pixel colors, where losing information is acceptable because the eye has finite precision too.

Put simply, the construction is **a field over $\mathbb{Z}[\sqrt D]$ that happens to be drawn on a screen**. It is not a field over $\mathbb{R}$ that happens to involve $\sqrt D$. The current implementation follows the second description, and every discrepancy catalogued above traces back to that choice.

## 11. Recommendations

1. Make the algebraic core exact:
   - Represent $\eta$ in $\mathbb{Z}[\sqrt D]$ using `BigInt` pairs, built from fractional parts of algebraic integers (§6).
   - Use the exact sign test for snapping (§4).
   - Keep the trig field as an optional "smooth preview" mode, clearly labelled.
2. Use `QuadField.conj`. Display the internal-space coordinate $\bar\xi$ as a mode. This is the canonical certificate of aperiodicity and shows the "window".
3. Fix `upsample` to supersample (§7). Show the Nyquist ratio.
4. Window and zero-pad the FFT and the autocorrelation (§8).
5. Offer analytic normalization in addition to min–max (§9).
6. Encode pan as an exact rational, or as an element of $\mathbb{Q}(\sqrt D)$, in the URL hash, so shared links reproduce snap decisions bit-for-bit.
7. Reword the README:
   - "deterministic" instead of "no floating-point drift";
   - "aperiodic in the ideal limit; the computed field has period $\gtrsim 2^{50}$" instead of an unqualified "provably aperiodic".

## 12. Claims versus implementation

| claim (README) | ideal mathematics | current code |
|---|---|---|
| coordinates in $\mathbb{Q}(\sqrt D)$ | yes, by design | pairs of transcendental doubles labeled $(a,b)$ |
| integer-arithmetic clean | achievable via $\mathbb{Z}[\sqrt D]$ | IEEE-754 doubles and `Math.sin` |
| provably aperiodic | yes (Kronecker/Weyl) | periodic with astronomically large period |
| rigidity from $\|b\|$ bounded below | needs an explicit construction | not enforced; `irrRMS` only observes |
| exact snapping | decidable in $K$ | float `Math.round`, ties resolved arbitrarily |
| pure-point spectrum | yes | leakage-smeared, aliased FFT estimate |
| deterministic and reproducible | yes | yes (per engine), modulo hash rounding |
| spectrally tunable | yes | yes, via $1/k$ weights; extendable |

The gaps in this table are not failures of the idea. Each row marks a place where the real numbers made a promise that only an algebraic number field can keep.
---
# Part III — Algebraic & Quasicrystalline Fractal Plots
The structures in Parts I and II—the snap residual, the unused Galois conjugate $\bar\xi = a - b\sqrt D$, the continued-fraction convergents, and the discontinuous rounding seams—naturally generate fractal geometry when iterated or rescaled. Below are six domain-specific fractal constructions implemented in `fractals.js`.
## 1. Iterated Snap Basin Boundaries ("Lattice Julia Set")
The single-step snap map $T_\varepsilon(x) = x + \varepsilon\eta(x)$ followed by rounding $S(x) = \operatorname{round}(T_\varepsilon(x))$ is piecewise smooth with jump discontinuities along the half-integer pullbacks $T_\varepsilon^{-1}(\mathbb{Z}^2 + \tfrac12)$.
Instead of rounding immediately, iterate the continuous-plus-quantized feedback loop:
$$z_0 = x,\qquad z_{n+1} = z_n + \varepsilon\,\eta(z_n) - \gamma\big(z_n - \operatorname{round}(z_n)\big)$$
where $\gamma \in (0,1]$ is a lattice-gravity parameter. Points eventually lock into a periodic orbit around a specific lattice site $L_\infty(x) \in \mathbb{Z}^2$.
- **Primary scalar $v_1$:** Smoothed escape/locking time $n + \frac{\ln \tau}{\ln \|z_n - \operatorname{round}(z_n)\|}$, or the cumulative Galois norm accumulated along the trajectory.
- **Secondary scalar $v_2$:** A spatial hash of the terminal lattice site $L_\infty(x) = (m, n)$, coloring each lattice site's basin of attraction distinctly. Because $\eta$ is incommensurate with $\mathbb{Z}^2$, every basin boundary is crumpled by all $K$ modes simultaneously, producing a quasicrystalline Julia set.
## 2. Galois Conjugate Hofstadter Carpet (Minkowski Internal-Space IFS)
In Part II §6 & §10, an element $\xi = a + b\sqrt D \in \mathbb{Q}(\sqrt D)$ has two embeddings:
$$\sigma_+(\xi) = a + b\sqrt D \quad(\text{physical space}),\qquad \sigma_-(\xi) = a - b\sqrt D \quad(\text{internal / conjugate space}).$$
And the field norm is $N_{K/\mathbb{Q}}(\xi) = a^2 - D b^2$.
At each iteration $n = 1, \dots, M$:
1. Evaluate $(a_x, b_x, a_y, b_y) = \eta(z_n)$.
2. Compute the conjugate displacement $\bar r = (a_x - b_x\sqrt D,\ a_y - b_y\sqrt D)$ and the norm product $\mathcal{N}_n = (a_x^2 - D b_x^2) + (a_y^2 - D b_y^2)$.
3. **Fold** into the fundamental domain using the conjugate while scaling by the fundamental unit $\varepsilon_D$ (or $\sqrt D$):
   $$z_{n+1} = \sqrt D\,\big(\langle z_n + \varepsilon\, r(z_n)\rangle\big) + \bar r(z_n).$$
4. Count how many iterations land inside the **quasicrystal acceptance window** $|\bar r_\bullet| < W_{\text{acc}}$ and have negative norm $\mathcal{N}_n < 0$ (hyperbolic sectors).
The norm hyperbolas $a^2 - D b^2 = c$ slice through the multi-scale torus foldings, creating hyperbolic lacy carpets whose symmetries are governed by the unit group of $\mathbb{Z}[\sqrt D]$.
## 3. Metallic Weierstrass–Beatty Cascade (Exact $\mathbb{Z}[\sqrt D]$ Sturmian Fractal)
Addressing Part II §6 directly, this plot replaces $\sin/\cos$ with exact algebraic fractional parts and sums them over a geometric self-similar hierarchy of scales $\lambda_m = (\sqrt D)^m$:
$$W(x) = \sum_{m=0}^{M-1} D^{-m H} \, \Psi\!\left( (\sqrt D)^m R_{m\theta}\, x \right)$$
where $H \in (0,1)$ is the Hurst exponent, $R_{m\theta}$ rotates by the golden/plastic angle, and $\Psi$ is the exact Beatty sawtooth distance to the nearest half-integer seam:
$$\Psi(u) = \Big|\,\big|\langle u_1\sqrt D\rangle\big| - \big|\langle u_2\varphi\rangle\big|\,\Big|.$$
Because the discontinuities of $\langle t\rangle$ at integers propagate across all powers of $\sqrt D$, the graph of $W(x)$ has box-counting dimension $3 - H$ and exhibits exact inflation symmetry under $x \mapsto \sqrt D\, x$.
## 4. Diophantine Resonance Flame (Continued-Fraction Arnold Tongues)
Part II §10(d) notes that on a machine, "irrational" means "quantitatively far from small-denominator rationals." At each pixel $x$, form the local ratio of the irrational and rational field magnitudes:
$$\omega_x = \frac{a_x(x)}{b_x(x)\sqrt D + \delta},\qquad \omega_y = \frac{a_y(x)}{b_y(x)\sqrt D + \delta}.$$
Run the Euclidean / continued-fraction Gauss map $G(t) = \{1/t\}$ simultaneously on $(\omega_x, \omega_y)$ for $M$ steps, tracking the partial quotients $a_n = \lfloor 1/t_n \rfloor$:
- **Lyapunov / Brjuno sum:** $\Lambda(x) = \sum_{n=1}^M \frac{\ln(a_n^{(x)} + a_n^{(y)})}{2^n}$.
- **Denominator explosion:** where $\omega(x)$ crosses a low-order rational $p/q$, the Gauss map terminates early, burning dark **Arnold tongues** into a glowing web of badly-approximable noble numbers (where $a_n = 1$ or $2$ forever).
## 5. Precision-Horizon Catastrophe Map (Float64 vs. Exact $\mathbb{Z}[\sqrt D]$ Divergence)
A visual proof of Part II §3–4 ("rounding is not computable on $\mathbb{R}$"). For each pixel $x$, we simulate panning out to an extreme horizon $P = 2^B$ (for bit-depth $B \in [20, 53]$) and compare:
1. The IEEE-754 float64 evaluation $\operatorname{round}\big((x + P) + \varepsilon\eta_{\text{float}}(x + P)\big) - P$.
2. The exact `BigInt` algebraic sign-test evaluation in $\mathbb{Z}[\sqrt D]$ (using `isqrt(n²D)` from Part II §6).
Sweeping $B$ from $1$ to $53$ across each pixel and summing the bit-depth at which the snap decision flips produces a **Moiré-Cantor catastrophe dust** along the snap seams—directly imaging the measure-zero set where IEEE-754 betrays the continuum.
## 6. Strange Attractor Density of the Quasicrystal Kick Map
Treat the displacement field as a symplectic-like kick map on the torus $[ -L, L )^2$:
$$x_{n+1} = x_n + \varepsilon\, r(x_n),\qquad y_{n+1} = y_n + \varepsilon\, \bar r(x_{n+1})$$
coupling the real embedding $r = a + b\sqrt D$ with the Galois conjugate $\bar r = a - b\sqrt D$. Accumulate a 2D histogram of $10^6$ orbit steps, coloring by local Lyapunov divergence $\ln \|\det D\Phi(x)\|$.