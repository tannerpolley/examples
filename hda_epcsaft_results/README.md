# HDA two-flash with ePC-SAFT properties: phase-root qualification (ePC-SAFT issue #73)

Command, run from this checkout's root with the Engine environment's Python. `PYTHONPATH` must
point at this checkout, and only the installed non-editable `epcsaft` wheel may be imported:

    OMP_NUM_THREADS=1 PYTHONPATH=. python -m idaes_examples.mod.hda.hda_epcsaft_qualification

`qualification_meta.json` records the wheel, the wheel and ASL-library SHA-256 values, the
input-file SHA-256 values, and the proxy stream states.

## Files

- `ideal_baseline.json`: the corrected ideal baseline, an adiabatic reactor at fixed 75 %
  conversion. It holds the fixed simulation and the optimization from the three fixed starts.
- `property_transport.csv`: IDAES ePC-SAFT expressions (density, total enthalpy, ln phi_i
  assembled through `log_fug_phase_comp_eq`) compared with `Mixture.state` at the issue's
  qualification states.
- `callback_domain.csv`: availability of the ASL callback versus temperature.
- `bubble_dew_idaes.csv`: the IDAES `LogBubbleDew` bubble and dew equations solved alone on
  single ePC-SAFT state blocks, from a fixed grid of starting temperatures, at the ideal-baseline
  stream compositions. Each converged row is audited with the public API: fractions and sums,
  distinct roots, composition distance of at least 1e-3, the declared liquid having the higher
  mass density, and fugacity residual <= 1e-8.
- `bubble_sum_at_stream_pressure.csv`: the bubble test S = sum K_i z_i at each hydrogen-rich
  stream's pressure, from 12 to 300 K, with the incipient vapor converged by successive
  substitution through the public API.

All files come from one fresh process on the final Engine wheel:
- wheel SHA-256 `60b15ce85a065fdbc98d40c9597f6e93775d39fbefd613866410804287dd9854`
- `libepcsaft_asl.so` SHA-256 `af60d4755193d3bf2b686d4b416d0acbbbbdc43da31fa4cc3d90342e9c0073d9`
- parameters SHA-256 `c41ac486...6822`, NASA9 records SHA-256 `d1341f07...c39b`

## Results

- **Property transport.** The IDAES ePC-SAFT package reproduces `Mixture.state` density, total
  enthalpy and all four ln phi_i at the four issue states, both phases. The maximum relative
  difference is 1.9e-12.
- **Callback domain.** Density and ln phi_i are available below 200 K. Enthalpy is NaN at
  T <= 200 K (NASA9 lower bound, `unavailable_derivative:OutsideInterval`).
- **Ideal baseline.**
  - Fixed simulation: reactor outlet 771.8571 K, purity 0.824296.
  - Optimization, all three starts optimal: 312674.2456 USD/yr at H101 500 K, F101 301.881 K,
    F102 362.935 K and 105000 Pa; reactor outlet 698.61 K.
  - KKT evidence was not computed.

## Supported negative at numerical-study step 4 (issue #73 stop rule)

With all four species in both phases, `SmoothVLE` needs a bubble temperature,
`_t1 = smooth_max(T, T_bub)`, for every stream. The hydrogen-rich HDA streams have none at their
pressure.

**IDAES bubble equations (`bubble_dew_idaes.csv`).** Starts ran from 80 to 400 K with a light
and a mixed incipient-vapor guess, for the gas feed, M101 outlet, R101 outlet (the F101 feed)
and F101 vapor.
- Every converged row is rejected.
- Almost all are the trivial solution: |x - y| <= 1e-8, with one shared root or two
  single-root states of equal density. Their temperatures are arbitrary, between 21 and
  1497 K, for example gas feed 21.04 / 238.7 / 1486-1497 K, M101 outlet 234.5 / 541.7 /
  569.8 K, R101 outlet 30.7 / 242.0 / 314.5 K, F101 vapor 187.3 / 219.7 K.
- The one non-trivial candidate (gas feed at 19.74 K) fails the fugacity residual (135). It is
  also cryogenic, below the methane triple point.
- The failure is not a 200 K floor: starts at 80-160 K reach the same outcomes.

**Bubble sum at the stream pressure (`bubble_sum_at_stream_pressure.csv`).** With the liquid
at the stream composition z and 0.35 MPa, the incipient vapor is converged by successive
substitution from a hydrogen-rich start, and S = sum K_i z_i is recorded from 12 to 300 K.
S > 1 means the liquid is unstable to that vapor; a bubble temperature needs S = 1.

| Stream | S where a liquid-like root exists | Higher temperatures |
|---|---|---|
| M101 outlet (28 % H2) | 78 to 2.9e9 (12-210 K) | only one gas-like root (240-300 K) |
| R101 outlet / F101 feed (16 % H2) | 34 to 5.8e9 (12-240 K) | only one gas-like root (270-300 K) |
| F101 vapor (19 % H2) | 28 to 8.7e4 (12-180 K) | only one gas-like root (210-300 K) |

So these process streams have no bubble temperature at their pressure: dissolving 16-28 % H2
needs tens of MPa of hydrogen (the retained H2/CH4 observations show x_H2 = 0.12-0.19 only at
8-10 MPa). The gas feed (94 % H2) is gas-like above 40 K; below it, hydrogen condensation near its
33 K critical point dominates. Its only non-trivial IDAES bubble candidate (19.74 K) is liquid
hydrogen holding about 6 % methane, below the methane triple point, with incipient aromatic
fractions at the 1e-12 floor and audit fugacity residual 135. It lies outside the model and HDA
domains and is an unavailable auxiliary root under the issue's rule. An independent reviewer
located a model bubble root at that temperature.

**What does exist.**
- Accepted dew rows (distinct roots, audit passed): gas feed 227.65 K, M101 outlet 366.284 K,
  R101 outlet 350.911 K, F101 vapor 324.08 K.
- Aromatic bubble and dew points converge with distinct Liquid/Vapor roots:
  - liquid feed: bubble 432.728 K, dew 433.614 K
  - F102 feed: bubble 406.604 K, dew 411.198 K
  - F102 vapor: bubble 371.091 K, dew 375.283 K
  - F102 liquid: bubble 375.487 K, dew 381.268 K
- These aromatic rows have IDAES equation residuals near 5e-11, and their benzene/toluene audit
  residuals are near 5e-11. The audit still rejects them, because the trace H2/CH4 residuals
  are 1e-5 to 0.25. The trace fractions are ~1e-11 (the 1e-9 floor on the ideal-proxy traces),
  and IDAES solves `exp(ln x) = x` to an absolute tolerance, so ln x of a trace species is
  inaccurate there. This is a scaling limit of trace species in the proxy compositions, not a
  missing root. The liquid-feed dew row also falls under the 1e-3 composition-distance rule,
  because near-pure toluene gives |x - y| = 9.6e-5.

## Claim limits and owner options

- The ePC-SAFT simulation, optimization, KKT evidence and smoothing sensitivity were not run.
  The formulation has no bubble temperature for the hydrogen-rich streams, so the issue stops
  before study steps 5-6.
- The compositions are proxies taken from the ideal baseline.
- Benzene/toluene k_ij = -0.003 with the Gross pure set is provisional (the rounding check
  failed on liquid ln phi). The hydrogen pairs are the final fitted values.
- There is no claim of global phase stability or of accuracy outside the observation ranges.
- Owner options:
  - Vapor-only H2/CH4, the standard IDAES treatment, requires the callback to accept
    structural zero weights, which the issue forbids.
  - A complementarity phase-equilibrium formulation with no bubble or dew temperature changes
    the issue's scientific question.
