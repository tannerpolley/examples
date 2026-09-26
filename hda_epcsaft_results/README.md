# HDA two-flash with ePC-SAFT properties: phase-root qualification (ePC-SAFT issue #73)

Command, run from this checkout's root with the Engine environment's Python. `PYTHONPATH` must
point at this checkout, and only the installed non-editable `epcsaft` wheel may be imported:

    OMP_NUM_THREADS=1 PYTHONPATH=. python -m idaes_examples.mod.hda.hda_epcsaft_qualification

`qualification_meta.json` records the wheel, the wheel and ASL-library SHA-256 values, the
input-file SHA-256 values, and the proxy stream states.

## Files

- `ideal_baseline.json`: the corrected ideal baseline, an adiabatic reactor at fixed 75 %
  conversion. It holds the fixed simulation and the optimization from the three fixed starts.
- `property_transport.csv` (the last valid run found max relative difference 1.9e-12): IDAES ePC-SAFT expressions (density, total enthalpy, ln phi_i
  assembled through `log_fug_phase_comp_eq`) compared with `Mixture.state` at the issue's
  qualification states.
- `callback_domain.csv`: availability of the ASL callback versus temperature.
- `bubble_dew_idaes.csv`: the IDAES `LogBubbleDew` bubble and dew equations solved alone on
  single ePC-SAFT state blocks, from a fixed grid of starting temperatures, at the ideal-baseline
  stream compositions. Each converged row is audited with the public API: fractions and sums,
  distinct roots, composition distance of at least 1e-3, the declared liquid having the higher
  mass density, and fugacity residual <= 1e-8.
- `bubble_pressure_engine_diagnostic.csv`: a diagnostic only, never used by the IDAES path.
  The Engine equilibrium driver computes the bubble pressure of the hydrogen-rich compositions
  at fixed temperature, and each result passes the same audit.

Only `ideal_baseline.json` is retained in this commit. The Engine-dependent files are not:
the Engine wheel was reinstalled during the last full run, so its rows mixed libraries. They
must be regenerated with the command above on the final wheel (SHA-256 55eff6a0...a31f4,
`libepcsaft_asl.so` 7aa66c1d...8f0f).

## Result: supported negative at numerical-study step 4 (issue #73 stop rule)

With all four species admitted to both phases, `SmoothVLE` builds a bubble temperature
(`_t1 = smooth_max(T, T_bub)`) for every stream. The hydrogen-rich HDA streams have no genuine
bubble temperature at their pressure:

- IDAES `LogBubbleDew` bubble equations, from every start at or above 210 K, converge only to
  the trivial solution y = x on one shared `Unspecified` root. The temperature of that solution
  is arbitrary: gas feed 200-201 K (1488 K from other guesses), M101 outlet 234.5 K, reactor
  outlet 241.8 K, F101 vapor 200.0 K; max |x - y| < 5e-9. Checked on wheel 55eff6a0. The audit
  rejects every such row on the distinct-root, composition-distance, mass-ordering and residual
  checks.
- Engine equilibrium-driver diagnostic (bubble pressure at fixed T): the reactor-outlet /
  F101-feed composition (H2 0.161, CH4 0.626) has accepted bubble pressures of 107 MPa at
  300 K and 78 MPa at 350 K, and no accepted root from 20 to 250 K. The gas feed has no accepted
  root from 20 to 350 K. A bubble temperature at 350 kPa would require P_bub = 0.35 MPa.
- Successive-substitution diagnostic: sum(x_i K_i) at 350 kPa is 74 to 1.7e8 for the F101 feed
  from 15 to 150 K (distinct liquid and vapor roots where two exist). Above that range only the
  trivial single root remains.
- The ASL callback returns NaN for every observable, including density and ln phi_i, at
  T <= 200 K. The causes are `unavailable_derivative:OutsideInterval:0` and, at 200 K,
  `ThermalBoundary:0`: the NASA9 lower bound, although density and ln phi need no ideal-gas
  record. Any cryogenic auxiliary root is therefore unavailable through the callback as well.

Dew temperatures exist, with distinct roots: gas feed 227.65 K, M101 outlet 366.3 K, F101 feed
350.9 K, F101 vapor 324.1 K, F102 streams 375-381 K. The aromatic bubble point also exists
(F102 feed 372.6 K, Liquid/Vapor roots). Only the hydrogen-rich bubble branch fails.

## Claim limits

- The compositions are proxies taken from the ideal baseline. The ePC-SAFT flowsheet was not
  solved, because study step 4 stops before steps 5 and 6.
- Benzene/toluene k_ij = -0.003 with the Gross pure set is provisional: the rounding check
  failed on liquid ln phi. The hydrogen pairs are the final fitted values.
- There is no claim of global phase stability or of accuracy outside the observation ranges.
