"""Independent KKT evidence for the HDA optimization in fixed scales (ePC-SAFT issue #73).

PyNumero (PyomoNLP) evaluates the unscaled objective gradient, constraint Jacobian and
Lagrangian Hessian at the Ipopt primal point. Variables and constraints are scaled by their
physical units (hda_two_flash_epcsaft.unit_scale: T 300 K, P 1e5 Pa, flow 1 mol/s, duty
1e4 W, dimensionless 1) and the objective by 1e5 USD/yr. Multipliers are not taken from
Ipopt: with the four decisions as independent variables, the equality Jacobian gives the
null-space basis Z, the active inequality multipliers solve Z^T grad f = (A Z)^T mu by least
squares, and the equality multipliers follow from the state rows. Inequalities are h(x) >= 0.
"""

import numpy as np
import scipy.linalg
import scipy.sparse.linalg
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.environ import value
from pyomo.network import Arc

from idaes_examples.mod.hda.hda_two_flash_epcsaft import OBJECTIVE_SCALE, unit_scale

TOL = 1e-6  # scaled KKT norms and active-set slack
RANK_RATIO, CURVATURE, EXCLUDE = 1e-8, 1e-8, 1e-10


def kkt(m, decisions):
    """KKT residual norms, active set, active-Jacobian singular values and reduced-Hessian
    eigenvalues at the current primal point; `decisions` are the independent variables."""
    nlp = PyomoNLP(m)
    variables, constraints = nlp.get_pyomo_variables(), nlp.get_pyomo_constraints()
    sx = [unit_scale(v) for v in variables]
    sc = [unit_scale(c.body) for c in constraints]
    unclassified = [o.name for o, s in zip(variables + constraints, sx + sc) if s is None]
    sx = np.array([s or 1.0 for s in sx])
    sc = np.array([s or 1.0 for s in sc])
    x, g = nlp.get_primals(), nlp.evaluate_constraints()
    J = nlp.evaluate_jacobian().toarray() * sx[None, :] / sc[:, None]
    grad = nlp.evaluate_grad_objective() * sx / OBJECTIVE_SCALE
    gl, gu, xl, xu = nlp.constraints_lb(), nlp.constraints_ub(), nlp.primals_lb(), nlp.primals_ub()
    eq = gl == gu

    # Inequalities h >= 0 (scaled): process constraints, then variable bounds.
    rows, h, names, source, a_var = [], [], [], [], {}
    for i in np.flatnonzero(~eq):
        for bound, sign in ((gl[i], 1.0), (gu[i], -1.0)):
            if np.isfinite(bound):
                rows.append(sign * J[i]); h.append(sign * (g[i] - bound) / sc[i])
                names.append(constraints[i].name + (" >= lb" if sign > 0 else " <= ub")); source.append((i, sign))
    for j in range(len(x)):
        for bound, sign in ((xl[j], 1.0), (xu[j], -1.0)):
            if np.isfinite(bound):
                e = np.zeros(len(x)); e[j] = sign
                rows.append(e); h.append(sign * (x[j] - bound) / sx[j])
                names.append(variables[j].name + (" >= lb" if sign > 0 else " <= ub")); source.append(None)
                a_var[len(rows) - 1] = j
    A, h = np.array(rows), np.array(h)

    # Null space of the equality Jacobian with the decisions independent.
    position = {id(v): j for j, v in enumerate(variables)}
    dec = [position[id(v)] for v in decisions]
    st = [j for j in range(len(x)) if j not in dec]
    Je = J[eq]
    assert Je.shape[0] == len(st), "equality system is not square in the state variables"
    lu = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(Je[:, st]))
    Z = np.zeros((len(x), len(dec)))
    Z[st] = -lu.solve(Je[:, dec])
    Z[dec] = np.eye(len(dec))

    # Active set (owner rule M2): every inequality and bound with scaled slack <= TOL is reported.
    # Row j leaves the rank test only when its projection onto the decisions is shown to be zero:
    # ||a_j Z||_2 + err_j <= 1e-10 ||a_j||_2, err_j the change of a_j Z over one refinement step.
    dZ = np.zeros_like(Z)
    dZ[st] = -lu.solve(Je @ Z)
    cond = np.abs(Je[:, st]).sum(axis=0).max() * scipy.sparse.linalg.onenormest(  # 1-norm condition estimate
        scipy.sparse.linalg.LinearOperator(lu.shape, matvec=lu.solve, rmatvec=lambda v: lu.solve(v, trans="T")))
    Z = Z + dZ
    reduced_grad = Z.T @ grad
    near = np.flatnonzero(h <= TOL)
    proj, err = A[near] @ Z, np.abs(A[near] @ dZ).max(axis=1, initial=0.0)
    norm_a = np.linalg.norm(A[near], axis=1)
    excluded = np.linalg.norm(proj, axis=1) + err <= EXCLUDE * norm_a
    active = near[~excluded]
    Ar = A[active] @ Z
    mu = np.linalg.lstsq(Ar.T, reduced_grad, rcond=None)[0] if len(active) else np.zeros(0)
    stationarity = reduced_grad - Ar.T @ mu
    Aa = A[active]
    lam = lu.solve(grad[st] - Aa[:, st].T @ mu, trans="T")  # scaled equality multipliers
    full_residual = grad - Je.T @ lam - Aa.T @ mu

    primal = max(np.abs((g[eq] - gl[eq]) / sc[eq]).max(), max(0.0, -h.min()))
    comp = np.abs(mu * h[active]).max(initial=0.0)
    sv = scipy.linalg.svdvals(np.vstack([Je, Aa]))
    multiplier = dict(zip(active, mu))
    partner = arc_partners(m)
    near_rows = []
    for k, a in enumerate(near):
        row = {"name": names[a], "slack": float(h[a]), "row_norm": float(norm_a[k]),
               "projection": proj[k].tolist(), "projection_norm": float(np.linalg.norm(proj[k])),
               "projection_error": float(err[k]), "excluded": bool(excluded[k]),
               "multiplier": float(multiplier[a]) if a in multiplier else None}
        if source[a] is None:
            v = variables[a_var[a]]
            row.update(variable_identity(v), arc_partner=partner.get(id(v)))
        near_rows.append(row)

    # Reduced Lagrangian Hessian on the null space of the active constraints (unscaled duals).
    y = np.zeros(len(constraints))
    y[np.flatnonzero(eq)] = -lam / sc[eq]
    for k, a in enumerate(active):
        if source[a] is not None:
            i, sign = source[a]
            y[i] += -sign * mu[k] / sc[i]
    nlp.set_duals(y)
    nlp.set_obj_factor(1.0 / OBJECTIVE_SCALE)
    H = nlp.evaluate_hessian_lag().toarray() * sx[:, None] * sx[None, :]  # full symmetric
    N = scipy.linalg.null_space(Ar) if len(active) else np.eye(len(dec))
    eig = np.linalg.eigvalsh(N.T @ Z.T @ H @ Z @ N) if N.shape[1] else np.zeros(0)

    norms = {"primal_feasibility": float(primal), "stationarity": float(np.abs(stationarity).max()),
             "full_stationarity": float(np.abs(full_residual).max()),
             "dual_sign": float(max(0.0, -mu.min(initial=0.0))), "complementarity": float(comp)}
    full_rank = bool(len(active) <= len(dec) and sv.min() >= RANK_RATIO * sv.max())
    stationary = all(v <= TOL for v in norms.values()) and full_rank
    # M4: with as many counted active rows as decisions and every multiplier > TOL the critical
    # cone is {0}; otherwise A5's reduced-Hessian eigenvalue test.
    empty_cone = len(active) == len(dec) and mu.min(initial=-np.inf) > TOL
    local_min = stationary and (empty_cone or (eig.size > 0 and eig.min() > CURVATURE))
    return {
        "norms": norms,
        "near_active": near_rows,
        "counted_active": len(active), "excluded": int(excluded.sum()),
        "state_jacobian_condition_1norm": float(cond),
        "singular_values": {"min": float(sv.min()), "max": float(sv.max()), "ratio": float(sv.min() / sv.max()),
                            "rows": int(Je.shape[0] + len(active)), "columns": int(Je.shape[1]),
                            "full_row_rank": full_rank},
        "reduced_hessian_eigenvalues": eig.tolist(), "reduced_gradient": reduced_grad.tolist(),
        "label": "local minimum" if local_min else ("stationary candidate" if stationary else "not stationary"),
        "unclassified_scales": unclassified,
        "equality_multiplier_max_abs": float(np.abs(lam).max()),
        "_primal_dual": {"variables": [(v.name, float(xv), float(s)) for v, xv, s in zip(variables, x, sx)],
                         "equalities": [(constraints[i].name, float(l), float(sc[i]))
                                        for i, l in zip(np.flatnonzero(eq), lam)]},
    }


def variable_identity(v):
    """Declared and effective bounds of a near-active variable; for a phase quantity both the raw
    phase flow and phase_frac x flow_mol (two stored values linked by one equality)."""
    out = {"domain": str(v.domain), "declared_bounds": [None if b is None else value(b) for b in (v.lower, v.upper)],
           "effective_bounds": list(v.bounds), "value": float(v.value)}
    sb, name = v.parent_block(), v.parent_component().local_name
    if name in ("flow_mol_phase", "phase_frac"):
        p = v.index()
        out["flow_mol_phase"] = float(value(sb.flow_mol_phase[p]))
        out["phase_frac_x_flow_mol"] = float(value(sb.phase_frac[p] * sb.flow_mol))
    return out


def arc_partners(m):
    """id(variable) -> name of the variable an arc equality identifies it with."""
    out = {}
    for arc in m.component_data_objects(Arc):
        for key, src in arc.source.vars.items():
            dst = arc.destination.vars[key]
            if src.is_indexed() and dst.is_indexed():
                for i in src:
                    if i in dst and src[i].is_variable_type() and dst[i].is_variable_type():
                        out[id(src[i])], out[id(dst[i])] = dst[i].name, src[i].name
    return out
