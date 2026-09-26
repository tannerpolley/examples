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

from idaes_examples.mod.hda.hda_two_flash_epcsaft import OBJECTIVE_SCALE, unit_scale

TOL = 1e-6  # scaled KKT norms and active-set slack
RANK_RATIO, CURVATURE = 1e-8, 1e-8


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
    rows, h, names, source = [], [], [], []
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
    reduced_grad = Z.T @ grad

    # Active set: scaled slack <= TOL, measured in the decision space (slack / max-norm of the
    # constraint's reduced row), which equals the plain scaled slack for decision bounds. A trace
    # quantity next to its bound (smoothing- or eps-sized phase amounts, vapor-only fractions,
    # 1e-12 to 1e-7) has a small plain slack but is set by the equalities and lies a finite decision
    # change away from its bound: it is inactive, with zero multiplier, and reported separately.
    near = np.flatnonzero(h <= TOL)
    reach = np.abs(A[near] @ Z).max(axis=1, initial=0)
    distance = np.divide(h[near], reach, out=np.full(near.size, np.inf), where=reach > 0)
    active = near[distance <= TOL]
    inactive_near = near[distance > TOL]
    Ar = A[active] @ Z
    mu = np.linalg.lstsq(Ar.T, reduced_grad, rcond=None)[0] if len(active) else np.zeros(0)
    stationarity = reduced_grad - Ar.T @ mu
    Aa = A[active]
    lam = lu.solve(grad[st] - Aa[:, st].T @ mu, trans="T")  # scaled equality multipliers

    primal = max(np.abs((g[eq] - gl[eq]) / sc[eq]).max(), max(0.0, -h.min()))
    comp = np.abs(mu * h[active]).max(initial=0.0)
    sv = scipy.linalg.svdvals(np.vstack([Je, Aa]))

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
             "dual_sign": float(max(0.0, -mu.min(initial=0.0))), "complementarity": float(comp)}
    full_rank = bool(sv.min() >= RANK_RATIO * sv.max())
    stationary = all(v <= TOL for v in norms.values()) and full_rank
    # With as many independent active constraints as decisions the critical cone is {0} under
    # strict complementarity, so positive multipliers replace the (empty) curvature test.
    local_min = stationary and (eig.min() > CURVATURE if eig.size else mu.min(initial=np.inf) > TOL)
    return {
        "norms": norms,
        "active_set": [{"name": names[a], "slack": float(h[a]), "multiplier": float(mu[k])}
                       for k, a in enumerate(active)],
        "near_bound_inactive": {"count": int(inactive_near.size),
                                "min_decision_distance": float(distance[distance > TOL].min(initial=np.inf)),
                                "names": [names[a] for a in inactive_near]},
        "singular_values": {"min": float(sv.min()), "max": float(sv.max()), "ratio": float(sv.min() / sv.max()),
                            "full_row_rank": full_rank},
        "reduced_hessian_eigenvalues": eig.tolist(), "reduced_gradient": reduced_grad.tolist(),
        "label": "local minimum" if local_min else ("stationary candidate" if stationary else "not stationary"),
        "unclassified_scales": unclassified,
        "equality_multiplier_max_abs": float(np.abs(lam).max()),
        "_primal_dual": {"variables": [(v.name, float(xv), float(s)) for v, xv, s in zip(variables, x, sx)],
                         "equalities": [(constraints[i].name, float(l), float(sc[i]))
                                        for i, l in zip(np.flatnonzero(eq), lam)]},
    }
