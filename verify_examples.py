#!/usr/bin/env python3
"""Reproduce the paper's mathematical examples. Python 3.10+, standard library only.
All outputs are synthetic mathematical calculations, not experimental observations.
"""
from __future__ import annotations
import json
import cmath
import math
from fractions import Fraction as F
from pathlib import Path
from typing import Sequence


def determinant(matrix: Sequence[Sequence[int | F]]) -> F:
    a = [[F(x) for x in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("A square matrix is required")
    det = F(1)
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            det = -det
        p = a[i][i]
        det *= p
        for j in range(i + 1, n):
            q = a[j][i] / p
            for k in range(i + 1, n):
                a[j][k] -= q * a[i][k]
    return det


def solve(matrix: Sequence[Sequence[int | F]], rhs: Sequence[int | F]) -> list[F]:
    n = len(matrix)
    if len(rhs) != n or any(len(row) != n for row in matrix):
        raise ValueError("Incompatible dimensions")
    a = [[F(x) for x in row] + [F(rhs[i])] for i, row in enumerate(matrix)]
    for i in range(n):
        p = next((j for j in range(i, n) if a[j][i]), None)
        if p is None:
            raise ValueError("Singular matrix")
        a[i], a[p] = a[p], a[i]
        scale = a[i][i]
        a[i] = [x / scale for x in a[i]]
        for j in range(n):
            if j != i:
                scale = a[j][i]
                a[j] = [x - scale*y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def laplacian(edges: Sequence[tuple[int, int]], n: int = 6) -> list[list[int]]:
    a = [[0]*n for _ in range(n)]
    for i, j in edges:
        if not (0 <= i < n and 0 <= j < n) or i == j or a[i][j]:
            raise ValueError("Edges must define a simple undirected graph")
        a[i][j] = a[j][i] = -1
    for i in range(n):
        a[i][i] = -sum(a[i])
    return a


def theta_f(a: float, exponent_cutoff: float = 45.0) -> tuple[float, float]:
    """Return Theta(a) and F(a)=a sum(n^2 exp(-a n^2))/Theta(a).
    The Poisson transform avoids long sums and small-a cancellation in Theta.
    """
    if not math.isfinite(a) or a <= 0:
        raise ValueError("a must be finite and positive")
    if a < 1:
        theta, f = theta_f(math.pi**2/a, exponent_cutoff)
        return math.sqrt(math.pi/a)*theta, 0.5-f
    nmax = max(1, math.ceil(math.sqrt(exponent_cutoff/a)))
    theta = 1.0 + 2.0*math.fsum(math.exp(-a*n*n) for n in range(1, nmax+1))
    moment = 2.0*math.fsum(n*n*math.exp(-a*n*n) for n in range(1, nmax+1))
    return theta, a*moment/theta


def spectral_dimension(t: float, epsilon: float, cutoff: float = 45.0) -> float:
    if not 0 < epsilon < 1:
        raise ValueError("epsilon must be in (0,1)")
    return 4.0*(theta_f(t, cutoff)[1] + theta_f(t/epsilon**2, cutoff)[1])


def generalized_dimension(p: float, q: float) -> float:
    if not 0 < p < 1:
        raise ValueError("p must be in (0,1)")
    if q == 1:
        return 3.0 - (p*math.log(p)+(1-p)*math.log(1-p))/math.log(2)
    return 3.0 + math.log(p**q+(1-p)**q)/((1-q)*math.log(2))


def resolvent_trace(l: Sequence[Sequence[int]], kappa: F) -> F:
    """Exact trace from rational linear solves, independently of spectral formula."""
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    n = len(l)
    matrix = [[F(l[i][j]) + kappa*(i == j) for j in range(n)] for i in range(n)]
    return sum((solve(matrix, [F(i == j) for i in range(n)])[j]
                for j in range(n)), F(0))


def inner(a: Sequence[complex], b: Sequence[complex]) -> complex:
    if len(a) != len(b):
        raise ValueError("Inner-product dimensions must agree")
    return sum((complex(x).conjugate()*y for x, y in zip(a, b)), 0j)


def finite_detector_vectors(xi: float, theta: float) -> tuple[list[complex], list[complex]]:
    """Unnormalized vectors from the single, entrywise nonnegative relation kernel."""
    if not math.isfinite(xi) or not 0.0 <= xi <= 1.0 or not math.isfinite(theta):
        raise ValueError("xi must be in [0,1] and theta finite")
    a, b = math.sqrt(xi)/2.0, math.sqrt(1.0-xi)/2.0
    return [a, a, b, 0j], [a, a*cmath.exp(1j*theta), 0j, b]


def verify_amendments(graphs: dict) -> dict:
    """Verify the first-amendment trace and joint-kernel constructions.

    Exact special-phase tests and direct floating-point inner products are separate.
    Nothing here is experimental data or a universal ECL derivation.
    """
    traces = []
    kappas = [F(1,10), F(1,2), F(1), F(2), F(4), F(7,3), F(10), F(100)]
    for kappa in kappas:
        r_k = resolvent_trace(graphs["K33"]["laplacian"], kappa)
        r_p = resolvent_trace(graphs["prism"]["laplacian"], kappa)
        assert r_k == 1/kappa + 4/(kappa+3) + 1/(kappa+6)
        assert r_p == 1/kappa + 1/(kappa+2) + 2/(kappa+3) + 2/(kappa+5)
        gap = F(12)/((kappa+2)*(kappa+3)*(kappa+5)*(kappa+6))
        assert r_p-r_k == gap > 0
        traces.append({"kappa": str(kappa), "trace_K": str(r_k),
                       "trace_prism": str(r_p), "gap": str(gap)})

    exact_count = 0
    rational_xis = sorted(set([F(n,16) for n in range(17)] + [F(1,3)]))
    quarter_phases = [(F(1),F(0)), (F(0),F(-1)), (F(-1),F(0)), (F(0),F(1))]
    exact_rows = []
    for xi in rational_xis:
        rows = [[xi,xi,1-xi,0], [xi,xi,0,1-xi]]
        eta = F(1,4)*sum(rows[0])
        assert eta == F(1,4)*sum(rows[1]) == (1+xi)/4
        assert F(1,4) <= eta <= F(1,2)
        assert 4*eta-1 == xi
        for phase_index, (real_e, imag_e) in enumerate(quarter_phases):
            gamma_re = xi*(1+real_e)/(1+xi)
            gamma_im = xi*imag_e/(1+xi)
            assert 4*eta*gamma_re == (4*eta-1)*(1+real_e)
            assert 4*eta*gamma_im == (4*eta-1)*imag_e
            v2 = gamma_re**2 + gamma_im**2
            assert 0 <= v2 <= 1
            if phase_index == 0:
                assert gamma_re == 2-1/(2*eta) and gamma_im == 0
            exact_count += 1
            if (xi,phase_index) in [(F(0),0),(F(1,3),0),(F(1,2),0),
                                    (F(1,2),1),(F(1,2),2),(F(1),0),(F(1),2)]:
                exact_rows.append({"xi":str(xi),"eta":str(eta),
                    "theta_in_units_pi_over_2":phase_index,
                    "gamma_real":str(gamma_re),"gamma_imag":str(gamma_im),
                    "visibility_squared":str(v2)})

    max_norm_error = max_overlap_error = max_isometry_error = max_port_error = 0.0
    grid_count = port_count = 0
    for n in range(21):
        xi = n/20.0
        eta = (1+xi)/4.0
        for m in range(17):
            theta = -math.pi + m*math.pi/8.0
            a1,a2 = finite_detector_vectors(xi,theta)
            max_norm_error=max(max_norm_error,abs(inner(a1,a1)-eta),abs(inner(a2,a2)-eta))
            d1,d2 = [[z/math.sqrt(eta) for z in a] for a in (a1,a2)]
            gamma = inner(d2,d1)
            expected = (1-1/(4*eta))*(1+cmath.exp(-1j*theta))
            max_overlap_error=max(max_overlap_error,abs(gamma-expected))
            assert abs(gamma) <= 1+1e-14
            j1 = a1+[0j]*4+[math.sqrt(1-eta),0j]
            j2 = [0j]*4+a2+[0j,math.sqrt(1-eta)]
            max_isometry_error=max(max_isometry_error,abs(inner(j1,j1)-1),
                                   abs(inner(j2,j2)-1),abs(inner(j1,j2)))
            permutation=(2,0,3,1)
            assert abs(inner([d2[i] for i in permutation],
                             [d1[i] for i in permutation])-gamma)<1e-14
            # Direct recombination is checked against the phase-referenced formula.
            for prior in (0.0,0.1,0.5,0.9,1.0):
                for phi in (0.0,0.25,math.pi/2,math.pi):
                    plus=[(math.sqrt(prior)*u+cmath.exp(1j*phi)*math.sqrt(1-prior)*v)/math.sqrt(2)
                          for u,v in zip(d1,d2)]
                    minus=[(math.sqrt(prior)*u-cmath.exp(1j*phi)*math.sqrt(1-prior)*v)/math.sqrt(2)
                           for u,v in zip(d1,d2)]
                    p_plus=inner(plus,plus).real
                    p_minus=inner(minus,minus).real
                    formula=(1+2*math.sqrt(prior*(1-prior))*(cmath.exp(-1j*phi)*gamma).real)/2
                    max_port_error=max(max_port_error,abs(p_plus-formula),abs(p_plus+p_minus-1))
                    assert -1e-14 <= p_plus <= 1+1e-14 and -1e-14 <= p_minus <= 1+1e-14
                    port_count += 1
            grid_count += 1
    assert max(max_norm_error,max_overlap_error,max_isometry_error,max_port_error)<1e-12

    # Exact calibration-to-prediction on the two observable channels.
    eta_calibration=F(3,8)
    xi_from_calibration=4*eta_calibration-1
    prediction=2-1/(2*eta_calibration)
    foil=F(9,10)
    assert xi_from_calibration==F(1,2) and prediction==F(2,3)
    assert foil-prediction==F(7,30)
    eta_error=F(1,1000)
    measurement_error=F(1,200)
    tolerance=8*eta_error+measurement_error
    assert foil-prediction > tolerance
    # The broader-quantum foil has a normalized pointer and the stated overlap.
    d1=[1+0j,0j]
    d2=[0.9+0j,math.sqrt(19)/10]
    assert abs(inner(d2,d2)-1)<1e-14 and abs(inner(d2,d1)-0.9)<1e-14

    return {"status":"All amendment checks passed",
        "resolvent_traces":traces,"joint_kernel_exact_rows":exact_rows,
        "exact_quarter_turn_cases":exact_count,"complex_grid_cases":grid_count,
        "port_probability_cases":port_count,"max_norm_error":max_norm_error,
        "max_overlap_error":max_overlap_error,"max_isometry_error":max_isometry_error,
        "max_port_error":max_port_error,
        "held_out_test":{"accessibility_calibration":str(eta_calibration),
            "xi_inferred_without_coherence":str(xi_from_calibration),
            "predicted_visibility_theta_zero":str(prediction),
            "hypothetical_foil_visibility":str(foil),"gap":str(foil-prediction),
            "example_error_allowance":str(tolerance),"foil_rejected":True},
        "scope":"Synthetic checks of a stipulated finite representation, not empirical ECL evidence"}



def _bounds(interval: tuple[F, F], name: str, *, positive: bool = False,
            unit: bool = False) -> tuple[F, F]:
    """Validate exact finite bounds before interval propagation."""
    lo, hi = map(F, interval)
    if lo > hi or lo < 0 or (positive and lo <= 0) or (unit and hi > 1):
        raise ValueError(f"Invalid {name} interval")
    return lo, hi


def absolute_accessibility_interval(q: tuple[F, F], efficiency: tuple[F, F]
                                   ) -> tuple[F, F] | None:
    """Conservative corrected interval; None means joint-package domain rejection.

    Uncertainty intervals must hold simultaneously. Unknown efficiency is not
    estimated from the target outcome. A zero efficiency lower bound is unresolved.
    """
    q_lo, q_hi = _bounds(q, "rate", unit=True)
    e_lo, e_hi = _bounds(efficiency, "efficiency", positive=True, unit=True)
    lo, hi = max(F(1,4), q_lo/e_hi), min(F(1,2), q_hi/e_lo)
    return None if lo > hi else (lo, hi)


def relative_accessibility_interval(q_all: tuple[F, F], q_shared: tuple[F, F],
                                    efficiency_ratio: tuple[F, F]
                                   ) -> tuple[F, F] | None:
    """Conservative same-kernel gate-ratio interval; no near-zero division."""
    a_lo, a_hi = _bounds(q_all, "all-gate rate", positive=True, unit=True)
    c_lo, c_hi = _bounds(q_shared, "shared-gate rate", unit=True)
    z_lo, z_hi = _bounds(efficiency_ratio, "efficiency ratio", positive=True)
    lo, hi = max(F(0), c_lo/(z_hi*a_hi)), min(F(1), c_hi/(z_lo*a_lo))
    return None if lo > hi else (lo, hi)


def verify_second_amendment() -> dict:
    """Independent exact and finite-grid checks of measurement-aware calibration.

    This is synthetic validation of the stated instrument and assumptions, not
    evidence that the reference gates or scalar efficiency apply to an experiment.
    """
    rational_xis = sorted(set([F(n,32) for n in range(33)] + [F(1,3)]))
    efficiencies = [F(1), F(9,10), F(1,5), F(1,1000)]
    eps = F(1,100000)
    common_cases = mismatch_cases = interval_cases = 0
    for xi in rational_xis:
        eta, shared = (1+xi)/4, xi/2
        rho = 2*xi/(1+xi)
        assert 0 <= rho <= 1 and rho/(2-rho) == xi
        assert 1/(2*(2-rho)) == eta and (2-1/(2*eta)) == rho
        for e_all in efficiencies:
            q_all, q_shared = e_all*eta, e_all*shared
            assert q_shared/q_all == rho and q_all/e_all == eta
            assert 4*q_all/e_all-1 == xi
            qi=(max(F(0),q_all-eps), min(F(1),q_all+eps))
            ei=(max(e_all/2,e_all-eps), min(F(1),e_all+eps))
            ai=absolute_accessibility_interval(qi,ei)
            assert ai is not None and ai[0] <= eta <= ai[1]
            common_cases += 1; interval_cases += 1
            for e_shared in efficiencies:
                zeta=e_shared/e_all
                qc=e_shared*shared
                corrected=qc/(zeta*q_all)
                assert corrected == rho
                assert abs(qc/q_all-rho) == abs(zeta-1)*rho <= abs(zeta-1)
                qci=(max(F(0),qc-eps), min(F(1),qc+eps))
                zi=(zeta*(1-eps),zeta*(1+eps))
                ri=relative_accessibility_interval(qi,qci,zi)
                assert ri is not None and ri[0] <= rho <= ri[1]
                mismatch_cases += 1; interval_cases += 1

    # Calibration domains are checked after correction, never from a raw low rate.
    domain_tests = [
        (absolute_accessibility_interval((F(1,5),F(1,5)),(F(1),F(1))), None),
        (absolute_accessibility_interval((F(3,5),F(3,5)),(F(1),F(1))), None),
        (absolute_accessibility_interval((F(1,20),F(1,20)),(F(1,10),F(1,10))), (F(1,2),F(1,2))),
        (absolute_accessibility_interval((F(6,25),F(13,50)),(F(1),F(1))), (F(1,4),F(13,50))),
        (relative_accessibility_interval((F(1,5),F(1,5)),(F(2,5),F(2,5)),(F(1),F(1))),None),
        (relative_accessibility_interval((F(27,80),F(27,80)),(F(9,40),F(9,40)),(F(1),F(1))), (F(2,3),F(2,3))),
        (relative_accessibility_interval((F(3,80),F(3,80)),(F(1,4),F(1,4)),(F(10),F(10))), (F(2,3),F(2,3))),
        (relative_accessibility_interval((F(1,4),F(1,4)),(F(0),F(0)),(F(1),F(1))), (F(0),F(0))),
    ]
    for actual, expected in domain_tests: assert actual == expected
    invalid_cases=[
        lambda: absolute_accessibility_interval((F(0),F(1,2)),(F(0),F(1))),
        lambda: relative_accessibility_interval((F(0),F(1,2)),(F(0),F(1,2)),(F(1),F(1))),
        lambda: relative_accessibility_interval((F(1,4),F(1,4)),(F(0),F(1,2)),(F(0),F(1))),
        lambda: absolute_accessibility_interval((F(1,3),F(1,4)),(F(1),F(1))),
    ]
    for operation in invalid_cases:
        try: operation()
        except ValueError: pass
        else: raise AssertionError("Invalid or unresolved interval was not rejected")

    # Exact 10% erasure calculation and nontrivial held-out prediction.
    eta=F(3,8); xi=F(1,2); e=F(9,10)
    q_all=e*eta; q_shared=e*xi/2; rho=q_shared/q_all
    naive_xi=4*q_all-1; naive_visibility=2-1/(2*q_all)
    assert q_all==F(27,80) and q_shared==F(9,40)
    assert naive_xi==F(7,20) and xi-naive_xi==F(3,20)
    assert naive_visibility==F(14,27) and rho-naive_visibility==F(4,27)
    assert rho==F(2,3) and 8*(eta-q_all)==F(3,10)
    foil=F(9,10)
    tolerance=F(1,100)+F(1,2)*F(1,100)+F(1,200)
    assert foil-rho==F(7,30)>tolerance

    # The four atomic distances realize the weighted Definition 6.1 exactly.
    path_points=[(-2,0,0,0),(2,0,0,0)]
    alternative_points=[(-1,1,0,0),(1,1,0,0),(-1,-1,0,0),(1,-1,0,0)]
    d2=lambda a,b: sum((F(x)-F(y))**2 for x,y in zip(a,b))
    distances=[[d2(p,a) for a in alternative_points] for p in path_points]
    assert distances==[[F(2),F(10),F(2),F(10)],[F(10),F(2),F(10),F(2)]]
    ell=F(7,2)
    assert F(10)<ell**2<F(16) and d2(*path_points)==16
    atomic_checks=0
    for xi in rational_xis:
        rows=[[xi,xi,1-xi,0],[xi,xi,0,1-xi]]
        for i,row in enumerate(rows):
            in_ball=[j for j,ds in enumerate(distances[i]) if ds<ell**2]
            value=sum((row[j] for j in in_ball),F(0))/len(in_ball)
            assert len(in_ball)==4 and value==(1+xi)/4
            atomic_checks += 1

    # Loss-scaled gate maps remain contractions and admit the stated isometries.
    max_isometry_error=max_rate_error=max_overlap_error=0.0
    instrument_cases=0
    for n in range(17):
        xi_float=n/16
        eta_float=(1+xi_float)/4
        rho_float=2*xi_float/(1+xi_float)
        for m in range(13):
            theta=-math.pi+m*math.pi/6
            a1,a2=finite_detector_vectors(xi_float,theta)
            for efficiency in (1.0,0.9,0.2,0.001):
                qg=[]
                for gate in ((1,1,1,1),(1,1,0,0)):
                    vectors=[[math.sqrt(efficiency)*a[k]*gate[k] for k in range(4)] for a in (a1,a2)]
                    ideal_eta=eta_float if gate[2] else xi_float/2
                    expected_rate=efficiency*ideal_eta
                    qg.append(inner(vectors[0],vectors[0]).real)
                    for vector in vectors:
                        max_rate_error=max(max_rate_error,abs(inner(vector,vector)-expected_rate))
                    j1=vectors[0]+[0j]*4+[math.sqrt(1-expected_rate),0j]
                    j2=[0j]*4+vectors[1]+[0j,math.sqrt(1-expected_rate)]
                    max_isometry_error=max(max_isometry_error,abs(inner(j1,j1)-1),abs(inner(j2,j2)-1),abs(inner(j1,j2)))
                    instrument_cases += 1
                computed_rho=qg[1]/qg[0]
                original_gamma=inner(a2,a1)/eta_float
                relative_gamma=computed_rho*(1+cmath.exp(-1j*theta))/2
                max_overlap_error=max(max_overlap_error,abs(relative_gamma-original_gamma),abs(computed_rho-rho_float))
    assert max(max_isometry_error,max_rate_error,max_overlap_error)<1e-12

    # Verify the relative Lipschitz estimate over deterministic independent pairs.
    stability_cases=0
    for n in range(1001):
        r=(17*n%1001)/1000
        rt=(83*n%1001)/1000
        th=-math.pi+2*math.pi*((37*n%1001)/1000)
        tht=-math.pi+2*math.pi*((71*n%1001)/1000)
        phase_distance=abs((th-tht+math.pi)%(2*math.pi)-math.pi)
        g=r*(1+cmath.exp(-1j*th))/2
        gt=rt*(1+cmath.exp(-1j*tht))/2
        bound=abs(r-rt)+phase_distance/2
        assert abs(g-gt)<=bound+1e-14
        assert abs(abs(g)-abs(gt))<=bound+1e-14
        stability_cases+=1

    # Non-scalar loss on the two shared alternatives invalidates phase cancellation.
    a1,a2=finite_detector_vectors(0.5,math.pi)
    effect=[1.0,0.25,1.0,1.0]
    a1d=[a*math.sqrt(t) for a,t in zip(a1,effect)]
    a2d=[a*math.sqrt(t) for a,t in zip(a2,effect)]
    rate1,rate2=inner(a1d,a1d).real,inner(a2d,a2d).real
    v_filtered=abs(inner(a2d,a1d))/(0.5*rate1+0.5*rate2)
    assert abs(v_filtered-1/3)<1e-14
    assert abs(inner(a2,a1))<1e-14
    # A broader model can match the ratio 2/3 and nevertheless have V=9/10.
    generic1=[math.sqrt(2/3),0j,math.sqrt(1/3),0j]
    generic2=[math.sqrt(2/3),0j,0.7*math.sqrt(1/3),math.sqrt(0.51/3)]
    assert abs(inner(generic1,generic1)-1)<1e-14
    assert abs(inner(generic2,generic2)-1)<1e-14
    assert abs(inner(generic2,generic1)-0.9)<1e-14
    assert abs(sum(abs(t)**2 for t in generic1[:2])-2/3)<1e-14
    assert abs(sum(abs(t)**2 for t in generic2[:2])-2/3)<1e-14

    return {"status":"All second-amendment checks passed",
        "common_efficiency_exact_cases":common_cases,
        "calibrated_efficiency_mismatch_exact_cases":mismatch_cases,
        "uncertainty_interval_containment_cases":interval_cases,
        "explicit_domain_tests":len(domain_tests),"invalid_input_guards":len(invalid_cases),
        "atomic_metric_ball_cases":atomic_checks,
        "atomic_squared_distances":[[str(x) for x in row] for row in distances],
        "loss_completed_instrument_cases":instrument_cases,
        "relative_stability_cases":stability_cases,
        "max_loss_isometry_error":max_isometry_error,
        "max_loss_rate_error":max_rate_error,"max_relative_overlap_error":max_overlap_error,
        "ten_percent_loss":{"intrinsic_eta":str(eta),"efficiency":str(e),
            "all_rate":str(q_all),"shared_rate":str(q_shared),"relative_accessibility":str(rho),
            "naive_xi":str(naive_xi),"naive_visibility":str(naive_visibility),
            "exact_visibility_bias":str(rho-naive_visibility),"original_bound":str(8*(eta-q_all)),
            "corrected_xi":str(rho/(2-rho)),"corrected_eta":str(1/(2*(2-rho)))},
        "relative_held_out_test":{"predicted_visibility":str(rho),"hypothetical_visibility":str(foil),
            "gap":str(foil-rho),"predeclared_example_tolerance":str(tolerance),"rejected":True},
        "nonscalar_filter_counterexample":{"xi":"1/2","theta":"pi", "effect_diagonal":["1","1/4","1","1"],
             "ideal_visibility":"0","filtered_visibility_exact":"1/3","calculated":v_filtered},
        "scope":"Conditional, synthetic instrument and measurement checks; no experimental-anchor implementation claimed"}


def generic_block_vectors(eta_a: float, rho: float, omega: float, theta: float,
                          private_overlap: float = 0.0
                          ) -> tuple[list[complex], list[complex]]:
    """A comparison class parametrized by norms, not by the stipulated kernel.

    Shared amplitudes agree at theta=0; theta acts on the second shared mode.
    Zero private_overlap implements the block lemma. Nonzero values supply
    an explicit foil to private orthogonality, without invalidating the instrument.
    """
    values = (eta_a, rho, omega, theta, private_overlap)
    if any(not math.isfinite(x) for x in values):
        raise ValueError("All parameters must be finite")
    if not (0 < eta_a <= 1 and 0 <= rho <= 1 and 0 <= omega <= 1
            and 0 <= private_overlap <= 1):
        raise ValueError("Invalid transmission, shared fraction, split, or overlap")
    ns, np = eta_a*rho, eta_a*(1-rho)
    c1 = math.sqrt(ns*omega)*cmath.exp(0.23j)
    c2 = math.sqrt(ns*(1-omega))*cmath.exp(-0.41j)
    # A fixed common phase on private coordinates does not affect their overlap.
    r = math.sqrt(np)*cmath.exp(0.37j)
    a1 = [c1,c2,r,0j]
    a2 = [c1,c2*cmath.exp(1j*theta),private_overlap*r,
          math.sqrt(max(0.0,1-private_overlap**2))*r]
    return a1,a2


def accessibility_residual_interval(eta_a: tuple[F,F], eta_c: tuple[F,F]
                                   ) -> tuple[F,F]:
    """Conservative closure residual from simultaneously valid marginal bounds.

    No independence is assumed. A joint observation model can give a smaller set.
    The observable test is 2*eta_A - eta_C - 1/2, not substitution of the kernel's
    own inferred absolute norm into both sides of its prediction.
    """
    al,ah = _bounds(eta_a,"all-gate accessibility",unit=True)
    cl,ch = _bounds(eta_c,"shared-gate accessibility",unit=True)
    return 2*al-ch-F(1,2), 2*ah-cl-F(1,2)


def verify_third_amendment() -> dict:
    """Verify hypothesis-specific attribution and accessibility-only closure.

    These calculations add no new observations. The old zero-phase rejection
    remains numerically true, but its hypothesis class is the broader shared /
    orthogonal-private class, not specifically the named incidence kernel.
    """
    transmissions = [F(1,10), F(1,4), F(3,10), F(3,8), F(1,2), F(9,10), F(1)]
    rhos = [F(n,8) for n in range(9)]
    splits = [F(n,6) for n in range(7)]
    quarters = [(F(1),F(0)),(F(0),F(-1)),(F(-1),F(0)),(F(0),F(1))]
    exact_zero = exact_split = interval_cases = 0
    for eta in transmissions:
        for rho in rhos:
            ns,np = eta*rho,eta*(1-rho)
            assert ns+np==eta and ns/eta==rho
            # This equality uses no xi and no 1/4 weight.
            exact_zero += 1
            err = F(1,1000)
            ai=(max(F(0),eta-err),min(F(1),eta+err))
            ci=(max(F(0),ns-err),min(F(1),ns+err))
            lo,hi=accessibility_residual_interval(ai,ci)
            assert lo <= 2*eta-ns-F(1,2) <= hi
            interval_cases += 1
            for w in splits:
                for re,im in quarters:
                    gr=rho*(w+(1-w)*re)
                    gi=rho*(1-w)*im
                    v2=rho*rho*(1-2*w*(1-w)*(1-re))
                    assert gr*gr+gi*gi==v2
                    if re==1:
                        assert gr==rho and gi==0
                    if re==-1:
                        assert gr==rho*(2*w-1)
                        if rho>0: assert (gr==0)==(w==F(1,2))
                    gr_swap=rho*((1-w)+w*re)
                    gi_swap=rho*w*im
                    assert gr_swap*gr_swap+gi_swap*gi_swap==v2
                    exact_split += 1

    complex_cases=scale_cases=0
    max_norm=max_overlap=max_scale=max_completion=0.0
    phases=(-math.pi,-0.6,0.0,0.17,math.pi/2,math.pi,2.7)
    for etaq in transmissions:
        eta=float(etaq)
        for rhoq in rhos:
            rho=float(rhoq)
            for wq in splits:
                w=float(wq)
                for th in phases:
                    a1,a2=generic_block_vectors(eta,rho,w,th)
                    norm1,norm2=inner(a1,a1).real,inner(a2,a2).real
                    shared1=inner(a1[:2],a1[:2]).real
                    shared2=inner(a2[:2],a2[:2]).real
                    gamma=inner(a2,a1)/math.sqrt(norm1*norm2)
                    expected=rho*(w+(1-w)*cmath.exp(-1j*th))
                    max_norm=max(max_norm,abs(norm1-eta),abs(norm2-eta),
                                 abs(shared1-eta*rho),abs(shared2-eta*rho))
                    max_overlap=max(max_overlap,abs(gamma-expected))
                    v2=rho*rho*(1-4*w*(1-w)*math.sin(th/2)**2)
                    assert abs(abs(gamma)**2-v2)<1e-13
                    if th==0.0: assert abs(gamma-shared1/norm1)<1e-13
                    j1=a1+[0j]*4+[math.sqrt(1-eta),0j]
                    j2=[0j]*4+a2+[0j,math.sqrt(1-eta)]
                    max_completion=max(max_completion,abs(inner(j1,j1)-1),
                                       abs(inner(j2,j2)-1),abs(inner(j1,j2)))
                    for t in (0.2,0.8,1.0,1.5):
                        if t*eta>1: continue
                        b1,b2=[[math.sqrt(t)*z for z in a] for a in (a1,a2)]
                        bn=inner(b1,b1).real
                        bg=inner(b2,b1)/math.sqrt(bn*inner(b2,b2).real)
                        br=inner(b1[:2],b1[:2]).real/bn
                        max_scale=max(max_scale,abs(bg-gamma),abs(br-rho),abs(bn-t*eta))
                        scale_cases += 1
                    complex_cases += 1
    assert max(max_norm,max_overlap,max_scale,max_completion)<1e-12

    # The same nominal rho and absolute norm can fail the block identity when
    # private components have overlap 7/10 instead of being orthogonal.
    a1,a2=generic_block_vectors(3/8,2/3,0.5,0.0,private_overlap=0.7)
    eta_calc=inner(a1,a1).real
    rho_calc=inner(a1[:2],a1[:2]).real/eta_calc
    gamma_calc=inner(a2,a1)/eta_calc
    assert abs(eta_calc-3/8)<1e-14 and abs(rho_calc-2/3)<1e-14
    assert abs(gamma_calc-0.9)<1e-14
    # Exact rational counterpart:
    ns,np,ov=F(1,4),F(1,8),F(7,80)
    assert (ns+ov)/(ns+np)==F(9,10)
    assert ov/np==F(7,10)
    assert 2*(ns+np)-ns-F(1,2)==0

    # Equal-split foil: both accessibility checks and V(0)=rho pass, V(pi) fails.
    a1,a2=generic_block_vectors(3/8,2/3,3/4,math.pi)
    unequal_split_v=abs(inner(a2,a1))/(3/8)
    assert abs(unequal_split_v-1/3)<1e-14

    # Absolute-normalization foil: the entire relative phase curve passes.
    norm_foil=F(3,10)
    c_foil=F(1,5)
    ratio_foil=c_foil/norm_foil
    residual=2*norm_foil-c_foil-F(1,2)
    assert ratio_foil==F(2,3) and residual==F(-1,10)
    assert 1/(2*(2-ratio_foil))==F(3,8)!=norm_foil
    assert 2-1/(2*norm_foil)==F(1,3)!=ratio_foil
    for th in phases:
        aa=[1/math.sqrt(10),1/math.sqrt(10),1/math.sqrt(10),0j]
        bb=[1/math.sqrt(10),cmath.exp(1j*th)/math.sqrt(10),0j,1/math.sqrt(10)]
        assert abs(inner(bb,aa)/float(norm_foil)-(1+cmath.exp(-1j*th))/3)<1e-13
    foil_interval=accessibility_residual_interval((F(299,1000),F(301,1000)),
                                                  (F(199,1000),F(201,1000)))
    assert foil_interval==(F(-103,1000),F(-97,1000))
    assert not foil_interval[0]<=0<=foil_interval[1]
    valid_interval=accessibility_residual_interval((F(3,8)-F(1,1000),F(3,8)+F(1,1000)),
                                                  (F(1,4)-F(1,1000),F(1,4)+F(1,1000)))
    assert valid_interval==(F(-3,1000),F(3,1000))

    # Exact count-space residual: the shared q_A in rho is never treated as
    # independent; test the algebra on the original rates and calibrated losses.
    count_cases=0
    for xi in [F(n,20) for n in range(21)]:
        eta,c=(1+xi)/4,xi/2
        rho=c/eta
        assert 2*eta-c==F(1,2)
        assert eta==1/(2*(2-rho))
        for ea in (F(1,10),F(1,2),F(9,10),F(1)):
            for ec in (F(1,5),F(2,5),F(4,5),F(1)):
                qa,qc,zeta=ea*eta,ec*c,ec/ea
                assert (2*qa-qc/zeta)/ea-F(1,2)==0
                assert qc/(zeta*qa)==rho
                count_cases += 1

    # Source weights cannot be recovered from ratio curves under common scaling.
    scale_exact=0
    for eta in transmissions:
        for rho in rhos:
            for tq in (F(1,5),F(4,5),F(1),F(3,2)):
                if tq*eta>1: continue
                assert tq*eta*rho/(tq*eta)==rho
                scale_exact += 1
    invalid=0
    for args in [(0.,0.5,0.5,0.),(1.1,0.5,0.5,0.),(0.5,-0.1,0.5,0.),
                 (0.5,0.5,1.1,0.),(0.5,0.5,0.5,float('nan'))]:
        try: generic_block_vectors(*args)
        except ValueError: invalid += 1
        else: raise AssertionError('Missing input guard')

    return {
        'status':'All third-amendment checks passed',
        'generic_zero_phase_exact_cases':exact_zero,
        'general_split_exact_quarter_turn_cases':exact_split,
        'generic_complex_vector_cases':complex_cases,
        'rescaling_complex_cases':scale_cases,
        'rescaling_exact_cases':scale_exact,
        'accessibility_residual_interval_containment_cases':interval_cases,
        'exact_corrected_count_closure_cases':count_cases,
        'invalid_input_guards':invalid,
        'max_norm_error':max_norm,'max_overlap_error':max_overlap,
        'max_rescaling_error':max_scale,'max_instrument_completion_error':max_completion,
        'reattribution':{
            'relative_held_out_test_from_second_amendment':
                'Numerical rejection unchanged; it follows already from the identical-shared/orthogonal-private block class with balanced preparation and calibrated gates/losses, not specifically the chosen kernel.',
            'zero_phase':'Generic block and measurement assumptions',
            'phase_scan':'Equal shared split within the declared block and phase-action class',
            'absolute_route':'Joint absolute-normalization consequence, not separate unique identification of the lift and reference weight',
            'accessibility_only':'2 eta_A - eta_C = 1/2; independent absolute efficiency required; no coherence data and no substitution of model-inferred eta_A',
            'independence':'Additional diagnostic, not statistically independent evidence; shared counts and calibration inputs must be propagated jointly.'},
        'block_identity_foil':{'eta_A':'3/8','eta_C':'1/4','rho':'2/3',
            'private_inner_product':'7/80','private_normalized_overlap':'7/10','gamma_zero':'9/10',
            'violates':'Private orthogonality','valid_broader_quantum_instrument':True},
        'unequal_split_foil':{'eta_A':'3/8','eta_C':'1/4','rho':'2/3','omega':'3/4',
            'accessibility_residual':'0','visibility_pi':'1/3','kernel_visibility_pi':'0',
            'violates':'Equal shared split only, among these three structural tests'},
        'absolute_normalization_foil':{'eta_A':'3/10','eta_C':'1/5','rho':'2/3',
            'gamma':'(1+exp(-i theta))/3','accessibility_residual':'-1/10',
            'kernel_eta_A':'3/8','route_A_visibility_zero':'1/3','actual_visibility_zero':'2/3',
            'residual_interval':[str(x) for x in foil_interval],
            'passes':'Block identity and entire equal-split relative phase curve',
            'violates':'Absolute-normalization closure'},
        'valid_closure_interval':[str(x) for x in valid_interval],
        'coefficient_convention':'b(eta)=1-1/(4 eta); full gamma(0)=2 b(eta)',
        'scope':'Synthetic checks and conditional proofs; no experimental verification or unique ECL identification.'}


def main() -> None:
    edges_k = [(i, j) for i in range(3) for j in range(3, 6)]
    edges_p = [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]
    specs = {"K33": [0,3,3,3,3,6], "prism": [0,2,3,3,5,5]}
    graphs = {}
    for name, edges in [("K33", edges_k), ("prism", edges_p)]:
        l = laplacian(edges)
        assert all(l[i][i] == 3 for i in range(6))
        # Equality of two degree-six monic polynomials at seven distinct points.
        for x in range(7):
            char_value = determinant([[x*(i == j)-l[i][j] for j in range(6)] for i in range(6)])
            assert char_value == math.prod(x-e for e in specs[name])
        triangles = sum((-l[i][j])*(-l[j][k])*(-l[k][i])
                        for i in range(6) for j in range(i+1,6) for k in range(j+1,6))
        response = {}
        for kappa in (1, 2, 4):
            a = [[l[i][j]+kappa*(i == j) for j in range(6)] for i in range(6)]
            b = [1,-1,0,0,0,0]
            x = solve(a,b)
            chi = sum(F(b[i])*x[i] for i in range(6))/2
            expected = F(1,kappa+3) if name == "K33" else (F(1,kappa+3)+F(1,kappa+5))/2
            assert chi == expected
            response[str(kappa)] = {"exact": str(chi), "decimal": float(chi)}
        graphs[name] = {"edges": edges, "edge_density": 0.6,
                        "laplacian": l, "spectrum": specs[name],
                        "triangles": triangles, "susceptibility": response}
    epsilon = 0.01
    times = [1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0]
    heat = []
    max_delta = 0.0
    for t in times:
        value = spectral_dimension(t, epsilon)
        delta = abs(value-spectral_dimension(t,epsilon,90.0))
        max_delta = max(max_delta, delta)
        assert delta < 1e-12
        heat.append({"auxiliary_scale": t, "epsilon": epsilon, "d_s": value})
    multifractal = [{"p": p, "D0": generalized_dimension(p,0),
                    "D1": generalized_dimension(p,1), "D2": generalized_dimension(p,2)}
                   for p in (0.5,0.65,0.8,0.95)]
    q0, q1 = [[-12,-5],[-5,-4]], [[-12,1],[1,-2]]
    assert determinant(q0) == determinant(q1) == 23
    assert q1[1][1] == -2
    # Finite enumeration supplements, but does not replace, the exact proof in the article.
    assert all(-12*m*m-10*m*n-4*n*n != -2
               for m in range(-50,51) for n in range(-50,51))
    out = {"status": "All exact and numerical checks passed",
           "data_type": "synthetic mathematical examples, not experimental data",
           "graphs": graphs, "spectral_example": heat,
           "theta_cutoff_doubling_max_difference": max_delta,
           "multifractal_example": multifractal,
           "fourball_comparison": {"Q0": q0,"Q1": q1,"determinants": [23,23]},
           "synthetic_response_gap_kappa1": str(F(1,4)-F(5,24)),
           "synthetic_response_gap_kappa2": str(F(1,5)-F(6,35))}
    out["amendment"] = verify_amendments(graphs)
    out["second_amendment"] = verify_second_amendment()
    out["third_amendment"] = verify_third_amendment()
    destination = Path(__file__).resolve().parent / "verified_examples.json"
    destination.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
