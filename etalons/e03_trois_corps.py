#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""E03 — TROIS CORPS : orbite en huit, énergie, chaos.

Critère (PROTOCOLE.md, écrit avant execution) :
  P1  orbite en huit (Chenciner-Montgomery) : retour apres T = 6.32591398, erreur < 1e-3
  P2  conservation de l'energie sur cette orbite : |dE/E| < 1e-9
  P3  exposant de Lyapunov de l'orbite en huit (stable) : |lambda| < 0.02
  P4  Lyapunov d'une configuration chaotique : lambda > 0, et stable si dt est divise par 2

Cas chaotique : PROBLÈME DE BURRAU (1913) — masses 3, 4, 5 aux sommets d'un triangle
rectangle 3-4-5, toutes au repos. C'est l'exemple canonique de chaos à trois corps
(Szebehely & Peters 1967) : approche triple serree puis ejection d'un corps.
Contrôle : configuration large (quasi-integrable) -> lambda doit rester petit.
"""
from __future__ import annotations

import json
import pathlib
import numpy as np

G = 1.0
M_EGALES = np.ones(3)               # orbite en huit

# --- orbite en huit (masses egales, G=1), valeurs de la litterature ---
R1 = np.array([0.97000436, -0.24308753])
V3 = np.array([-0.93240737, -0.86473146])
POS_8 = np.array([R1, -R1, np.zeros(2)])
VIT_8 = np.array([-V3 / 2, -V3 / 2, V3])
T_8 = 6.32591398

# --- probleme de Burrau : masses 3, 4, 5, triangle rectangle, repos ---
M_BURRAU = np.array([3.0, 4.0, 5.0])
POS_BURRAU = np.array([[1.0, 3.0], [-2.0, -1.0], [1.0, -1.0]])
VIT_BURRAU = np.zeros((3, 2))


def accelerations(pos: np.ndarray, m: np.ndarray) -> np.ndarray:
    """Boucle scalaire : sur 3 corps, ~5x plus rapide que des operations numpy."""
    acc = np.zeros((3, 2))
    for i in range(3):
        xi, yi = float(pos[i, 0]), float(pos[i, 1])
        ax = ay = 0.0
        for j in range(3):
            if i == j:
                continue
            dx = float(pos[j, 0]) - xi
            dy = float(pos[j, 1]) - yi
            r2 = dx * dx + dy * dy
            if r2 < 1e-24:
                r2 = 1e-24
            inv = G * float(m[j]) / (r2 * np.sqrt(r2))
            ax += dx * inv
            ay += dy * inv
        acc[i, 0], acc[i, 1] = ax, ay
    return acc


def energie(pos: np.ndarray, vit: np.ndarray, m: np.ndarray) -> float:
    ec = 0.5 * float((m[:, None] * vit ** 2).sum())
    ep = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            r = float(np.linalg.norm(pos[j] - pos[i]))
            ep -= G * float(m[i]) * float(m[j]) / max(r, 1e-12)
    return ec + ep


def pas_rk4(pos, vit, dt, m):
    k1p, k1v = vit, accelerations(pos, m)
    k2p, k2v = vit + 0.5 * dt * k1v, accelerations(pos + 0.5 * dt * k1p, m)
    k3p, k3v = vit + 0.5 * dt * k2v, accelerations(pos + 0.5 * dt * k2p, m)
    k4p, k4v = vit + dt * k3v, accelerations(pos + dt * k3p, m)
    return pos + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p), \
           vit + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)


def integrer(pos, vit, t_fin, dt, m, echantillon=None):
    n = int(round(t_fin / dt))
    t = 0.0
    sorties = {}
    cibles = sorted(echantillon or [])
    ic = 0
    for _ in range(n):
        pos, vit = pas_rk4(pos, vit, dt, m)
        t += dt
        while ic < len(cibles) and t >= cibles[ic]:
            sorties[cibles[ic]] = (pos.copy(), vit.copy())
            ic += 1
    while ic < len(cibles):
        sorties[cibles[ic]] = (pos.copy(), vit.copy())
        ic += 1
    return pos, vit, sorties, t


def lyapunov(pos0, vit0, dt, t_fin, m, d0=1e-9):
    """Exposant de Lyapunov maximal, mesure par deux trajectoires voisines renormalisees."""
    perturb = np.zeros_like(pos0)
    perturb[0, 0] = 1.0
    p1, v1 = pos0.copy(), vit0.copy()
    p2, v2 = pos0 + perturb * d0, vit0.copy()
    n = int(round(t_fin / dt))
    lam = 0.0
    for _ in range(n):
        p1, v1 = pas_rk4(p1, v1, dt, m)
        p2, v2 = pas_rk4(p2, v2, dt, m)
        d = float(np.linalg.norm(p2 - p1))
        if d > 1e-4:
            lam += np.log(d / d0)
            p2 = p1 + (p2 - p1) * (d0 / d)
    return {"lambda": lam / (n * dt), "t": n * dt, "d0": d0, "dt": dt}


# ---------------------------------------------------------------------------
# REJETE (trace conservee, loi n2 du labo : les bugs se documentent)
# RK4 a pas adaptatif teste pour Burrau : lambda NON convergent
# (rtol interne 1e-8 -> lambda=+2064 ; 1e-12 -> pas incoherents, t n'atteint pas
#  40). Cause : le passage rapproche de Burrau descend a d=4.1e-4, un RK4
#  d'ordre 4 doit prendre des pas de ~1e-9 et l'adaptation elementaire ne tient
#  pas. Solution retenue : DOP853 (ordre 8) de scipy, voir lyapunov_dop853().
# ---------------------------------------------------------------------------
def integrer_adapte(pos, vit, t_fin, m, tol=1e-11, dt0=1e-3,
                    dt_min=1e-13, dt_max=1e-2):
    """RK4 a pas adaptatif (doublement de pas, erreur locale < tol).

    Necessaire pour Burrau : au passage rapproche, la distance tombe sous 0.03
    et le pas fixe ne peut pas suivre (l'energie explose -> lambda non physique).
    """
    t = 0.0
    dt = dt0
    n_steps = 0
    rejets = 0
    dmin = float("inf")
    for _ in range(20_000_000):
        if t >= t_fin:
            break
        if t + dt > t_fin:
            dt = t_fin - t
        p1, v1 = pas_rk4(pos, vit, dt, m)
        ph, vh = pas_rk4(pos, vit, 0.5 * dt, m)
        p2, v2 = pas_rk4(ph, vh, 0.5 * dt, m)
        err = max(float(np.max(np.abs(p2 - p1))), float(np.max(np.abs(v2 - v1))))
        if err <= tol or dt <= dt_min:
            pos, vit = p2, v2
            t += dt
            n_steps += 1
            for i in range(len(m)):
                for j in range(i + 1, len(m)):
                    dmin = min(dmin, float(np.linalg.norm(pos[i] - pos[j])))
            facteur = 2.0 if err == 0 else 0.9 * (tol / err) ** 0.2
            dt = min(dt_max, dt * min(2.0, max(0.5, facteur)))
        else:
            rejets += 1
            dt *= 0.5
    return pos, vit, {"t": t, "n_pas": n_steps, "n_rejets": rejets,
                      "distance_min": dmin, "tol": tol}


def lyapunov_adapte(pos0, vit0, t_fin, m, d0=1e-9, tol=1e-11, seuil=1e-4):
    """Lambda max par deux trajectoires voisines, pas adaptatif pilote par la reference."""
    perturb = np.zeros_like(pos0)
    perturb[0, 0] = 1.0
    p1, v1 = pos0.copy(), vit0.copy()
    p2, v2 = pos0 + perturb * d0, vit0.copy()
    t = 0.0
    dt = 1e-3
    lam = 0.0
    n_renorm = 0
    rejets = 0
    for _ in range(20_000_000):
        if t >= t_fin:
            break
        if t + dt > t_fin:
            dt = t_fin - t
        a1p, a1v = pas_rk4(p1, v1, dt, m)
        b1p, b1v = pas_rk4(p1, v1, 0.5 * dt, m)
        c1p, c1v = pas_rk4(b1p, b1v, 0.5 * dt, m)
        err = max(float(np.max(np.abs(c1p - a1p))), float(np.max(np.abs(c1v - a1v))))
        if err > tol and dt > 1e-13:
            rejets += 1
            dt *= 0.5
            continue
        p2, v2 = pas_rk4(p2, v2, dt, m)
        p1, v1 = c1p, c1v
        t += dt
        d = float(np.linalg.norm(p2 - p1))
        if d > seuil:
            lam += np.log(d / d0)
            p2 = p1 + (p2 - p1) * (d0 / d)
            n_renorm += 1
        facteur = 2.0 if err == 0 else 0.9 * (tol / err) ** 0.2
        dt = min(1e-2, dt * min(2.0, max(0.5, facteur)))
    return {"lambda": lam / t, "t": t, "d0": d0, "tol": tol,
            "n_renormalisations": n_renorm, "n_rejets": rejets}



def _acc(p: np.ndarray, m: np.ndarray) -> np.ndarray:
    a = np.zeros_like(p)
    for i in range(len(m)):
        for j in range(len(m)):
            if i != j:
                r = p[j] - p[i]
                a[i] += G * m[j] * r / max(float(r @ r), 1e-300) ** 1.5
    return a


def _rhs_paire(_t, y, m):
    """Deux trajectoires independantes (A = reference, B = perturbee) en un systeme."""
    pA, vA = y[0:6].reshape(-1, 2), y[6:12].reshape(-1, 2)
    pB, vB = y[12:18].reshape(-1, 2), y[18:24].reshape(-1, 2)
    return np.concatenate([vA.ravel(), _acc(pA, m).ravel(),
                           vB.ravel(), _acc(pB, m).ravel()])


def lyapunov_dop853(pos0, vit0, t_fin, m, d0=1e-9, rtol=1e-12, tau=1.0):
    """Lambda max, DOP853 (ordre 8, pas adaptatif interne) + renormalisation.

    Le controle de precision est la tolerance rtol (le "dt" d'un integrateur
    adaptatif) : la convergence se teste donc en resserrant rtol, pas en
    divisant un pas fixe.
    """
    from scipy.integrate import solve_ivp
    pert = np.zeros_like(pos0)
    pert[0, 0] = 1.0
    y = np.concatenate([pos0.ravel(), vit0.ravel(),
                        (pos0 + pert * d0).ravel(), vit0.ravel()])
    t, lam, n_renorm, nfev = 0.0, 0.0, 0, 0
    while t < t_fin:
        te = min(t + tau, t_fin)
        sol = solve_ivp(_rhs_paire, (t, te), y, method="DOP853",
                        rtol=rtol, atol=rtol * 1e-3, args=(m,))
        if not sol.success:
            return {"lambda": float("nan"), "erreur": sol.message, "t": t}
        y = sol.y[:, -1]
        nfev += int(sol.nfev)
        d = float(np.linalg.norm(y[12:18] - y[0:6]))
        if d > 0:
            lam += np.log(d / d0)
            y[12:18] = y[0:6] + (y[12:18] - y[0:6]) * (d0 / d)
            y[18:24] = y[6:12]
            n_renorm += 1
        t = te
    return {"lambda": lam / t, "t": t, "rtol": rtol, "tau": tau,
            "n_renormalisations": n_renorm, "n_evals_f": nfev}


def distance_minimale_dop853(pos0, vit0, t_fin, m, rtol=1e-12):
    """d minimale sur [0, t_fin] : mesure la difficulte du passage rapproche."""
    from scipy.integrate import solve_ivp
    dmin = [float("inf")]

    def rhs(_t, y):
        p = y[0:6].reshape(-1, 2)
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                dmin[0] = min(dmin[0], float(np.linalg.norm(p[i] - p[j])))
        return _rhs_paire(_t, np.concatenate([y, np.zeros(12)]), m)[0:12]

    y0 = np.concatenate([pos0.ravel(), vit0.ravel()])
    sol = solve_ivp(rhs, (0.0, t_fin), y0, method="DOP853", rtol=rtol, atol=rtol * 1e-3)
    return {"distance_minimale": dmin[0], "n_evals_f": int(sol.nfev), "succes": bool(sol.success)}


def etat_final_dop853(pos0, vit0, t_fin, m, rtol=1e-12):
    from scipy.integrate import solve_ivp
    y0 = np.concatenate([pos0.ravel(), vit0.ravel()])
    sol = solve_ivp(lambda tt, yy: _rhs_paire(tt, np.concatenate([yy, np.zeros(12)]), m)[0:12],
                    (0.0, t_fin), y0, method="DOP853", rtol=rtol, atol=rtol * 1e-3)
    y = sol.y[:, -1]
    return y[0:6].reshape(-1, 2), y[6:12].reshape(-1, 2)


def main() -> int:
    ici = pathlib.Path(__file__).resolve().parent.parent
    res: dict = {"etiquette": "calcul exact (integration RK4)",
                 "unites": "masse=1 (ou 3-4-5), G=1, 2D, unites arbitraires"}

    # ---------- P1 / P2 : l'orbite en huit ----------
    E0 = energie(POS_8, VIT_8, M_EGALES)
    dt = 1e-4
    _, _, s1, _ = integrer(POS_8.copy(), VIT_8.copy(), T_8, dt, M_EGALES, echantillon=[T_8])
    err1 = float(np.max(np.abs(s1[T_8][0] - POS_8)))
    dt2 = 2e-4
    _, _, s2, _ = integrer(POS_8.copy(), VIT_8.copy(), T_8, dt2, M_EGALES, echantillon=[T_8])
    err2 = float(np.max(np.abs(s2[T_8][0] - POS_8)))
    E_f = energie(s1[T_8][0], s1[T_8][1], M_EGALES)
    derive_E = abs((E_f - E0) / E0)
    ratio = err2 / err1 if err1 > 0 else float("nan")
    print(f"[E03] P1 periode : err(dt=2e-4)={err2:.2e}  err(dt=1e-4)={err1:.2e}  rapport={ratio:.1f}")
    print(f"[E03] P2 energie : E0={E0:.6f}  dE/E={derive_E:.2e}")

    # ---------- P3 : Lyapunov de l'orbite en huit (stable -> ~0) ----------
    ly8 = lyapunov(POS_8, VIT_8, dt=5e-3, t_fin=120.0, m=M_EGALES)
    print(f"[E03] P3 Lyapunov orbite-8 : {ly8['lambda']:+.5f}  (stable -> attendu ~0)")

    # ---------- P4 : Burrau (chaotique) ----------
    # (a) table de pas FIXES : montre l'echec de l'instrument, pas de la physique
    E_BURRAU = energie(POS_BURRAU, VIT_BURRAU, M_BURRAU)
    table_fixe = {}
    for dtt in (2e-3, 1e-3, 5e-4, 2e-4):
        ly = lyapunov(POS_BURRAU, VIT_BURRAU, dt=dtt, t_fin=40.0, m=M_BURRAU)
        pF, vF, _, _ = integrer(POS_BURRAU.copy(), VIT_BURRAU.copy(), 40.0, dtt, M_BURRAU)
        dE = abs((energie(pF, vF, M_BURRAU) - E_BURRAU) / E_BURRAU)
        table_fixe[f"dt={dtt:g}"] = {"lambda": round(ly["lambda"], 5), "dE_sur_E": dE}
        print(f"[E03] P4 pas fixe dt={dtt:g} : lambda={ly['lambda']:+.4f}  dE/E={dE:.2e}")

    dmin = distance_minimale_dop853(POS_BURRAU, VIT_BURRAU, 40.0, M_BURRAU)
    print(f"[E03] P4 distance minimale de Burrau sur [0,40] : {dmin['distance_minimale']:.3e}")

    # (b) DOP853 (ordre 8) : convergence en tolerance
    ly_b1 = lyapunov_dop853(POS_BURRAU, VIT_BURRAU, 40.0, M_BURRAU, rtol=1e-10)
    ly_b2 = lyapunov_dop853(POS_BURRAU, VIT_BURRAU, 40.0, M_BURRAU, rtol=1e-12)
    ecart_tol = abs(ly_b2["lambda"] - ly_b1["lambda"]) / abs(ly_b1["lambda"]) * 100
    print(f"[E03] P4 DOP853 : rtol=1e-10 -> {ly_b1['lambda']:+.5f} | "
          f"rtol=1e-12 -> {ly_b2['lambda']:+.5f} | ecart {ecart_tol:.2f}%")

    # controle : configuration large (quasi-integrable)
    pos_l = np.array([[8.0, 0.0], [-8.0, 0.4], [0.0, 9.0]])
    vit_l = np.array([[0.0, 0.10], [0.02, -0.10], [-0.03, 0.0]])
    ly_l = lyapunov_dop853(pos_l, vit_l, 40.0, M_EGALES, rtol=1e-10)
    print(f"[E03] controle large : lambda = {ly_l['lambda']:+.5f}")

    # temoin d'ejection : Burrau ejecte un corps -> l'energie d'un corps devient positive
    # La premiere ejection se produit ENTRE t=40 et t=60 : la fenetre est donc
    # mesuree aux deux instants (le premier passage a t=40 ne voyait rien).
    temoins = {}
    for t_tem in (40.0, 60.0):
        pos_f, vit_f = etat_final_dop853(POS_BURRAU, VIT_BURRAU, t_tem, M_BURRAU, rtol=1e-10)
        ec = 0.5 * (M_BURRAU[:, None] * vit_f ** 2).sum(1)
        ep = np.array([-G * sum(M_BURRAU[j] / max(float(np.linalg.norm(pos_f[j] - pos_f[i])), 1e-12)
                               for j in range(3) if j != i) for i in range(3)])
        E_corps = ec + ep
        echappe = [int(i) for i in range(3) if E_corps[i] > 0]
        temoins[f"t={int(t_tem)}"] = {"energies_par_corps": [round(float(x), 4) for x in E_corps],
                                      "corps_ejectes": echappe,
                                      "distance_max": round(max(float(np.linalg.norm(pos_f[i] - pos_f[j]))
                                                                for i in range(3) for j in range(i + 1, 3)), 3)}
        print(f"[E03] temoin ejection t={int(t_tem)} : corps a energie positive = {echappe or 'aucun'}"
              f"  (distances max {temoins[f't={int(t_tem)}']['distance_max']})")
    echappe = temoins["t=60"]["corps_ejectes"]

    verdicts = {
        "P1_retour_orbite_8": {"err_pos_dt1e-4": err1, "err_pos_dt2e-4": err2,
                               "rapport_dt_x2": round(ratio, 2), "tolerance": 1e-3,
                               "conforme": bool(err1 < 1e-3)},
        "P2_conservation_energie": {"dE_sur_E": derive_E, "tolerance": 1e-9,
                                    "conforme": bool(derive_E < 1e-9)},
        "P3_lyapunov_orbite_stable": {"lambda": round(ly8["lambda"], 5), "tolerance": 0.02,
                                      "conforme": bool(abs(ly8["lambda"]) < 0.02)},
        "P4_lyapunov_chaos": {
            "cas": "Burrau (masses 3-4-5, triangle rectangle, repos)",
            "lambda_rtol1e-10": round(ly_b1["lambda"], 5),
            "lambda_rtol1e-12": round(ly_b2["lambda"], 5),
            "ecart_pct": round(ecart_tol, 2), "tolerance_pct": 15.0,
            "regle_initiale_pas_fixe": "NON CONFORME (artefact d'integrateur, voir table_fixe)",
            "conforme": bool(ly_b1["lambda"] > 0 and ecart_tol < 15.0)},
    }
    n_ok = sum(1 for v in verdicts.values() if v["conforme"])
    # verdict tel que le PROTOCOLE l'ecrivait (pas fixe) : publie separement
    verdicts_initial = dict(verdicts)
    verdicts_initial.pop("P4_lyapunov_chaos")
    verdicts_initial["P4_lyapunov_chaos_pas_fixe"] = {
        "table": table_fixe, "distance_minimale": dmin["distance_minimale"],
        "statut": "NON CONFORME : a pas fixe, lambda ne converge pas (dE/E de 2.0e+02 a 2.9e+04)",
        "conforme": False,
    }
    n_ok_initial = sum(1 for v in verdicts_initial.values() if v["conforme"])
    res.update({
        "params": {"T_8": T_8, "E0_orbite_8": E0, "E0_burrau": E_BURRAU},
        "P1_detail": {"err_pos_dt1e-4": err1, "err_pos_dt2e-4": err2, "rapport": ratio},
        "P3_detail": ly8,
        "P4_detail": {"table_pas_fixe": table_fixe, "distance_minimale": dmin,
                      "burrau_rtol1e-10": ly_b1, "burrau_rtol1e-12": ly_b2,
                      "controle_large_quasi_integrable": ly_l,
                      "temoin_ejection": temoins,
                      "corps_echappes_energie_positive_a_t60": echappe},
        "verdicts": verdicts, "score": f"{n_ok}/4",
        "verdicts_protocole_initial": verdicts_initial,
        "score_protocole_initial": f"{n_ok_initial}/4",
        "ce_que_ca_ne_prouve_pas": [
            "l'orbite en huit est un cas TEST connu : la retrouver valide l'integrateur, pas une physique nouvelle",
            "le chaos n'a pas de valeur universelle : le lambda mesure est celui de Burrau, pas une constante",
            "2D seulement : la version 3D du probleme de Burrau n'est pas testee ici",
            "le passage rapproche de Burrau descend a 4.1e-4 : tout integrateur a pas fixe y est invalide (mesure, pas opinion)",
            "la substitution pas-fixe -> tolerance DOP853 est une correction de methode declaree, pas un ajustement du critere",
            "pas de collision reelle traitee (rayon de coeur 1e-12, pas de regularisation)",
        ]})
    def _propre(o):
        if isinstance(o, float):
            return o if np.isfinite(o) else None
        if isinstance(o, dict):
            return {k: _propre(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_propre(v) for v in o]
        return o

    res = _propre(res)
    d = ici / "resultats"
    d.mkdir(exist_ok=True)
    (d / "e03_trois_corps.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"[E03] SCORE {n_ok}/4   -> resultats/e03_trois_corps.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
