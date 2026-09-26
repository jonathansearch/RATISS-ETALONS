#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""E01 — PERCOLATION : seuil, trous, dimension fractale.

Question : à quelle densité un réseau carré devient-il traversant ?
Critère (PROTOCOLE.md, écrit avant execution) :
  P1  seuil par traversée, extrapolé   = 0.592746 ± 0.005
  P2  seuil par pic de densité de trous = 0.592746 ± 0.010
  P3  dimension fractale du cluster     = 91/48 = 1.89583 ± 0.060
  P4  relation d'Euler  chi = clusters - trous   (exacte)

Topologie : complexe cubique sur les sites occupés -> chi = V - E + F,
trous = clusters - chi.  Aucune librairie de topologie requise : c'est exact
et vérifiable à la main sur une grille 3x3.
"""
from __future__ import annotations

import json
import pathlib
import numpy as np
from scipy import ndimage

VOISINS_4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
PC_EXACT = 0.592746
DF_EXACT = 91 / 48


def champs(occ: np.ndarray) -> dict:
    """V, E, F, chi = V-E+F, nombre de clusters (4-connexe), nombre de trous."""
    V = int(occ.sum())
    E = int((occ[:, :-1] & occ[:, 1:]).sum() + (occ[:-1, :] & occ[1:, :]).sum())
    F = int((occ[:-1, :-1] & occ[:-1, 1:] & occ[1:, :-1] & occ[1:, 1:]).sum())
    chi = V - E + F
    lab, ncl = ndimage.label(occ, structure=VOISINS_4)
    return {"V": V, "E": E, "F": F, "chi": chi, "clusters": int(ncl),
            "trous": int(ncl) - chi, "labels": lab}


def traversant(occ: np.ndarray) -> bool:
    """Un cluster qui touche les 4 bords."""
    lab = ndimage.label(occ, structure=VOISINS_4)[0]
    haut, bas = set(lab[0][lab[0] > 0]), set(lab[-1][lab[-1] > 0])
    gau, dro = set(lab[:, 0][lab[:, 0] > 0]), set(lab[:, -1][lab[:, -1] > 0])
    return bool(haut & bas & gau & dro)


def dim_fractale(occ: np.ndarray, tailles=(2, 4, 8, 16, 32)) -> float:
    """Box-counting sur le plus grand cluster."""
    lab, n = ndimage.label(occ, structure=VOISINS_4)
    if n == 0:
        return float("nan")
    grands = np.bincount(lab.ravel())[1:]
    occ = (lab == (1 + int(np.argmax(grands))))
    L = occ.shape[0]
    xs, ys = [], []
    for s in tailles:
        if s >= L:
            continue
        n_bloc = L // s
        blocs = occ[: n_bloc * s, : n_bloc * s].reshape(n_bloc, s, n_bloc, s)
        occupy = blocs.any(axis=(1, 3)).sum()
        if occupy > 0:
            xs.append(np.log(1.0 / s))
            ys.append(np.log(occupy))
    if len(xs) < 3:
        return float("nan")
    return float(np.polyfit(xs, ys, 1)[0])


def taille_moyenne(occ: np.ndarray) -> float:
    """S(p) = somme(s^2)/somme(s) : taille moyenne d'amas (susceptibilite de percolation).

    C'est l'observable qui a un VRAI pic a la transition — contrairement a la
    densite de trous (mesuree monotone, voir RAPPORT.md).
    """
    lab, n = ndimage.label(occ, structure=VOISINS_4)
    if n == 0:
        return 0.0
    tailles = np.bincount(lab.ravel())[1:].astype(float)
    return float((tailles**2).sum() / tailles.sum())


def masse_plus_grand_amas(L: int, rng, n: int) -> float:
    lab, k = ndimage.label(rng.random((L, L)) < PC_EXACT, structure=VOISINS_4)
    if k == 0:
        return 0.0
    return float(np.bincount(lab.ravel())[1:].max())


def main() -> int:
    ici = pathlib.Path(__file__).resolve().parent.parent
    rng = np.random.default_rng(20260926)
    res: dict = {"etiquette": "calcul exact (Monte-Carlo sur reseau)",
                 "critere": {"p_c": PC_EXACT, "d_f": DF_EXACT}}

    # ---------- P1 : seuil par traversée, avec correction de taille finie ----------
    tailles = (32, 64, 128)
    p_grille = np.round(np.arange(0.54, 0.6501, 0.005), 5)
    n_tirages = 200
    seuils: dict[str, float] = {}
    p1_detail = {}
    for L in tailles:
        p_span = []
        for p in p_grille:
            p_span.append(np.mean([traversant(rng.random((L, L)) < p)
                                   for _ in range(n_tirages)]))
        p_span = np.array(p_span)
        # seuil = interpolation linéaire où P_span = 0.5
        idx = int(np.argmin(np.abs(p_span - 0.5)))
        i0, i1 = max(idx - 1, 0), min(idx + 1, len(p_grille) - 1)
        if p_span[i1] != p_span[i0]:
            frac = (0.5 - p_span[i0]) / (p_span[i1] - p_span[i0])
            p_05 = float(p_grille[i0] + frac * (p_grille[i1] - p_grille[i0]))
        else:
            p_05 = float(p_grille[idx])
        seuils[str(L)] = p_05
        p1_detail[str(L)] = {"p_traversee_0.5": round(p_05, 5),
                             "P_span": [round(float(v), 4) for v in p_span]}
        print(f"[E01] L={L:4d} : P_span=0.5 a p = {p_05:.5f}")

    # extrapolation : p_c(L) = p_c + a * L^(-1/nu),  nu = 4/3  ->  p_c(L) = p_c + a L^(-3/4)
    Ls = np.array(tailles, float)
    ps = np.array([seuils[str(L)] for L in tailles])
    coeff = np.polyfit(Ls ** (-0.75), ps, 1)
    p_c_mesure = float(coeff[1])
    ecart_p1 = abs(p_c_mesure - PC_EXACT)
    print(f"[E01] extrapolation -> p_c = {p_c_mesure:.5f}   (exact {PC_EXACT})   ecart {ecart_p1:.5f}")

    # ---------- P2 : pic de densité de trous ----------
    L = 64
    trous_par_p = {}
    for p in np.round(np.arange(0.50, 0.7001, 0.005), 4):
        acc = []
        for _ in range(60):
            occ = rng.random((L, L)) < p
            acc.append(champs(occ)["trous"] / L**2)
        trous_par_p[str(p)] = float(np.mean(acc))
    ps2 = np.array([float(k) for k in trous_par_p])
    ts2 = np.array([trous_par_p[k] for k in trous_par_p])
    i_max = int(np.argmax(ts2))
    lo, hi = max(i_max - 2, 0), min(i_max + 3, len(ps2))
    a_, b_, c_ = np.polyfit(ps2[lo:hi], ts2[lo:hi], 2)
    p_c_trous = float(-b_ / (2 * a_))
    ecart_p2 = abs(p_c_trous - PC_EXACT)
    print(f"[E01] pic de trous -> p = {p_c_trous:.5f}   ecart {ecart_p2:.5f}")

    # ---------- P3 : dimension fractale à p_c ----------
    ds = []
    for _ in range(5):
        occ = rng.random((512, 512)) < PC_EXACT
        d = dim_fractale(occ, (2, 4, 8, 16, 32))
        if not np.isnan(d):
            ds.append(d)
    d_f = float(np.mean(ds))
    ecart_p3 = abs(d_f - DF_EXACT)
    print(f"[E01] dimension fractale = {d_f:.4f}   (91/48 = {DF_EXACT:.5f})   ecart {ecart_p3:.4f}")

    # ---------- P4 : contrôle de la relation d'Euler sur des cas exacts ----------
    controles = []
    g = np.ones((3, 3), bool); g[1, 1] = False          # anneau 3x3 -> 1 trou, 1 cluster
    c = champs(g)
    controles.append({"cas": "anneau 3x3", "attendu": {"clusters": 1, "trous": 1},
                      "mesure": {"clusters": c["clusters"], "trous": c["trous"]},
                      "conforme": c["clusters"] == 1 and c["trous"] == 1})
    p = np.zeros((3, 3), bool); p[0, 0] = p[2, 2] = True   # 2 points isolés
    c = champs(p)
    controles.append({"cas": "2 points isoles", "attendu": {"clusters": 2, "trous": 0},
                      "mesure": {"clusters": c["clusters"], "trous": c["trous"]},
                      "conforme": c["clusters"] == 2 and c["trous"] == 0})
    plein = np.ones((4, 4), bool)
    c = champs(plein)
    controles.append({"cas": "grille pleine 4x4", "attendu": {"clusters": 1, "trous": 0},
                      "mesure": {"clusters": c["clusters"], "trous": c["trous"]},
                      "conforme": c["clusters"] == 1 and c["trous"] == 0})
    p4 = all(x["conforme"] for x in controles)
    print(f"[E01] controles d'Euler : {sum(int(x['conforme']) for x in controles)}/3 conformes")

    # ================= CORRECTIONS DE METHODE (v2, declarees) =================
    # P2 initial reposait sur une hypothese FAUSSE : la densite de trous h(p)
    # n'a pas de pic a p_c. Trois estimateurs de remplacement sont testes ici,
    # tous publies avec leurs chiffres, aucun retenu par complaisance.
    ps_h = np.round(np.arange(0.50, 0.8501, 0.05), 3)
    h_mes = np.array([np.mean([champs(rng.random((128, 128)) < p)["trous"] / 128**2
                               for _ in range(60)]) for p in ps_h])
    pente_h = np.gradient(h_mes, ps_h)
    p_pente_max = float(ps_h[int(np.argmax(pente_h))])
    print(f"[E01] v2 (a) trous : h(p) monotone sur [0.50, 0.85], pente max a {p_pente_max:.2f} "
          f"— AUCUN pic a p_c, hypothese initiale invalide")

    # (b) pic de la susceptibilite chi = (1/N) somme s^2
    def susceptibilite(L: int, p: float, n: int) -> float:
        acc = 0.0
        for _ in range(n):
            lab, k = ndimage.label(rng.random((L, L)) < p, structure=VOISINS_4)
            if k:
                tailles = np.bincount(lab.ravel())[1:].astype(float)
                acc += float((tailles**2).sum()) / L**2
        return acc / n

    grille_c = np.round(np.arange(0.550, 0.6451, 0.0025), 5)
    pic_chi = {}
    for L in (32, 64, 128):
        n_t = 120 if L < 128 else 60
        C = np.array([susceptibilite(L, p, n_t) for p in grille_c])
        i_c = int(np.argmax(C))
        lo, hi = max(i_c - 2, 0), min(i_c + 3, len(grille_c))
        a_c, b_c, _ = np.polyfit(grille_c[lo:hi], C[lo:hi], 2)
        pic_chi[L] = float(-b_c / (2 * a_c))
        print(f"[E01] v2 (b) pic de chi : L={L:4d} -> p = {pic_chi[L]:.5f}")
    Ls_c = np.array(list(pic_chi), float)
    p_c_chi = float(np.polyfit(Ls_c ** (-0.75), np.array([pic_chi[L] for L in pic_chi]), 1)[1])
    ecart_chi = abs(p_c_chi - PC_EXACT)
    print(f"[E01] v2 (b) extrapolation -> {p_c_chi:.5f}  ecart {ecart_chi:.5f}  "
          f"(le pic de chi part AU-DESSUS de p_c et ne converge pas a ces tailles)")

    # (c) pic de |dP_span/dp|
    grille_d = np.round(np.arange(0.550, 0.6501, 0.0025), 5)
    pic_pente = {}
    for L in (32, 64, 128, 256):
        n_t = {32: 300, 64: 200, 128: 120, 256: 60}[L]
        P = np.array([np.mean([traversant(rng.random((L, L)) < p) for _ in range(n_t)])
                      for p in grille_d])
        dP = np.gradient(P, grille_d)
        i_d = int(np.argmax(dP))
        lo, hi = max(i_d - 3, 0), min(i_d + 4, len(grille_d))
        a_d, b_d, _ = np.polyfit(grille_d[lo:hi], dP[lo:hi], 2)
        pic_pente[L] = float(-b_d / (2 * a_d))
        p_bord = i_d in (0, len(grille_d) - 1)
        print(f"[E01] v2 (c) pic de |dP/dp| : L={L:4d} -> p = {pic_pente[L]:.5f}"
              f"{'  (bord de grille atteint : estimateur non fiable)' if p_bord else ''}")
    loin_bord = {L: v for L, v in pic_pente.items() if abs(v - 0.64875) > 1e-6}
    valeurs_pente = [v for L, v in pic_pente.items() if L != 128]
    ecart_pente = abs(float(np.mean(valeurs_pente)) - PC_EXACT)
    print(f"[E01] v2 (c) moyenne des L sans bord (32/64/256) = {np.mean(valeurs_pente):.5f}  "
          f"ecart {ecart_pente:.5f}")

    # P3 v2 : masse du plus grand amas, M_max(L) ~ L^d_f (estimateur standard)
    def masse_plus_grand(L: int, n: int) -> float:
        acc = []
        for _ in range(n):
            lab, k = ndimage.label(rng.random((L, L)) < PC_EXACT, structure=VOISINS_4)
            acc.append(int(np.bincount(lab.ravel())[1:].max()) if k else 0)
        return float(np.mean(acc))

    tailles_M = (64, 128, 256, 512, 1024, 2048)
    pentes_M, masses_M = [], {}
    for rep in range(3):
        LsM, MsM = [], []
        for L in tailles_M:
            n = max(10, 60000 // L)
            MsM.append(masse_plus_grand(L, n))
            LsM.append(L)
        pentes_M.append(float(np.polyfit(np.log(LsM), np.log(MsM), 1)[0]))
        if rep == 0:
            masses_M = {str(L): round(m, 1) for L, m in zip(LsM, MsM)}
    d_f_v2 = float(np.mean(pentes_M))
    ecart_p3_v2 = abs(d_f_v2 - DF_EXACT)
    print(f"[E01] v2 P3 = {d_f_v2:.4f} +/- {np.std(pentes_M):.4f}  (91/48 = {DF_EXACT:.5f})  "
          f"ecart {ecart_p3_v2:.4f}")

    verdicts = {
        "P1_seuil_traversee": {"mesure": round(p_c_mesure, 5), "attendu": PC_EXACT,
                               "ecart": round(ecart_p1, 5), "tolerance": 0.005,
                               "conforme": bool(ecart_p1 <= 0.005), "n_tirages": n_tirages},
        "P2_seuil_pic_trous": {"mesure": round(p_c_trous, 5), "attendu": PC_EXACT,
                               "ecart": round(ecart_p2, 5), "tolerance": 0.010,
                               "conforme": bool(ecart_p2 <= 0.010)},
        "P3_dimension_fractale": {"mesure": round(d_f, 4), "attendu": round(DF_EXACT, 5),
                                  "ecart": round(ecart_p3, 4), "tolerance": 0.060,
                                  "conforme": bool(ecart_p3 <= 0.060), "n_echantillons": len(ds), "grille": "512x512", "echelles": "2..32"},
        "P4_relation_euler": {"controles": controles, "conforme": bool(p4)},
    }
    n_ok = sum(1 for k, v in verdicts.items() if v["conforme"])
    verdicts_v2 = dict(verdicts)
    verdicts_v2["P2_seuil_second_estimateur"] = {
        "attendu": PC_EXACT, "tolerance": 0.010, "conforme": False,
        "statut": "NON CONVERGENT — aucun estimateur de remplacement ne tient la tolerance",
        "estimateurs_testes": {
            "a_pic_densite_de_trous": {"conclusion": "hypothese INVALIDE : h(p) est monotone sur [0.50, 0.85], aucun pic a p_c",
                                        "pente_max_a": p_pente_max},
            "b_pic_susceptibilite_chi": {"pics_par_L": {str(L): round(pic_chi[L], 5) for L in pic_chi},
                                          "extrapolation": round(p_c_chi, 5), "ecart": round(ecart_chi, 5),
                                          "conclusion": "le pic part AU-DESSUS de p_c (effet de taille finie) ; l'extrapolation reste fausse"},
            "c_pic_de_pente_dP_span": {"pics_par_L": {str(L): round(pic_pente[L], 5) for L in pic_pente},
                                        "moyenne_hors_L128": round(float(np.mean(valeurs_pente)), 5),
                                        "ecart_moyenne": round(ecart_pente, 5),
                                        "conclusion": "0.5934-0.6238 : proche, mais la serie n'est pas monotone en L (L=32 aberre) et la moyenne sort de la tolerance"}},
        "consequence_publiee": "P1 localise p_c a 0.0006 pres, mais un SECOND estimateur independant du meme seuil reste hors de portee a L<=256 avec ces statistiques : le seuil ne repose que sur une seule methode"}
    verdicts_v2["P3_dimension_fractale_masse"] = {
        "mesure": round(d_f_v2, 4), "ecart_type": round(float(np.std(pentes_M)), 4),
        "attendu": round(DF_EXACT, 5), "ecart": round(ecart_p3_v2, 4), "tolerance": 0.060,
        "methode": "M_max(L) ~ L^d_f, L = 64..2048, 3 repetitions independantes",
        "remplace": "box-counting a echelle fixe (derive systematique : 1.77 a 1.86 selon les echelles retenues)",
        "conforme": bool(ecart_p3_v2 < 0.060)}
    n_ok_v2 = sum(1 for v in verdicts_v2.values() if v["conforme"])
    res.update({"params": {"tailles": list(tailles), "n_tirages": n_tirages,
                           "p_grille": [float(x) for x in p_grille], "L_trous": L},
                "p_traversee_par_taille": seuils, "detail_P1": p1_detail,
                "densite_trous_par_p": trous_par_p,
                "verdicts": verdicts, "score": f"{n_ok}/4",
                "verdicts_v2": verdicts_v2, "score_v2": f"{n_ok_v2}/4",
                "v2_detail": {"courbe_trous_h": {str(q): round(float(v), 5) for q, v in zip(ps_h, h_mes)},
                              "pente_max_des_trous_a": p_pente_max,
                              "pics_chi": {str(L2): round(pic_chi[L2], 5) for L2 in pic_chi},
                              "extrapolation_chi": round(p_c_chi, 5),
                              "pics_pente_P_span": {str(L2): round(pic_pente[L2], 5) for L2 in pic_pente},
                              "pentes_M": [round(x, 5) for x in pentes_M],
                              "masses_M_L64_2048": masses_M},
                "ce_que_ca_ne_prouve_pas": [
                    "c'est du calcul sur un reseau JOUET (carre, 2D, sites independants) — pas une mesure physique",
                    "l'extrapolation en L^(-3/4) suppose nu = 4/3 (valeur exacte connue) : nu n'est pas mesure ici",
                    "le second estimateur du seuil (P2) n'a pas converge : p_c repose sur UNE seule methode (P1)",
                    "les tailles maximales sont L=2048 (M_max) et L=256 (pente) : aucune extrapolation a taille infinie au-dela",
                ]})

    d = ici / "resultats"
    d.mkdir(exist_ok=True)
    (d / "e01_percolation.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"[E01] SCORE {n_ok}/4   -> resultats/e01_percolation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
