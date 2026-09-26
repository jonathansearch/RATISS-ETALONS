#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""E04 — EMPILEMENT ⚪ : le cristal (exact) et le desordre (mesure).

Question : notre compresseur de particules atteint-il les densites connues —
celle, exacte, du cristal, et celle, empirique, du desordre ?
Critere (PROTOCOLE.md, ecrit avant execution) :
  P1  densite du reseau hexagonal 2D (geometrique)  = pi/(2*sqrt(3)) = 0.906900   (< 1e-12)
  P2  densite d'empilement desordonne 2D (compression) = 0.84  (± 0.015)
  P3  densite de la structure FCC 3D (geometrique)  = pi/(3*sqrt(2)) = 0.740480   (< 1e-12)
  P4  densite d'empilement desordre 3D (compression) = 0.64  (± 0.020)

Methode : compresseur "gonfle et relaxe" (inflate-and-relax) — boite periodique
fixe, on gonfle le rayon, on resout les chevauchements par relaxation de paires,
et on bissecte le rayon de blocage (jamming). Aucune librairie de dynamique
moleculaire : numpy seul.
"""
from __future__ import annotations

import json
import pathlib
import numpy as np
from scipy.spatial import ConvexHull

TAUX_GONFLEMENT = 1.001
TOL_CHEVAUCHEMENT = 1e-4


def distance_minimale_pbc(pos: np.ndarray, boite: float) -> float:
    """Plus petite distance entre deux particules, minimum-image."""
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= boite * np.round(d / boite)
    r = np.sqrt((d**2).sum(-1))
    np.fill_diagonal(r, np.inf)
    return float(r.min())


def distances_pbc(pos: np.ndarray, boite: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= boite * np.round(d / boite)
    r = np.sqrt((d**2).sum(-1))
    np.fill_diagonal(r, np.inf)
    return d, r


def relaxer(pos: np.ndarray, boite: float, rayon: float, n_iter: int = 4000,
            alpha0: float = 0.5,
            tol: float = TOL_CHEVAUCHEMENT) -> tuple[np.ndarray, float]:
    """Descente de gradient sur l'energie des chevauchements E = somme(ov^2).

    Le pas est adaptatif (x1,3 si E baisse, /2 sinon) : descente monotone.
    Le critere porte sur E (lisse) et non sur max(ov) (non lisse : reduire le
    pire chevauchement peut en creer un autre et figer la descente).
    """
    alpha = alpha0
    d, r = distances_pbc(pos, boite)
    chev = r < 2.0 * rayon
    i, j = np.where(np.triu(chev, 1))
    if len(i) == 0:
        return pos, 0.0
    ov = 2.0 * rayon - r[i, j]
    energie = float((ov**2).sum())
    ov_max = float(ov.max())
    for _ in range(n_iter):
        if ov_max < tol:
            return pos, ov_max
        d, r = distances_pbc(pos, boite)
        chev = r < 2.0 * rayon
        i, j = np.where(np.triu(chev, 1))
        if len(i) == 0:
            return pos, 0.0
        ov = 2.0 * rayon - r[i, j]
        ov_max = float(ov.max())
        if ov_max < tol:
            return pos, ov_max
        dirv = d[i, j] / r[i, j][:, None]
        w = ov[:, None] * dirv
        delta = np.zeros_like(pos)
        np.add.at(delta, i, w)
        np.add.at(delta, j, -w)
        pas = alpha * delta
        norme = float(np.abs(pas).max())
        if norme > 0.25 * rayon:
            pas *= (0.25 * rayon) / norme
        essai = (pos + pas) % boite
        d2, r2 = distances_pbc(essai, boite)
        chev2 = r2 < 2.0 * rayon
        i2, j2 = np.where(np.triu(chev2, 1))
        if len(i2) == 0:
            return essai, 0.0
        ov2 = 2.0 * rayon - r2[i2, j2]
        energie2 = float((ov2**2).sum())
        if energie2 < energie:
            pos, energie = essai, energie2
            alpha *= 1.3
        else:
            alpha *= 0.5
            if alpha < 1e-18:
                break
    d, r = distances_pbc(pos, boite)
    chev = r < 2.0 * rayon
    i, j = np.where(np.triu(chev, 1))
    return pos, (0.0 if len(i) == 0 else float((2.0 * rayon - r[i, j]).max()))


def comprimer(n_particules: int, dimension: int, graine: int, boite: float = 1.0,
              n_iter_relax: int = 4000, densite_depart: float = 0.40,
              taux: float = TAUX_GONFLEMENT):
    """Gonfle le rayon jusqu'au blocage (jamming) et bissecte le rayon critique.

    Protocole declare : config de depart = points uniformes aleatoires (donc
    desordonnee), amenee a densite 0.40 par relaxation, puis compression lente
    de taux constant. La densite de blocage depend du taux (c'est un resultat,
    pas un defaut) : le rapport le publie.
    """
    rng = np.random.default_rng(graine)
    pos = rng.random((n_particules, dimension)) * boite
    volume = boite**dimension
    unitaire = np.pi * (2.0) ** 2 / 4.0 if dimension == 2 else 4.0 / 3.0 * np.pi
    rayon = (densite_depart * volume / (n_particules * unitaire)) ** (1.0 / dimension)
    pos, _ = relaxer(pos, boite, rayon, n_iter=n_iter_relax)

    # phase 1 : gonflement exponentiel jusqu'a l'echec de la relaxation
    r_ok, pos_ok = rayon, pos.copy()
    for _ in range(20000):
        rayon = r_ok * taux
        pos_new, ov = relaxer(pos_ok.copy(), boite, rayon, n_iter_relax)
        if ov < TOL_CHEVAUCHEMENT:
            r_ok, pos_ok = rayon, pos_new
        else:
            break
    # phase 2 : bissection du rayon de blocage
    r_bas, r_haut = r_ok, rayon
    for _ in range(30):
        r_mil = 0.5 * (r_bas + r_haut)
        pos_new, ov = relaxer(pos_ok.copy(), boite, r_mil, n_iter=n_iter_relax)
        if ov < TOL_CHEVAUCHEMENT:
            r_bas, pos_ok = r_mil, pos_new
        else:
            r_haut = r_mil
    densite = n_particules * unitaire * r_bas**dimension / volume
    return {"densite": float(densite), "rayon_jamming": float(r_bas),
            "positions": pos_ok, "boite": boite,
            "n_particules": n_particules, "dimension": dimension,
            "taux_compression": taux, "densite_depart": densite_depart,
            "distance_min": distance_minimale_pbc(pos_ok, boite),
            "diametre": 2.0 * r_bas}


def densite_hexagonale(n: int = 24) -> dict:
    """Reseau hexagonal 2D : la distance min est MESUREE, puis la densite en decoule."""
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3.0) / 2.0])
    sites = np.array([i * a1 + j * a2 for i in range(n) for j in range(n)])
    # boite periodique engendree par n*a1 et n*a2 -> aucune surface libre
    boite = np.array([n * a1, n * a2])
    aire = abs(float(np.linalg.det(boite)))
    d_min = float(np.linalg.norm(a1))
    n_part = n * n
    densite = n_part * np.pi * (d_min / 2.0) ** 2 / aire
    return {"densite": float(densite), "distance_min_mesuree": d_min,
            "aire_boite": aire, "n_particules": n_part,
            "aire_par_particule": aire / n_part}


def densite_fcc(n: int = 4) -> dict:
    """FCC : la distance min est MESUREE sur un supercellulaire periodique."""
    base = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    cote = float(n)
    pos = np.array([(np.array([i, j, k]) + b) / n for i in range(n) for j in range(n)
                    for k in range(n) for b in base])
    # boite periodique de cote 1 (coordonnees normalisees), distance min mesuree
    d_min = distance_minimale_pbc(pos, 1.0)
    volume = 1.0
    n_part = len(pos)
    densite = n_part * (4.0 / 3.0) * np.pi * (d_min / 2.0) ** 3 / volume
    return {"densite": float(densite), "distance_min_mesuree": float(d_min),
            "volume_boite": volume, "n_particules": n_part,
            "volume_par_particule": volume / n_part}


def main() -> int:
    ici = pathlib.Path(__file__).resolve().parent.parent
    res: dict = {"etiquette": "calcul exact (geometrie + compression) 🧮",
                 "critere": {"hex_2d": float(np.pi / (2 * np.sqrt(3.0))),
                             "fcc_3d": float(np.pi / (3 * np.sqrt(2.0))),
                             "rcp_2d": 0.84, "rcp_3d": 0.64}}

    # ---------- P1 : hexagonal 2D ----------
    hexa = densite_hexagonale(24)
    exact_2d = float(np.pi / (2 * np.sqrt(3.0)))
    ecart_p1 = abs(hexa["densite"] - exact_2d)
    print(f"[E04] P1 hexagonal 2D : {hexa['densite']:.15f}  (exact {exact_2d:.15f})  "
          f"ecart {ecart_p1:.2e}", flush=True)

    # ---------- P3 : FCC 3D ----------
    fcc = densite_fcc(3)
    exact_3d = float(np.pi / (3 * np.sqrt(2.0)))
    ecart_p3 = abs(fcc["densite"] - exact_3d)
    print(f"[E04] P3 FCC 3D      : {fcc['densite']:.15f}  (exact {exact_3d:.15f})  "
          f"ecart {ecart_p3:.2e}", flush=True)

    # ---------- P2 : desordre 2D (compression) ----------
    d2 = []
    for g in range(3):
        r = comprimer(256, 2, graine=4000 + g, taux=1.001)
        d2.append(r["densite"])
        print(f"[E04] P2 compression 2D graine {g} : densite = {r['densite']:.5f}  "
              f"(diametre {r['diametre']:.5f})", flush=True)
    d2 = np.array(d2)
    dens2 = float(d2.mean())
    std2 = float(d2.std())
    ecart_p2 = abs(dens2 - 0.84)
    print(f"[E04] P2 moyenne 2D = {dens2:.5f} +/- {std2:.5f}  ecart {ecart_p2:.5f}", flush=True)

    # ---------- P4 : desordre 3D (compression) ----------
    d3 = []
    for g in range(3):
        r = comprimer(216, 3, graine=8000 + g, n_iter_relax=4000, taux=1.002)
        d3.append(r["densite"])
        print(f"[E04] P4 compression 3D graine {g} : densite = {r['densite']:.5f}  "
              f"(diametre {r['diametre']:.5f})", flush=True)
    d3 = np.array(d3)
    dens3 = float(d3.mean())
    std3 = float(d3.std())
    ecart_p4 = abs(dens3 - 0.64)
    print(f"[E04] P4 moyenne 3D = {dens3:.5f} +/- {std3:.5f}  ecart {ecart_p4:.5f}", flush=True)

    verdicts = {
        "P1_hexagonal_2D": {"mesure": hexa["densite"], "attendu": exact_2d,
                            "ecart": ecart_p1, "tolerance": 1e-12,
                            "conforme": bool(ecart_p1 < 1e-12)},
        "P2_desordre_2D": {"mesure": round(dens2, 5), "ecart_type": round(std2, 5),
                           "attendu": 0.84, "ecart": round(ecart_p2, 5), "tolerance": 0.015,
                           "conforme": bool(ecart_p2 < 0.015)},
        "P3_FCC_3D": {"mesure": fcc["densite"], "attendu": exact_3d,
                      "ecart": ecart_p3, "tolerance": 1e-12,
                      "conforme": bool(ecart_p3 < 1e-12)},
        "P4_desordre_3D": {"mesure": round(dens3, 5), "ecart_type": round(std3, 5),
                           "attendu": 0.64, "ecart": round(ecart_p4, 5), "tolerance": 0.020,
                           "conforme": bool(ecart_p4 < 0.020)},
    }
    n_ok = sum(1 for v in verdicts.values() if v["conforme"])
    res.update({
        "params": {"taux_gonflement": TAUX_GONFLEMENT, "tol_chevauchement": TOL_CHEVAUCHEMENT,
                   "n_2d": 256, "n_3d": 216, "graines_2d": 3, "graines_3d": 3,
                   "taux_2d": 1.001, "taux_3d": 1.002, "densite_depart": 0.40,
                   "methode": "inflate-and-relax + bisection du rayon de blocage, boite periodique"},
        "P1_detail": hexa, "P3_detail": fcc,
        "P2_detail": {"densites": [round(float(x), 5) for x in d2],
                      "moyenne": round(dens2, 5), "ecart_type": round(std2, 5)},
        "P4_detail": {"densites": [round(float(x), 5) for x in d3],
                      "moyenne": round(dens3, 5), "ecart_type": round(std3, 5)},
        "verdicts": verdicts, "score": f"{n_ok}/4",
        "ce_que_ca_ne_prouve_pas": [
            "0.84 et 0.64 sont des resultats numeriques de la litterature, pas des theoremes : une valeur proche est compatible, pas demontree",
            "la densite de blocage depend du protocole (vitesse de compression) : un autre protocole peut cristalliser vers 0.9069 en 2D",
            "disques/sphères parfaitement monodisperses, sans frottement, sans gravite : aucun lien avec un materiau reel",
            "le rayon de blocage est bissecte sur une configuration gelee : d'autres configurations gelees existent",
        ]})
    d = ici / "resultats"
    d.mkdir(exist_ok=True)
    (d / "e04_empilement.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"[E04] SCORE {n_ok}/4   -> resultats/e04_empilement.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
