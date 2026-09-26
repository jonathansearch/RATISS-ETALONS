#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Ordre local d'un empilement gele : parametre de liaison psi6 (2D).

Sert a trancher une question posee par E04 : deux configurations produites par
LE MEME protocole de compression donnent 0.83 et 0.886. Sont-elles dans le meme
etat ? psi6 sur les plus proches voisins repond :
  psi6 ~ 0.3-0.5 -> desordre (verre)   |   psi6 > 0.8 -> ordre hexagonal (cristal)

Usage: python3 outils/ordre_local.py            (les trois cas du rapport)
       python3 outils/ordre_local.py 256 4000   (un cas)
"""
from __future__ import annotations

import pathlib
import sys
import time

import numpy as np

RACINE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RACINE / "etalons"))
from e04_empilement import comprimer  # noqa: E402


def psi6(pos: np.ndarray, diametre: float, boite: float, tolerance: float = 1.25) -> dict:
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= boite * np.round(d / boite)
    r = np.sqrt((d**2).sum(-1))
    np.fill_diagonal(r, np.inf)
    voisins = r < tolerance * diametre
    n_voisins = voisins.sum(1)
    angles = np.arctan2(d[:, :, 1], d[:, :, 0])
    z = np.zeros(n, dtype=complex)
    for i in range(n):
        j = np.where(voisins[i])[0]
        if len(j):
            z[i] = np.exp(6j * angles[i, j]).mean()
    coordonnes = n_voisins > 0
    return {"psi6_moyen_individuel": float(np.abs(z[coordonnes]).mean()) if coordonnes.any() else 0.0,
            "coordination_moyenne": float(n_voisins.mean()),
            "fraction_coordination_6": float((n_voisins == 6).mean())}


def main() -> int:
    cas = [(256, 4000), (256, 12256), (144, 12144)]
    if len(sys.argv) == 3:
        cas = [(int(sys.argv[1]), int(sys.argv[2]))]
    for n, graine in cas:
        t0 = time.time()
        r = comprimer(n, 2, graine=graine, n_iter_relax=4000, taux=1.002)
        o = psi6(r["positions"], r["diametre"], r["boite"])
        etat = "CRISTAL" if o["psi6_moyen_individuel"] > 0.8 else "desordre (verre)"
        print(f"N={n} graine={graine} : densite={r['densite']:.5f}  "
              f"psi6={o['psi6_moyen_individuel']:.3f}  coord={o['coordination_moyenne']:.2f}  "
              f"frac6={o['fraction_coordination_6']:.2f}  -> {etat}  ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
