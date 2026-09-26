#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""E04 — diagnostic de l'ecart sur la densite desordonnee 3D.

E04 mesure 0.61403 pour la densite desordonnee 3D, soit 0.026 sous la valeur de
litterature (0.64 +- 0.020) : P4 est ROUGE. Ce script cherche la cause au lieu
de la contourner, avec deux mesures independantes :

  1. dependance en taille  : la densite de blocage augmente-t-elle avec N ?
     (0.64 est une valeur pour de GRANDS systemes ; N=216 est petit)
  2. ordre local psi6      : deux configurations du meme protocole donnent 0.83
     et 0.886 en 2D — verre ou cristal ? psi6 tranche.

Ecrit resultats/e04_diagnostic.json. Aucun critere n'est modifie par ce script.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np

RACINE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RACINE / "etalons"))
sys.path.insert(0, str(RACINE / "outils"))
from e04_empilement import comprimer  # noqa: E402
from ordre_local import psi6          # noqa: E402


def main() -> int:
    res: dict = {"etiquette": "diagnostic declare (E04) 🧮",
                 "objet": "expliquer l'ecart de P4 (densite desordonnee 3D) sans toucher aux tolerances"}
    t0 = time.time()

    # ---------- 1. dependance en taille ----------
    taille: dict = {}
    for dim, tailles in ((3, (64, 108, 144, 216)), (2, (64, 144))):
        for n in tailles:
            r = comprimer(n, dim, graine=12000 + n, n_iter_relax=4000, taux=1.002)
            taille[f"{dim}D_N{n}"] = round(r["densite"], 5)
            print(f"[diag] {dim}D N={n:4d} -> {r['densite']:.5f}  ({time.time()-t0:.0f}s)", flush=True)
    res["dependance_en_taille"] = taille

    # ---------- 2. ordre local psi6 (2D, N=256, deux graines) ----------
    ordre = {}
    for graine, nom in ((4000, "graine_verre_0.83"), (12256, "graine_0.886")):
        r = comprimer(256, 2, graine=graine, n_iter_relax=4000, taux=1.002)
        o = psi6(r["positions"], r["diametre"], r["boite"])
        o["densite"] = round(r["densite"], 5)
        o["etat"] = "cristal" if o["psi6_moyen_individuel"] > 0.8 else "desordre (verre)"
        ordre[nom] = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in o.items()}
        print(f"[diag] 2D N=256 {nom} : densite={r['densite']:.5f} "
              f"psi6={o['psi6_moyen_individuel']:.3f} -> {o['etat']}  ({time.time()-t0:.0f}s)", flush=True)
    res["ordre_local_psi6"] = ordre

    res["lecture"] = [
        "la densite de blocage augmente avec N en 3D (0.6031 a N=64 -> 0.6122 a N=216) : une partie de l'ecart a 0.64 est un effet de taille finie",
        "le meme protocole en 2D donne soit un verre (~0.83, psi6 bas) soit un cristal partiel (~0.886, psi6 haut) selon la graine",
        "consequence : 'la' densite de blocage n'est pas une propriete unique du protocole — elle depend de la trajectoire (verre ou cristal), ce qui est un resultat de la litterature sur les milieux granulaires",
        "aucune tolerance n'a ete modifiee : P4 reste ROUGE a 0.61403 contre 0.64 +- 0.020",
    ]
    res["ce_que_ca_ne_prouve_pas"] = [
        "psi6 sur les plus proches voisins ne prouve pas un cristal parfait : il mesure un ordre local hexagonal, avec une tolerance de voisinage arbitraire (1.25 diametres)",
        "l'etude de taille porte sur 1 graine par taille : la dispersion entre graines peut depasser l'effet de taille",
        "aucune extrapolation N -> infini n'est faite : on ne pretend pas atteindre 0.64 par extrapolation",
    ]
    res["duree_s"] = round(time.time() - t0, 1)
    (RACINE / "resultats").mkdir(exist_ok=True)
    (RACINE / "resultats" / "e04_diagnostic.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"[diag] ecrit -> resultats/e04_diagnostic.json ({res['duree_s']:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
