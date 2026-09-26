#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sceau SHA-256 de la campagne : empreinte de chaque fichier du depot.

Usage: python3 outils/manifeste.py [--verifier]
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys

RACINE = pathlib.Path(__file__).resolve().parent.parent
IGNORES = {".git", "__pycache__", "MANIFESTE.json"}


def fichiers() -> list[pathlib.Path]:
    return sorted(p for p in RACINE.rglob("*")
                  if p.is_file() and not any(x in p.parts for x in IGNORES))


def sha256(p: pathlib.Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    cible = RACINE / "MANIFESTE.json"
    if "--verifier" in sys.argv:
        if not cible.exists():
            print("MANIFESTE.json absent")
            return 1
        attendu = json.loads(cible.read_text())["fichiers"]
        ok = manquants = faux = 0
        for nom, empreinte in attendu.items():
            p = RACINE / nom
            if not p.exists():
                print(f"  MANQUANT {nom}")
                manquants += 1
            elif sha256(p) != empreinte:
                print(f"  MODIFIE  {nom}")
                faux += 1
            else:
                ok += 1
        print(f"{ok}/{len(attendu)} conformes · {manquants} manquants · {faux} modifies")
        return 0 if (manquants == 0 and faux == 0) else 1

    manifeste = {str(p.relative_to(RACINE)): sha256(p) for p in fichiers()}
    cible.write_text(json.dumps({"algorithme": "sha256", "fichiers": manifeste},
                                indent=1, ensure_ascii=False, sort_keys=True))
    print(f"{len(manifeste)} fichiers sceelles -> MANIFESTE.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
