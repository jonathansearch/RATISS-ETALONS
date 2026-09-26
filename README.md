# RATISS-ETALONS 🧮

**Quatre étalons scientifiques pour vérifier nos instruments avant de leur confier un vrai sujet.**

Un étalon, c'est un problème dont la réponse est déjà connue par la science établie. On écrit le critère **avant** de coder, on lance, et on publie — y compris quand ça rate. Si l'instrument retrouve les valeurs connues, il a le droit de servir ensuite. Sinon, ce sont les instruments qu'on répare.

Campagne `RATISS-ETALONS` — RATISS Labs (Yaoundé) — sous licence MIT.
Étiquette de terrain : **🧮 calcul**, pas une mesure sur QPU. On ne mélange jamais les deux.

---

## Les quatre étalons

| # | Étalon | Ce qu'il teste | Valeur de référence |
|---|---|---|---|
| **E01** | Percolation 2D 🕳️ | seuil critique, dimension fractale, topologie des trous | `p_c = 0.592746` · `d_f = 91/48` |
| **E02** | Ising 2D 🧲 | Monte-Carlo, classe d'universalité | `T_c = 2/ln(1+√2) = 2.269185` · `β/ν = 1/8` · `γ/ν = 7/4` |
| **E03** | Trois corps 🪐 | intégration numérique, chaos | orbite en huit, période `T = 6.32591398` |
| **E04** | Empilement ⚪ | géométrie exacte + compression | `π/(2√3) = 0.906900` (2D) · `π/(3√2) = 0.740480` (3D) · RCP `0.84` / `0.64` |

---

## Résultats

| Étalon | Protocole initial | Après corrections de méthode (déclarées, voir RAPPORT.md) |
|---|---|---|
| E01 percolation 🕳️ | **2/4** | **3/4** — *P2 reste rouge* |
| E02 Ising 2D 🧲 | **3/4** | **4/4** |
| E03 trois corps 🪐 | **3/4** | **4/4** |
| E04 empilement ⚪ | **3/4** | **3/4** — *P4 reste rouge* |

**11 tests sur 16** dans le protocole tel qu'écrit, **14 sur 16** après corrections déclarées. Les deux tests qui restent rouges (E01-P2 et E04-P4) ne sont pas cachés : ils ont chacun leur section dans le rapport, avec leur cause mesurée.

Ce que le premier passage a révélé : un bug franc (E01), un témoin non thermalisé (E02), un intégrateur inadapté à une rencontre serrée (E03), et un compresseur qui sous-estime la densité désordonnée 3D (E04).

Détail complet, échecs compris : **[RAPPORT.md](RAPPORT.md)**.

---

## Rejouer

```bash
python3 etalons/e01_percolation.py     # ~30 s
python3 etalons/e02_ising.py           # ~10 min
python3 etalons/e03_trois_corps.py     # ~3 min
python3 etalons/e04_empilement.py      # ~30 min

python3 outils/diagnostic_e04.py       # diagnostic declare de l'ecart P4 (~15 min)
python3 outils/manifeste.py            # sceau SHA-256 -> MANIFESTE.json
python3 outils/manifeste.py --verifier # verification du sceau
```

Dépendances : `numpy` + `scipy` uniquement. Chaque script écrit son JSON dans `resultats/`.
Aucune valeur n'est ajustée après coup ; les critères sont figés dans `PROTOCOLE.md` (écrit le 26/09/2026 avant exécution).

---

## Les règles appliquées ici

1. **Aucune valeur ajustée après coup.** Les tolérances sont figées avant exécution.
2. **Témoins obligatoires.** Chaque étalon contient au moins un cas limite connu.
3. **Ce qui échoue est publié.** Un test rouge reste dans le rapport, avec son écart.
4. **Corrections déclarées, jamais silencieuses.** Quand une méthode change, l'ancienne valeur et la nouvelle cohabitent dans le même JSON (`verdicts` vs `verdicts_v2`).
5. **🧮 calcul, jamais 🛰️ QPU.** Cette campagne est purement numérique.
6. **Un échec expliqué vaut mieux qu'un succès tu.** Les deux tests rouges sont diagnostiqués, pas contournés.

---

*RATISS Labs · Yaoundé · 2026 · MIT*
