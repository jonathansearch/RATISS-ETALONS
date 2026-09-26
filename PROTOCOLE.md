# 🎯 PROTOCOLE — QUATRE ÉTALONS

**RATISS Labs · 26 septembre 2026**
**Statut : critères écrits AVANT toute execution. Les valeurs attendues sont les valeurs publiées par la science, pas les nôtres.**

> Méthode : on prend **quatre problèmes dont la réponse exacte est connue publiquement**, on les fait tourner
> sur notre moteur, et on publie l'écart. Si l'écart est grand, on le publie aussi.
> **Un étalon qui sort faux est plus instructif qu'un étalon qu'on n'a pas testé.**

---

## Pourquoi ces quatre-là

| # | Étalon | Ce qu'il teste dans nos moteurs | Le chiffre connu |
|---|---|---|---|
| **E01** | Percolation | topologie (clusters, **trous**, dimension fractale) | `p_c = 0.592746` · `d_f = 91/48 ≈ 1.8958` |
| **E02** | Ising 2D | thermodynamique, Monte-Carlo, universalité | `T_c = 2/ln(1+√2) ≈ 2.2692` · `β = 1/8` · `γ = 7/4` |
| **E03** | Trois corps | intégration, chaos, exposant de Lyapunov | période de l'orbite en huit : `T = 6.32591398` |
| **E04** | Empilement | géométrie, densité, compression | `π/(2√3) = 0.9069` (2D) · `π/(3√2) = 0.74048` (3D) · RCP `≈ 0.84` / `≈ 0.64` |

**Deux étalons portent sur des valeurs EXACTES** (dérivables à la main) → tolérance très serrée.
**Deux portent sur des valeurs MESURÉES par simulation** (RCP) → tolérance plus large, et c'est honnête de le dire.

---

## E01 — PERCOLATION 🕳️

**Question :** à quelle densité de sites occupés un réseau carré devient-il traversant — et nos trous topologiques mesurent-ils quelque chose de juste ?

**Critère :**
| Test | Attendu | Tolérance |
|---|---|---|
| P1 · seuil par traversée (`P_span = 0.5`) après extrapolation en taille | `0.592746` | ± 0.005 |
| P2 · seuil par **pic de densité de trous** | `0.592746` | ± 0.010 |
| P3 · **dimension fractale** du cluster traversant (box-counting, à `p_c`) | `91/48 = 1.89583` | ± 0.060 |
| P4 · relation d'Euler `χ = N_clusters − N_trous` vérifiée | exacte | 0 écart |

**Note importante :** ici le `β1` est **légitime** (complexe cubique construit sur les sites). Sur les puces IBM il avait été réfuté parce qu'il mesurait la platitude des murs — ici il mesure de vrais trous.

---

## E02 — ISING 2D 🧲

**Question :** une grille qui retourne ses spins par Monte-Carlo retrouve-t-elle la classe d'universalité exacte d'Onsager ?

**Critère :**
| Test | Attendu | Tolérance |
|---|---|---|
| P1 · température critique (croisement du cumulant de Binder, L = 16/32/64) | `2.269185` | ± 0.030 |
| P2 · exposant `β/ν` (loi `m ~ L^(−β/ν)` à `T_c`) | `0.125` | ± 0.030 |
| P3 · exposant `γ/ν` (loi `χ ~ L^(γ/ν)` à `T_c`) | `1.75` | ± 0.200 |
| P4 · témoins : `T = 1.0` → `m ≈ 1` ; `T = 3.5` → `m ≈ 0` | — | ± 0.10 / ± 0.05 |

---

## E03 — TROIS CORPS 🪐

**Question :** notre intégrateur tient-il une orbite périodique connue à 8 décimales, et mesure-t-il un chaos qui ne dépend pas du pas de temps ?

**Critère :**
| Test | Attendu | Tolérance |
|---|---|---|
| P1 · **orbite en huit** (Chenciner–Montgomery) : retour après une période | `T = 6.32591398` | erreur de position < 1e-3 |
| P2 · conservation de l'énergie sur cette orbite | `ΔE/E ≈ 0` | < 1e-9 |
| P3 · exposant de Lyapunov de l'orbite en huit (orbite **stable**) | `λ ≈ 0` | `|λ| < 0.02` |
| P4 · exposant de Lyapunov d'une configuration **chaotique** | `λ > 0` | et **stable quand `dt` est divisé par 2** (écart < 15 %) |

---

## E04 — EMPILEMENT ⚪

**Question :** notre compresseur de particules atteint-il les densités connues — celle, exacte, du cristal, et celle, empirique, du désordre ?

**Critère :**
| Test | Attendu | Tolérance |
|---|---|---|
| P1 · densité du réseau hexagonal 2D (calcul géométrique) | `π/(2√3) = 0.906900` | < 1e-12 |
| P2 · densité d'empilement **désordonné** 2D (compression) | `0.84` | ± 0.015 |
| P3 · densité de la structure **FCC** 3D (calcul géométrique) | `π/(3√2) = 0.740480` | < 1e-12 |
| P4 · densité d'empilement **désordonné** 3D (compression) | `0.64` | ± 0.020 |

⚠️ Honnêteté : la valeur `0.64` n'est **pas** un théorème — c'est un résultat numérique de la littérature (monodisperse, grand système). Une valeur mesurée entre 0.62 et 0.66 est **compatible**, et le rapport le dira ainsi.

---

## Règles de la campagne

1. **Aucune valeur n'est ajustée après coup.** Les tolérances ci-dessus sont figées.
2. **Témoins obligatoires** : chaque étalon en contient au moins un (limite connue, cas dégénéré).
3. **Ce qui échoue est publié.** Un test rouge n'est pas retiré du rapport : il est signalé, avec son écart.
4. **Étiquette de terrain** : cette campagne est du **calcul** 🧮, pas une mesure sur QPU 🛰️. Jamais de mélange.
5. **Chaque script est rejouable seul**, écrit son JSON, et n'utilise que `numpy` + `scipy`.

---

*Écrit le 26/09/2026 avant exécution. Les valeurs attendues sont publiques et vérifiables hors du labo.* 🔒
