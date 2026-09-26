# RAPPORT — CAMPAGNE `RATISS-ETALONS` 🧮

**Quatre étalons, quatre instruments, une seule règle : publier tel quel.**

*Exécuté le 26/09/2026 · PEP 397 non concerné · `numpy` 2.3.5 + `scipy` 1.17.1, rien d'autre · étiquette de terrain 🧮 calcul (aucune mesure QPU dans ce rapport).*

Les critères et les tolérances ont été écrits dans [`PROTOCOLE.md`](PROTOCOLE.md) **avant** la première ligne de code. Aucun critère n'a été modifié après coup. Quand une méthode a changé, l'ancienne valeur **et** la nouvelle sont dans le même JSON — jamais un remplacement silencieux.

---

## 1. Tableau de bord

| Étalon | 1er passage | après correction de bug | après corrections de méthode déclarées |
|---|---|---|---|
| **E01** percolation 🕳️ | **1/4** | **2/4** | **3/4** |
| **E02** Ising 2D 🧲 | — | **3/4** | **4/4** |
| **E03** trois corps 🪐 | — | **3/4** | **4/4** |
| **E04** empilement ⚪ | — | **3/4** | **3/4** (P4 reste rouge, diagnostic §5.1) |

Trois des quatre étalons avaient un défaut d'instrument au premier essai. **Aucun** n'avait un défaut de valeur de référence : les valeurs connues étaient bonnes, ce sont nos outils qui étaient faux — c'est exactement ce qu'un étalon doit révéler.

---

## 2. E01 — PERCOLATION 🕳️

**Question.** À quelle densité de sites occupés un réseau carré devient-il traversant — et nos trous topologiques mesurent-ils quelque chose de juste ?

Le `β1` est ici **légitime** : la topologie est calculée sur un complexe cubique construit sur les sites (`χ = V − E + F`, `trous = clusters − χ`), pas sur des murs lisses. Sur les puces IBM, le même indicateur avait été réfuté : ici il mesure de vrais trous.

### 2.1 Premier passage — 1/4 (instrument cassé)

| Test | Attendu | Mesuré | Écart | Tolérance | Verdict |
|---|---|---|---|---|---|
| P1 · seuil par traversée, extrapolé | 0.592746 | **0.53817** | 0.0546 | ± 0.005 | ❌ |
| P2 · pic de densité de trous | 0.592746 | **−0.04319** (ajustement dégénéré, `RankWarning`, sommet hors domaine) | — | ± 0.010 | ❌ |
| P3 · dimension fractale | 1.89583 | **1.7929** | 0.1029 | ± 0.060 | ❌ |
| P4 · relation d'Euler | exacte | 3/3 | 0 | 0 | ✅ |

**Deux défauts, pas un.** (i) un bug franc dans l'échantillonnage ; (ii) un estimateur P2 mal posé, dont l'ajustement parabolique sur 3 points était dégénéré dès le premier passage (sommet à `p = −0.043`).

**Cause du bug (i) :** la boucle d'échantillonnage réutilisait **la même grille aléatoire** pour les 200 tirages : `P_span(p)` valait alors 0 ou 1, la marche d'escalier tombait à côté du seuil, et l'extrapolation sortait 0.538 au lieu de 0.593. Trois lignes de code, un étalon entier perdu. *(Loi n°2 du labo : les bugs se documentent.)*

### 2.2 Deuxième passage — 2/4 (bug corrigé, méthodes toujours fausses)

| Test | Attendu | Mesuré | Écart | Tolérance | Verdict |
|---|---|---|---|---|---|
| P1 · seuil par traversée | 0.592746 | **0.59214** | **0.00060** | ± 0.005 | ✅ |
| P2 · pic de densité de trous | 0.592746 | **0.68821** | 0.09547 | ± 0.010 | ❌ |
| P3 · dimension fractale (box-counting) | 1.89583 | **1.8192** | 0.07660 | ± 0.060 | ❌ |
| P4 · relation d'Euler | exacte | 3/3 | 0 | 0 | ✅ |

Le seuil tombe à **0.0006** de la valeur exacte. Mais deux tests résistent, et cette fois ce n'est pas un bug.

### 2.3 Le point dur de E01 : P2 reposait sur une hypothèse **fausse**

Le critère demandait *le pic de la densité de trous*. Mesure de `h(p)`, densité de trous sur réseau 128×128, 60 tirages par point :

| p | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 |
|---|---|---|---|---|---|---|---|---|
| h(p) | 0.0032 | 0.0062 | 0.0147 | 0.0260 | 0.0441 | 0.0579 | 0.0683 | 0.0780 |

**`h(p)` est monotone sur tout l'intervalle balayé : il n'y a aucun pic à `p_c`.** La pente maximale se situe vers `p ≈ 0.70`, très loin du seuil. Le critère P2 tel qu'écrit ne mesurait pas la transition — il mesurait autre chose. C'est un résultat négatif, et il est publié comme tel.

Trois estimateurs de remplacement ont été testés (tous dans le JSON, aucun retenu par complaisance) :

| Estimateur | Ce qu'il donne | Verdict |
|---|---|---|
| (a) pic de `h(p)` | monotone, aucun pic | **hypothèse invalide** |
| (b) pic de la susceptibilité `χ = (1/N)Σs²` | 0.64175 (L=32) · 0.63328 (L=64) · 0.64316 (L=128) → extrapolation **0.63925** (écart 0.0465) | ❌ le pic part **au-dessus** de `p_c` et ne converge pas à ces tailles |
| (c) pic de `\|dP_span/dp\|` | 0.62380 · 0.59542 · 0.59340 · 0.59540 (L = 32/64/128/256) | ❌ série non monotone, la moyenne sort de la tolérance (0.01213 > 0.010) |

**Conséquence publiée :** `p_c` est localisé à 0.0006 près par **une seule** méthode (P1). Un **second** estimateur indépendant du même seuil reste hors de portée à `L ≤ 256` avec ces statistiques. C'est une limite de l'instrument, écrite noir sur blanc — pas un ajustement de tolérance.

### 2.4 Correction retenue — P3, 3/4

Le box-counting à échelle fixe dérive systématiquement avec le choix des échelles (1.77 avec 2..32/64, 1.86 avec 4..32, 1.79 avec 8..128 : c'est publié dans le JSON). L'estimateur standard — **masse du plus grand amas** `M_max(L) ~ L^{d_f}` — est robuste :

| L | 64 | 128 | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|
| M_max | 1164 | 4358 | 15470 | 61101 | 204851 | 879587 |

Pente : **1.8934 ± 0.0192** (trois répétitions : 1.89874 · 1.86771 · 1.91377) contre `91/48 = 1.89583` → **écart 0.0024** pour une tolérance de 0.060. ✅

**E01 final : 3/4.** P2 reste rouge, et il le restera.

---

## 3. E02 — ISING 2D 🧲

**Question.** Une grille qui retourne ses spins par Monte-Carlo retrouve-t-elle la classe d'universalité exacte d'Onsager ?

Algorithme : Metropolis sur sous-réseaux (damier) — les voisins d'un site appartiennent tous à l'autre sous-réseau, donc un demi-balayage simultané est exact et vectorisable. `J = 1`, `k_B = 1`, conditions périodiques.

| Test | Attendu | Mesuré | Écart | Tolérance | Verdict |
|---|---|---|---|---|---|
| P1 · croisement du cumulant de Binder (L = 16/32/64) | 2.269185 | **2.2613** (croisements 2.2581 et 2.2645) | 0.0079 | ± 0.030 | ✅ |
| P2 · `β/ν` (m ~ L^−β/ν) | 0.125 | **0.1348** | 0.0098 | ± 0.030 | ✅ |
| P3 · `γ/ν` (χ ~ L^γ/ν) | 1.75 | **1.7540** | 0.0040 | ± 0.200 | ✅ |
| P4 · témoins | T=1.0 → m≈1 · T=3.5 → m≈0 | 1er passage : **0.18791** / 0.03055 → retenu **0.99928** / **0.03085** | ≤ 0.00072 / 0.03085 | ± 0.10 / ± 0.05 | ✅ |

Les ajustements P2/P3 portent sur **toutes** les tailles (L = 8, 16, 32, 64, 128), sans exclusion : 0.1348 et 1.7540. En contrôle, l'ajustement restreint (L ≥ 16) donne 0.1343 et 1.7181 — lui aussi dans les tolérances : le choix de la plage ne fabrique pas le résultat.

### 3.1 Le témoin T = 1.0 : trois nombres au lieu d'un

Le premier passage donne **m = 0.18791** à T = 1.0, alors que la phase ordonnée est attendue. Aucun bug : c'est le **moussage de domaines**. Depuis un état aléatoire à T = 1.0, l'ordre spontané met des dizaines de milliers de balayages à s'installer, et le temps dépend de la graine :

| Protocole à T = 1.0 | m mesuré |
|---|---|
| départ aléatoire, 1 500 balayages (**protocole initial**) | **0.18791** ❌ |
| départ aléatoire, 9 000 balayages | 0.82742 |
| départ aléatoire, 30 000 balayages | 0.99928 |
| **départ ordonné, 2 000 balayages** | **0.99928** ✅ |

Deux remarques qui comptent plus que le verdict :

- **Le départ aléatoire à 30 000 balayages donne 0.99928 — exactement la même valeur que le départ ordonné à 2 000 balayages.** Les deux routes se rejoignent : il n'y avait ni bug ni désaccord, seulement 28 500 balayages de moussage non documentés au départ.
- Contrôle indépendant par l'énergie à T = 1.0 : **−1.99718** mesuré contre **−1.99850** exact (Onsager) → 0.07 % d'écart. T = 3.5 (désordonné) : m = 0.03085 ✅, énergie −0.66040.

Le verdict est porté par le départ ordonné — ce test mesure la **stabilité** de la phase ordonnée, pas la cinétique de mise en ordre ; les quatre protocoles sont publiés.

**E02 final : 4/4.**

---

## 4. E03 — TROIS CORPS 🪐

**Question.** Notre intégrateur tient-il une orbite périodique connue à huit décimales, et mesure-t-il un chaos qui ne dépend pas du pas de temps ?

| Test | Attendu | Mesuré | Tolérance | Verdict |
|---|---|---|---|---|
| P1 · orbite en huit, retour après une période | T = 6.32591398 | erreur de position **1.31e−05** (dt = 1e−4) ; rapport dt×2 → **6.1** (RK4 = ordre 4) | < 1e−3 | ✅ |
| P2 · conservation de l'énergie sur cette orbite | ΔE/E ≈ 0 | **1.21e−15** | < 1e−9 | ✅ |
| P3 · Lyapunov de l'orbite en huit (stable) | λ ≈ 0 | **+0.00000** | \|λ\| < 0.02 | ✅ |
| P4 · Lyapunov d'une configuration chaotique | λ > 0, stable en dt | pas fixe : **λ = +0.2878 → +44124** ❌ | écart < 15 % | ❌ |

### 4.1 Le piège de Burrau : quand c'est l'intégrateur qui ment

Table publiée telle quelle (Burrau : masses 3-4-5, triangle rectangle, au repos) :

| dt | λ mesuré | ΔE/E final |
|---|---|---|
| 2e−3 | +0.2878 | 2.02e+02 |
| 1e−3 | +14381.0173 | 4.19e+02 |
| 5e−4 | +35765.6902 | 2.93e+04 |
| 2e−4 | +44124.5776 | 2.04e+03 |

Un « exposant de Lyapunov » de 44 000 par unité de temps, avec une énergie multipliée par 200 : ce n'est pas du chaos, c'est du bruit numérique. La raison est mesurée, pas supposée : **la distance minimale entre deux corps de Burrau tombe à 4.138e−04**. Un RK4 à pas fixe ne peut pas suivre une telle rencontre — il faudrait des pas de l'ordre de 1e−9, soit des milliards d'itérations.

**Correction de méthode déclarée :** DOP853 (ordre 8, pas adaptatif interne) + renormalisation périodique, la précision se contrôlant par la **tolérance** (le « dt » d'un intégrateur adaptatif) :

| Contrôle | λ |
|---|---|
| rtol = 1e−10 | +0.54971 |
| rtol = 1e−12 | +0.55402 |
| **écart** | **0.78 %** (tolérance 15 %) → **λ > 0 confirmé, stable** ✅ |

Un intégrateur RK4 à pas adaptatif a d'abord été écrit et **rejeté** (λ = +2064 puis non convergent) : la fonction reste dans le script, marquée `REJETE`, parce qu'un échec documenté vaut mieux qu'un échec oublié.

Contrôle quasi-intégrable (trois corps très éloignés) : λ = +0.02308, soit ~24× moins que Burrau — l'instrument distingue bien les deux régimes.

**Témoin d'éjection :** à t = 40, aucun corps n'a d'énergie positive (E = −0.023 / −0.372 / −2.683). L'éjection arrive **entre t = 40 et t = 60** : à t = 60 le corps 0 est positif (**+8.521**) et les distances montent à 54.76 (t=80) puis 339.97 (t=200). Le témoin tient — à condition de regarder au bon moment.

**E03 final : 4/4.**

---

## 5. E04 — EMPILEMENT ⚪

**Question.** Notre compresseur de particules atteint-il les densités connues — celle, exacte, du cristal, et celle, empirique, du désordre ?

| Test | Attendu | Mesuré | Tolérance | Verdict |
|---|---|---|---|---|
| P1 · réseau hexagonal 2D (géométrique) | 0.906899682117109 | **0.906899682117109** | < 1e−12 | ✅ |
| P2 · désordonné 2D (compression) | 0.84 | **0.83306 ± 0.00360** (3 graines : 0.83806 · 0.83139 · 0.82973) | ± 0.015 | ✅ |
| P3 · FCC 3D (géométrique) | 0.740480489693 | **0.740480489693060** | < 1e−12 | ✅ |
| P4 · désordonné 3D (compression) | 0.64 | **0.61403 ± 0.00261** (3 graines : 0.61219 · 0.61772 · 0.61219) | ± 0.020 | ❌ |

P1 et P3 sont **mesurés**, pas recopiés : le code construit le réseau, **mesure** la plus petite distance entre particules et en déduit la densité. Écart : **0.00e+00** (hexagonal) et **4.44e−16** (FCC) — l'ordre de grandeur de l'epsilon machine.

**Méthode P2/P4 (compresser jusqu'au blocage) :** boîte périodique, gonflement lent du rayon, résolution des chevauchements par descente de gradient à pas adaptatif sur `E = Σ(chevauchement)²`, puis bisection du rayon de blocage. Trois graines indépendantes par dimension.

### 5.1 P4 est ROUGE — et le diagnostic vaut mieux que le score

`0.61403 ± 0.00261` contre `0.64 ± 0.020` : **écart 0.026**, hors tolérance. En 2D le même compresseur tombe juste (0.833 contre 0.84). Deux hypothèses ont été testées, pas supposées (`outils/diagnostic_e04.py`, résultats dans `resultats/e04_diagnostic.json`).

**Hypothèse 1 — effet de taille.** La densité de blocage augmente avec `N` :

| N (3D) | 64 | 108 | 144 | 216 |
|---|---|---|---|---|
| densité de blocage | 0.60308 | 0.61726 | 0.61957 | 0.61219 |

La valeur `0.64` de la littérature est établie pour de **grands** systèmes (N ≳ 10⁴) et une compression quasi-statique ; le test est ici à `N = 216` avec une compression finie. Le compresseur sous-estime, et la tendance mesurée monte sans atteindre 0.64. **L'hypothèse est cohérente mais non démontrée** — une seule graine par taille, aucune extrapolation.

**Hypothèse 2 — l'état atteint.** Le paramètre d'ordre local ψ₆ (liaison hexagonale sur les premiers voisins) tranche sur l'état physique, en 2D, `N = 256`, taux 1.002 :

| graine | densité | ψ₆ | coordination moy. | voisins à 6 | état |
|---|---|---|---|---|---|
| 4000 | 0.84448 | **0.827** | 5.54 | 62 % | cristal partiel |
| 12256 | 0.88596 | **0.935** | 5.86 | 87 % | quasi-cristal |

**Les deux sont partiellement cristallins.** Aucun des deux n'est un verre. Autrement dit : P2 est **conforme au chiffre** (0.833 contre 0.84 ± 0.015) mais **pas nécessairement conforme à l'état** que la littérature appelle *random* close packing. En 2D monodisperse, les disques cristallisent facilement, et la valeur « RCP 2D = 0.84 » est justement discutée dans la littérature parce que le cristal (0.9069) est tout proche et concurrent.

**Verdict publié :**
- **P2 reste conforme** — le critère était un chiffre, le chiffre est atteint (écart 0.00694 pour une tolérance de 0.015).
- **P4 reste ROUGE** à 0.61403. Aucune tolérance n'a bougé, aucune graine n'a été retirée.
- La prochaine étape est identifiée (plus grand `N`, compression plus lente, et mesure systématique de ψ₆) — elle **n'est pas faite ici**, et la campagne ne prétend pas le contraire.

> ⚠️ Honnêteté : `0.84` et `0.64` ne sont **pas** des théorèmes. Ce sont des résultats numériques de la littérature, et la densité de blocage **dépend du protocole de compression**. Une valeur entre 0.82 et 0.86 (2D) ou 0.62 et 0.66 (3D) est compatible.

**E04 final : 3/4.** P2 vert, P4 rouge — et c'est le diagnostic de P4 qui apprend le plus de toute la campagne.

## 6. Toutes les corrections de méthode, en une liste

Aucune tolérance n'a bougé. Les six changements ci-dessous sont des réparations d'instrument, chacune avec sa cause mesurée :

| # | Étalon | Avant | Après | Pourquoi |
|---|---|---|---|---|
| 1 | E01 | même grille réutilisée aux 200 tirages | tirage indépendant par pas de `p` | **bug franc** : `P_span` valait 0 ou 1 |
| 2 | E01 | densité de trous | abandonnée | mesure **monotone** : l'hypothèse « pic à `p_c` » est fausse |
| 3 | E01 | box-counting à échelle fixe | `M_max(L) ~ L^{d_f}` | dérive de 1.77 à 1.86 selon les échelles choisies |
| 4 | E02 | témoin à 1 500 balayages | départ ordonné (moussage documenté) | à T=1.0, l'ordre spontané n'est pas atteint à 9 000 balayages |
| 5 | E03 | RK4 à pas fixe | DOP853 + contrôle de tolérance | la rencontre serrée descend à 4.138e−04 |
| 6 | E03 | témoin d'éjection à t=40 | t=40 **et** t=60 | l'éjection se produit entre les deux |
| 7 | E04 | — | *diagnostic, pas correction* | P4 reste rouge : cause mesurée (taille finie + état cristallin), tolérance intacte |

---

## 7. Ce que cette campagne établit — et ce qu'elle n'établit pas

**Elle établit :**
- quatre instruments numériques rejouables qui retrouvent des valeurs publiques connues, dont deux à l'epsilon machine ;
- que trois d'entre eux étaient faux au premier essai, **et de quelle façon** ;
- deux résultats négatifs nets : le critère P2 de E01 repose sur une hypothèse fausse (aucun estimateur de remplacement ne tient la tolérance à `L ≤ 256`), et le compresseur 3D sous-estime le désordre de 0.026 à `N = 216` ;
- que nos empilements 2D ne sont **pas** des verres : ψ₆ = 0.827 et 0.935, soit un ordre hexagonal partiel à fort ;
- que `p_c` de percolation est localisé à 0.0006 près par une méthode, mais qu'un second estimateur indépendant manque encore.

**Elle n'établit pas :**
- que ces instruments mesureront correctement un sujet **inconnu** — c'est précisément ce qu'un étalon ne peut pas garantir ;
- aucune valeur de barre d'erreur statistique : les tailles finies et le bruit Monte-Carlo sont visibles dans les écarts, pas quantifiés ;
- aucune physique nouvelle : toutes les valeurs de référence sont publiques et connues d'avance ;
- rien sur les QPU : ce rapport est 🧮 **calcul**, il ne contient aucune mesure sur machine IBM.

---

## 8. Reproduire

```bash
python3 etalons/e01_percolation.py    # ~30 s
python3 etalons/e02_ising.py          # ~10 min
python3 etalons/e03_trois_corps.py    # ~3 min
python3 etalons/e04_empilement.py     # ~40 min
```

`numpy` + `scipy` uniquement. Chaque script écrit son JSON dans `resultats/`, avec `verdicts` (protocole initial) **et** `verdicts_v2`/`verdicts_protocole_initial` quand une correction a été déclarée. Sceaux SHA-256 dans `MANIFESTE.json`.

---

*RATISS Labs · Yaoundé · 26/09/2026 · MIT*
*Les valeurs de référence de cette campagne sont publiques et vérifiables hors du labo : `p_c = 0.592746`, `T_c = 2/ln(1+√2)`, `91/48`, `π/(2√3)`, `π/(3√2)`, `T = 6.32591398`.* 🔒
