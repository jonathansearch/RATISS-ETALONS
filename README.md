<div align="center">

<img src="assets/banner.jpg" width="100%" alt="RATISS-ETALONS — executable scientific audit · RATISS Labs">

# 🧮 RATISS-ETALONS

**Four scientific standards to verify our instruments before entrusting them with a real subject.**

**RATISS Labs** campaign (Yaoundé) · field tag **🧮 computation** · MIT

`14/16 tests green` · `2 failures published` · `7 declared corrections` · `SHA-256 seal: 17/17`

</div>

---

A standard is a problem whose answer is already known to established science. The criterion is written **before** the first line of code — [`PROTOCOLE.md`](PROTOCOLE.md), frozen on 09/26/2026 — then we run, and we publish **as is**, including when it fails. If the instrument recovers the known values, it earns the right to serve afterwards. Otherwise, it is the instruments that get repaired.

**Result in one sentence:** the four reference values were good — **three instruments out of four were wrong on the first try**, and [the report](RAPPORT.md) says precisely in what way.

---

## 📌 Opening this repository without context? Read in this order

| Order | File | Why |
|---|---|---|
| **1** | [`PROTOCOLE.md`](PROTOCOLE.md) | The criteria and tolerances, **frozen before execution**. Zero value adjusted after the fact. |
| **2** | [`RAPPORT.md`](RAPPORT.md) | The full story: the plain bug, the non-thermalized control, the lying integrator, and the two owned reds. |
| **3** | [`resultats/`](resultats/) | The raw JSONs — `verdicts` (initial protocol) and `verdicts_v2` **cohabit** in the same file. |
| **4** | [`MANIFESTE.json`](MANIFESTE.json) | The 17 SHA-256 fingerprints that seal it all. |

## 🧭 The four standards

| # | Standard | What it tests | Reference value |
|---|---|---|---|
| **E01** | 2D percolation 🕳️ | critical threshold, fractal dimension, hole topology | `p_c = 0.592746` · `d_f = 91/48` |
| **E02** | 2D Ising 🧲 | Monte-Carlo, universality class | `T_c = 2/ln(1+√2) = 2.269185` · `β/ν = 1/8` · `γ/ν = 7/4` |
| **E03** | Three-body 🪐 | numerical integration, chaos | figure-eight orbit, period `T = 6.32591398` |
| **E04** | Packing ⚪ | exact geometry + compression | `π/(2√3) = 0.906900` (2D) · `π/(3√2) = 0.740480` (3D) · RCP `0.84` / `0.64` |

## 📊 Dashboard

| Standard | 1st pass | after bug fix | after declared method corrections |
|---|---|---|---|
| **E01** percolation 🕳️ | **1/4** | **2/4** | **3/4** — *P2 stays red* |
| **E02** 2D Ising 🧲 | — | **3/4** | **4/4** |
| **E03** three-body 🪐 | — | **3/4** | **4/4** |
| **E04** packing ⚪ | — | **3/4** | **3/4** — *P4 stays red* |

**11/16 → 14/16.** The two tests that stay red are not hidden: each has its section in the report, with its **measured** cause — not worked around, no tolerance moved.

## 🎯 The numbers that count

| Quantity | Reference | Measured | Gap |
|---|---|---|---|
| Percolation threshold `p_c` | 0.592746 | 0.59214 | **0.00060** (tol. 0.005) |
| Fractal dimension `d_f` (M_max) | 1.89583 | 1.8934 ± 0.0192 | 0.0024 (tol. 0.060) |
| Ising `T_c` (Binder crossing) | 2.269185 | 2.2613 | 0.0079 (tol. 0.030) |
| Ising `γ/ν` | 1.75 | 1.7540 | 0.0040 (tol. 0.200) |
| Figure-eight orbit — return at `T` | 6.32591398 | position error 1.31e−05 | < 1e−3 |
| E03 energy conservation | ΔE/E ≈ 0 | **1.21e−15** | < 1e−9 |
| Burrau Lyapunov (DOP853) | λ > 0, stable in tolerance | 0.55402 (vs 0.54971 at rtol 1e−10) | 0.78% (tol. 15%) |
| 2D hexagonal | 0.906899682117109 | 0.906899682117109 | **0.00e+00** |
| 3D FCC | 0.740480489693 | 0.740480489693060 | 4.44e−16 — *machine epsilon* |

## 📈 The figures — plotted from the sealed JSONs

R7 applied to the pixel: **no decorative image**. Every curve comes from a file in `resultats/` (sealed by the manifest), and everything regenerates with one command:

```bash
python3 outils/figures.py    # → assets/fig_*.png (numpy + scipy + matplotlib)
```

<div align="center">

**E01 — a plain bug, a wrong hypothesis, published as is**

<img src="assets/fig_e01_percolation.png" width="100%" alt="E01 percolation: spanning transition, monotone h(p), M_max ~ L^d_f">

**E02 — Onsager's universality class recovered**

<img src="assets/fig_e02_ising.png" width="100%" alt="E02 2D Ising: Binder crossing, chi ~ L^(gamma/nu)">

**E03 — the integrator was caught by chaos**

<img src="assets/fig_e03_trois_corps.png" width="100%" alt="E03 three-body: figure-eight orbit, fixed-step vs DOP853 Lyapunov">

**E04 — P4 stays red, diagnosis first**

<img src="assets/fig_e04_empilement.png" width="100%" alt="E04 packing: size dependence N, local order psi6">

</div>

One honest precision, written into the script: on figure E03, the **trajectory** of the figure-eight orbit is re-integrated on the fly with the repository's initial conditions (DOP853, rtol 1e-12 — return error 5.03e-08); the **verdicts**, on their side, stay the sealed JSON's ones. This is the only recomputed value of all the figures.

## 🔴 The two red tests — owned

### E01-P2 · it was the hypothesis that was wrong, not the instrument

The criterion asked for the peak of the hole density at `p_c`. Published measurement: `h(p)` is **monotone** over the whole scanned interval — there is no peak. Three replacement estimators tested (susceptibility peak, peak of `|dP_span/dp|`), **none holds the tolerance** at `L ≤ 256`. Published consequence, black on white: `p_c` is localized to within 0.0006 by **one single** method; the second independent estimator remains out of reach. *(RAPPORT §2.3)*

### E04-P4 · the diagnosis is worth more than the score

`0.61403 ± 0.00261` against `0.64 ± 0.020`: out of tolerance, period. Two hypotheses tested, not assumed (`outils/diagnostic_e04.py`): finite size (blocking density 0.603 → 0.620 from `N = 64` to `216`, rising trend) and the state reached — **ψ₆ = 0.827 and 0.935: our packings are partially crystalline, not glasses**. P2 matches the number (0.833 against 0.84) but not necessarily the *random* state described by the literature. The continuation (larger `N`, slower compression, systematic ψ₆) is identified — it **is not done here**, and the campaign does not pretend otherwise. *(RAPPORT §5.1)*

> ⚠️ And a reminder: `0.84` and `0.64` are **not** theorems — they are numerical results that depend on the compression protocol.

## 🔧 The method corrections — all declared

No tolerance moved. Seven changes, each with its measured cause:

| # | Standard | Before | After | Why |
|---|---|---|---|---|
| 1 | E01 | same grid reused over 200 draws | independent draw | **plain bug**: `P_span` was 0 or 1 |
| 2 | E01 | hole density | abandoned | monotone measurement: invalid hypothesis |
| 3 | E01 | box-counting at fixed scale | `M_max(L) ~ L^{d_f}` | drift 1.77 → 1.86 depending on scales |
| 4 | E02 | control, 1,500 sweeps | ordered start | documented domain foaming: 0.18791 → 0.99928 |
| 5 | E03 | fixed-step RK4 | DOP853 + tolerance control | close encounter at 4.138e−04: λ = +44124, noise, not chaos |
| 6 | E03 | ejection control at t=40 | t=40 **and** t=60 | the ejection happens between the two |
| 7 | E04 | — | *diagnosis, not correction* | P4 stays red, tolerance intact |

The old value and the new one cohabit in the same JSON. An adaptive RK4 integrator was even written then **rejected** for E03 — the function stays in the script, tagged `REJETE` (REJECTED): a documented failure is worth more than a forgotten one. *(Lab law #2: bugs get documented.)*

## ⚖️ What the campaign establishes — and what it does not

**It establishes:**
- four replayable numerical instruments that recover known public values, two of them at machine epsilon;
- that three of them were wrong on the first try, **and in what way**;
- two clean negative results: the P2 criterion of E01 rests on a wrong hypothesis, and the 3D compressor underestimates disorder by 0.026 at `N = 216`;
- that our 2D packings are **not** glasses (ψ₆ = 0.827 and 0.935);
- that `p_c` is localized to within 0.0006 by one method — a second independent estimator is still missing.

**It does not establish:**
- that these instruments will correctly measure an **unknown** subject — that is precisely what a standard cannot guarantee;
- any complete statistical error bar: finite sizes and Monte-Carlo noise are visible in the gaps, not quantified;
- any new physics: all reference values are public and known in advance;
- **anything about QPUs**: this report is 🧮 **computation**, no measurement on an IBM machine. The two are never mixed.

## ▶️ Replay

```bash
git clone https://github.com/jonathansearch/RATISS-ETALONS.git
cd RATISS-ETALONS

python3 etalons/e01_percolation.py     # ~30 s
python3 etalons/e02_ising.py           # ~10 min
python3 etalons/e03_trois_corps.py     # ~3 min
python3 etalons/e04_empilement.py      # ~40 min

python3 outils/diagnostic_e04.py       # declared diagnosis of the P4 gap (~15 min)
python3 outils/manifeste.py --verifier # SHA-256 seal verification
```

Dependencies: `numpy` + `scipy`, nothing else *(campaign executed with numpy 2.3.5 + scipy 1.17.1)*. Each script writes its JSON into `resultats/`.

> **R7**: every claim of this README is replayable with a command above. A stranger must be able to verify without asking permission.

## 🔒 Verify integrity

```bash
python3 outils/manifeste.py --verifier   # 17/17 SHA-256 fingerprints compliant
```

Any modification of a repository file breaks the seal — that is on purpose. After any legitimate edit: `python3 outils/manifeste.py` to re-seal, and the modification goes to the journal.

## 📜 The rules applied here

1. **No value adjusted after the fact.** Tolerances are frozen before execution (`PROTOCOLE.md`, 09/26/2026).
2. **Mandatory controls.** Each standard contains at least one known edge case.
3. **What fails is published.** A red test stays in the report, with its gap.
4. **Declared corrections, never silent.** Old and new value in the same JSON.
5. **🧮 computation, never 🛰️ QPU.** Purely numerical campaign.
6. **An explained failure is worth more than a killed success.** The two reds are diagnosed, not worked around.

## 🧬 RATISS Labs ecosystem

| Repository / link | Role |
|---|---|
| [`RATISS-ARCHIVES`](https://github.com/jonathansearch/RATISS-ARCHIVES) | The lab's memory: evidence, registry of the 86 IBM Quantum tasks, identity |
| [`DISCORD-RATISS`](https://github.com/jonathansearch/DISCORD-RATISS) | The agent that posts the evidence in the lab's server |
| [`RATISS-QVM`](https://github.com/jonathansearch/RATISS-QVM) · [`RATISS-NAVIER`](https://github.com/jonathansearch/RATISS-NAVIER) · [`RATISS-FUSION`](https://github.com/jonathansearch/RATISS-FUSION) | The lab's engines — allowed to serve after passing the standard |
| [Official website](https://jonathansearch.github.io/ratiss-labs-site/) · [ORCID](https://orcid.org/0009-0000-4092-5313) | The executable scientific audit, off GitHub |

---

*RATISS Labs · Yaoundé · campaign of 09/26/2026 · MIT*
*The reference values of this campaign are public and verifiable outside the lab: `p_c = 0.592746`, `T_c = 2/ln(1+√2)`, `91/48`, `π/(2√3)`, `π/(3√2)`, `T = 6.32591398`.* 🔒
