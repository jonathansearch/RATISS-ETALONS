#!/usr/bin/env python3
"""FIGURES — RATISS-ETALONS : toutes les figures du README, régénérées ici.

R7 appliqué au pixel : aucune figure décorative. Chaque courbe est tracée
depuis `resultats/*.json` (données scellées par MANIFESTE.json). La seule
chose recalculée est la trajectoire de l'orbite en huit (E03, panneau
gauche), réintégrée avec les conditions initiales du dépôt — son verdict,
lui, reste celui du JSON scellé.

Usage :
    python3 outils/figures.py            # écrit assets/fig_*.png

Dépendances : numpy, scipy, matplotlib (aucune autre).
Style : identité du labo — fond #060813, turquoise/cyan, texte clair.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

RACINE = pathlib.Path(__file__).resolve().parent.parent
RESULTATS = RACINE / "resultats"
ASSETS = RACINE / "assets"

# --- identité visuelle du labo -------------------------------------------
FOND = "#060813"
PANNEAU = "#0b1122"
GRILLE = "#1e2a45"
TEXTE = "#e2e8f0"
MUTED = "#8ba3c7"
TURQUOISE = "#2dd4bf"
CYAN = "#38bdf8"
MENTHE = "#99f6e4"
BLANC = "#f1f5f9"
VERT = "#4ade80"
ROUGE = "#f87171"


def style_labo() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": FOND,
            "axes.facecolor": PANNEAU,
            "savefig.facecolor": FOND,
            "axes.edgecolor": GRILLE,
            "axes.labelcolor": TEXTE,
            "axes.titlecolor": BLANC,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "text.color": TEXTE,
            "grid.color": GRILLE,
            "font.size": 10.5,
            "axes.titlesize": 11.5,
            "axes.titleweight": "bold",
            "legend.facecolor": PANNEAU,
            "legend.edgecolor": GRILLE,
        }
    )


def nouvelle_fig(nrows: int, ncols: int, titre: str, **kwargs):
    """Figure avec layout contraint + bande de pied de page réservée."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows, ncols, layout="constrained", **kwargs)
    fig.get_layout_engine().set(rect=(0, 0.030, 1, 0.970))
    fig.suptitle(titre, color=BLANC, fontsize=13, fontweight="bold")
    return fig, ax


def pied(fig, source: str) -> None:
    fig.text(
        0.995,
        0.008,
        f"RATISS Labs · RATISS-ETALONS · campagne du 26/09/2026 · calcul — données : {source}",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=MUTED,
    )


def charger(nom: str) -> dict:
    with open(RESULTATS / nom, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------- E01 --
def fig_e01() -> None:
    import matplotlib.pyplot as plt

    d = charger("e01_percolation.json")
    p_c = d["critere"]["p_c"]
    p_grille = np.array(d["params"]["p_grille"])
    couleurs = {32: MENTHE, 64: TURQUOISE, 128: CYAN}

    fig, ax = nouvelle_fig(
        1, 3, figsize=(15.6, 4.9),
        titre="E01 · PERCOLATION 2D — un bug franc, une hypothèse fausse, publiés tels quels",
    )

    # -- A : P_span(p) par taille, transition traversante ----------------
    a = ax[0]
    for L, coul in couleurs.items():
        detail = d["detail_P1"][str(L)]
        a.plot(p_grille, detail["P_span"], color=coul, lw=1.8, label=f"L = {L} (mesuré)")
        p_t = d["p_traversee_par_taille"][str(L)]
        a.plot(p_t, 0.5, "o", ms=7, mfc=coul, mec=BLANC, mew=0.8, zorder=5)
    for i, (L, coul) in enumerate(couleurs.items()):
        p_t = d["p_traversee_par_taille"][str(L)]
        a.text(0.035, 0.965 - i * 0.075, f"L = {L} :  p_t = {p_t:.4f}",
               transform=a.transAxes, ha="left", va="top", fontsize=8.5, color=coul,
               fontweight="bold",
               bbox=dict(facecolor=FOND, edgecolor=coul, alpha=0.9, boxstyle="round,pad=0.25"))
    a.axhline(0.5, color=MUTED, lw=0.8, ls=":")
    a.axvline(p_c, color=BLANC, lw=1.1, ls="--", alpha=0.85)
    a.text(p_c + 0.003, 1.03, f"référence p_c = {p_c}", fontsize=8.5, color=BLANC)
    a.set_xlim(0.53, 0.665)
    a.set_xlabel("p — fraction de sites occupés")
    a.set_ylabel("P_traversée")
    a.set_title("P1 · transition traversante")
    a.legend(loc="lower right", fontsize=8.5)
    a.set_ylim(-0.04, 1.06)
    a.grid(alpha=0.5)

    # -- B : h(p) monotone -> hypothèse P2 invalide (négatif publié) -----
    a = ax[1]
    trous = d["densite_trous_par_p"]
    ps = np.array([float(k) for k in trous])
    hs = np.array([trous[k] for k in trous])
    a.plot(ps, hs, color=TURQUOISE, lw=2)
    a.axvline(p_c, color=BLANC, lw=1.1, ls="--", alpha=0.85)
    a.text(p_c + 0.003, max(hs) * 0.97, "p_c", color=BLANC, fontsize=9)
    a.annotate(
        "aucun pic : h(p) est monotone\n→ hypothèse P2 INVALIDE\n(résultat négatif publié)",
        xy=(p_c, hs[np.argmin(abs(ps - p_c))]),
        xytext=(p_c + 0.015, max(hs) * 0.40),
        color=ROUGE, fontsize=9,
        arrowprops=dict(arrowstyle="->", color=ROUGE, lw=1.2),
    )
    a.set_xlabel("p — fraction de sites occupés")
    a.set_ylabel("densité de trous h(p)")
    a.set_title("P2 · le critère ne mesurait pas la transition")
    a.grid(alpha=0.5)

    # -- C : M_max(L) ~ L^d_f, la correction retenue ---------------------
    a = ax[2]
    masses = d["v2_detail"]["masses_M_L64_2048"]
    tailles = np.array([float(k) for k in masses])
    m = np.array([masses[k] for k in masses])
    a.loglog(tailles, m, "o", color=TURQUOISE, ms=7, mec=BLANC, mew=0.7, label="M_max (mesuré)")
    ref = m[0] * (tailles / tailles[0]) ** (91.0 / 48.0)
    a.loglog(tailles, ref, "--", color=BLANC, lw=1.1, alpha=0.8, label="pente 91/48 = 1.89583")
    pentes = d["v2_detail"]["pentes_M"]
    a.set_title(f"P3 · M_max ~ L^(d_f) — pente {np.mean(pentes):.4f} ± {np.std(pentes, ddof=1):.4f}")
    a.set_xlabel("L (taille du réseau)")
    a.set_ylabel("M_max — masse du plus grand amas")
    a.legend(loc="upper left", fontsize=8.5)
    a.grid(alpha=0.5, which="both")

    pied(fig, "resultats/e01_percolation.json")
    fig.savefig(ASSETS / "fig_e01_percolation.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- E02 --
def fig_e02() -> None:
    import matplotlib.pyplot as plt

    d = charger("e02_ising.json")
    T_c = d["critere"]["T_c"]
    p1 = d["P1_detail"]
    T = np.array(p1["T"])

    fig, ax = nouvelle_fig(
        1, 2, figsize=(12.8, 4.9),
        titre="E02 · ISING 2D — 3/4 → 4/4 après thermalisation déclarée du témoin",
    )

    # -- A : croisement du cumulant de Binder ----------------------------
    a = ax[0]
    couleurs = {16: MENTHE, 32: TURQUOISE, 64: CYAN}
    for L, coul in couleurs.items():
        a.plot(T, p1["u4"][str(L)], color=coul, lw=1.8, label=f"L = {L} (mesuré)")
    a.axvline(T_c, color=BLANC, lw=1.1, ls="--", alpha=0.85)
    a.text(T_c - 0.004, 0.34, f"T_c = 2/ln(1+√2)\n= {T_c:.6f}", color=BLANC,
           fontsize=8.5, ha="right")
    a.set_xlabel("T (J = k_B = 1)")
    a.set_ylabel("cumulant de Binder u₄")
    a.set_title("P1 · croisement de Binder → T_c ≈ 2.2613")
    a.legend(loc="upper right", fontsize=8.5)
    a.grid(alpha=0.5)

    # -- B : chi ~ L^(gamma/nu) ------------------------------------------
    a = ax[1]
    det = d["P2_P3_detail"]["chi_par_L"]
    tailles = np.array([float(k) for k in det])
    chi = np.array([det[k] for k in det])
    pente = np.polyfit(np.log(tailles), np.log(chi), 1)[0]
    a.loglog(tailles, chi, "o", color=TURQUOISE, ms=8, mec=BLANC, mew=0.7,
             label="χ (mesuré, L = 8 → 128)")
    xs = np.linspace(tailles.min() * 0.9, tailles.max() * 1.08, 50)
    a.loglog(xs, np.exp(np.polyval(np.polyfit(np.log(tailles), np.log(chi), 1), np.log(xs))),
             "--", color=BLANC, lw=1.1, label=f"ajustement : γ/ν = {pente:.4f}")
    a.text(0.03, 0.94, "référence Onsager : γ/ν = 7/4 = 1.75",
           transform=a.transAxes, color=VERT, fontsize=8.5)
    a.set_xlabel("L (taille de la grille)")
    a.set_ylabel("χ — susceptibilité magnétique")
    a.set_title("P3 · χ ~ L^(γ/ν) → 1.7540 (écart 0.0040)")
    a.legend(loc="lower right", fontsize=8.5)
    a.grid(alpha=0.5, which="both")

    pied(fig, "resultats/e02_ising.json")
    fig.savefig(ASSETS / "fig_e02_ising.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- E03 --
def fig_e03() -> None:
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    d = charger("e03_trois_corps.json")

    # conditions initiales de l'orbite en huit — identiques au script E03
    r1 = np.array([0.97000436, -0.24308753])
    v3 = np.array([-0.93240737, -0.86473146])
    pos0 = np.array([r1, -r1, np.zeros(2)])
    vit0 = np.array([-v3 / 2, -v3 / 2, v3])
    T8 = d["params"]["T_8"]

    def rhs(t, y):
        p = y[:6].reshape(3, 2)
        dv = np.zeros((3, 2))
        for i in range(3):
            for j in range(3):
                if i != j:
                    diff = p[j] - p[i]
                    dv[i] += diff / np.linalg.norm(diff) ** 3
        return np.concatenate([y[6:], dv.ravel()])

    y0 = np.concatenate([pos0.ravel(), vit0.ravel()])
    sol = solve_ivp(rhs, (0.0, T8), y0, method="DOP853",
                    rtol=1e-12, atol=1e-15, dense_output=True, max_step=0.01)
    err_retour = float(np.linalg.norm(sol.y[:6, -1].reshape(3, 2) - pos0))

    e0 = d["params"]["E0_orbite_8"]
    p2 = d.get("P2_detail") or {}
    dE = p2.get("dE_sur_E")
    if dE is None:
        m = np.ones(3)

        def ener(y):
            p = y[:6].reshape(3, 2)
            v = y[6:].reshape(3, 2)
            ec = 0.5 * np.sum(m[:, None] * v**2)
            ep = 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    ep -= m[i] * m[j] / np.linalg.norm(p[j] - p[i])
            return ec + ep

        dE = abs((ener(sol.y[:, -1]) - e0) / e0)
        src_dE = "réintégration locale"
    else:
        src_dE = "JSON scellé"

    ts = np.linspace(0.0, T8, 900)
    pts = sol.sol(ts)[:6].reshape(3, 2, -1)

    fig, ax = nouvelle_fig(
        1, 2, figsize=(12.8, 5.1),
        titre="E03 · TROIS CORPS — l'intégrateur prenait le chaos en défaut",
    )

    # -- A : l'orbite en huit --------------------------------------------
    a = ax[0]
    for i, (coul, nom) in enumerate(zip([TURQUOISE, CYAN, MENTHE], ["corps 1", "corps 2", "corps 3"])):
        a.plot(pts[i][0], pts[i][1], color=coul, lw=1.6, label=nom)
        a.plot(pos0[i][0], pos0[i][1], "o", ms=7, color=coul, mec=BLANC, mew=0.8, zorder=5)
    a.set_aspect("equal")
    a.set_ylim(-0.80, 0.62)
    a.set_title(f"orbite en huit — une période T = {T8:.8f}")
    a.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
    a.grid(alpha=0.4)
    a.set_xlabel("x")
    a.set_ylabel("y")
    a.text(0.5, -0.74,
           f"réintégrée ici avec les CI du dépôt (DOP853, rtol 1e-12)\n"
           f"erreur de retour {err_retour:.2e} · ΔE/E {dE:.2e} ({src_dE}) · verdicts : JSON scellé",
           transform=a.transData, ha="center", va="bottom", fontsize=8, color=MUTED)

    # -- B : l'intégrateur mentait, pas le chaos -------------------------
    a = ax[1]
    table = d["P4_detail"]["table_pas_fixe"]
    ordre = ["dt=0.002", "dt=0.001", "dt=0.0005", "dt=0.0002"]
    dts = np.array([float(k.split("=")[1]) for k in ordre])
    lam = np.array([table[k]["lambda"] for k in ordre])
    a.loglog(dts, lam, "o-", color=ROUGE, lw=1.6, ms=7, mec=BLANC, mew=0.6,
             label="RK4 à pas fixe (Burrau) : λ explose → bruit numérique")
    b10 = d["P4_detail"]["burrau_rtol1e-10"]["lambda"]
    b12 = d["P4_detail"]["burrau_rtol1e-12"]["lambda"]
    a.axhspan(0.55 * 0.85, 0.55 * 1.15, color=VERT, alpha=0.12)
    a.axhline(0.55, color=VERT, lw=1.0, ls="--", alpha=0.8)
    a.plot([1e-10], [b10], "D", color=VERT, ms=8, mec=BLANC, mew=0.6,
           label="DOP853 rtol 1e-10 → λ = 0.54971")
    a.plot([1e-12], [b12], "D", color=MENTHE, ms=8, mec=BLANC, mew=0.6,
           label="DOP853 rtol 1e-12 → λ = 0.55402")
    quasi = d["P4_detail"]["controle_large_quasi_integrable"]["lambda"]
    a.plot([2e-11], [quasi], "s", color=MUTED, ms=7, mec=BLANC, mew=0.6,
           label=f"contrôle quasi-intégrable → λ = {quasi:.5f}")
    dmin = d["P4_detail"]["distance_minimale"]["distance_minimale"]
    a.annotate(f"distance minimale entre corps :\n{dmin:.3e} — un pas fixe ne peut pas suivre",
               xy=(2e-4, 4.4e4), xytext=(3e-10, 8e3), color=ROUGE, fontsize=8.5,
               arrowprops=dict(arrowstyle="->", color=ROUGE, lw=1.1))
    a.set_xlabel("pas / tolérance de l'intégrateur (échelle log)")
    a.set_ylabel("exposant de Lyapunov λ")
    a.set_title("P4 · λ > 0 confirmé (écart rtol 0.78 %, tol 15 %)")
    a.legend(loc="lower left", fontsize=8)
    a.grid(alpha=0.5, which="both")

    pied(fig, "resultats/e03_trois_corps.json")
    fig.savefig(ASSETS / "fig_e03_trois_corps.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- E04 --
def fig_e04() -> None:
    import matplotlib.pyplot as plt

    diag = charger("e04_diagnostic.json")
    emp = charger("e04_empilement.json")

    fig, ax = nouvelle_fig(
        1, 2, figsize=(12.8, 4.9),
        titre="E04 · EMPILEMENT — P4 reste ROUGE, le diagnostic d'abord",
    )

    # -- A : dépendance en taille + bandes de référence ------------------
    a = ax[0]
    dep = diag["dependance_en_taille"]
    n3 = [64, 108, 144, 216]
    d3 = [dep[f"3D_N{n}"] for n in n3]
    n2 = [64, 144]
    d2 = [dep[f"2D_N{n}"] for n in n2]

    a.axhspan(0.64 - 0.020, 0.64 + 0.020, color=VERT, alpha=0.10)
    a.axhline(0.64, color=VERT, lw=1.0, ls="--", alpha=0.85)
    a.text(60, 0.644, "référence 3D : 0.64 ± 0.020 (tolérance)", color=VERT, fontsize=8.3)
    a.plot(n3, d3, "o-", color=TURQUOISE, lw=1.7, ms=7, mec=BLANC, mew=0.6,
           label="désordonné 3D (diagnostic, 1 graine/taille)")

    a.axhspan(0.84 - 0.015, 0.84 + 0.015, color=CYAN, alpha=0.10)
    a.axhline(0.84, color=CYAN, lw=1.0, ls="--", alpha=0.7)
    a.text(60, 0.8455, "référence 2D : 0.84 ± 0.015", color=CYAN, fontsize=8.3)
    a.plot(n2, d2, "s--", color=CYAN, lw=1.4, ms=6, mec=BLANC, mew=0.6,
           label="désordonné 2D (diagnostic)")

    p4 = emp["P4_detail"]
    a.errorbar([216], [p4["moyenne"]], yerr=[p4["ecart_type"]], fmt="D", color=ROUGE,
               ms=9, mec=BLANC, mew=0.8, capsize=4, zorder=6,
               label=f"P4 campagne : {p4['moyenne']:.5f} ± {p4['ecart_type']:.5f} → ROUGE")
    a.set_xlabel("N — nombre de particules")
    a.set_ylabel("densité de blocage")
    a.set_title("P4 · dépendance en taille N (diagnostic)")
    a.legend(loc="upper right", fontsize=8)
    a.set_ylim(0.585, 0.875)
    a.grid(alpha=0.5)

    # -- B : psi6 — pas des verres ---------------------------------------
    a = ax[1]
    psi = diag["ordre_local_psi6"]
    graines = list(psi.keys())
    vals = [psi[g]["psi6_moyen_individuel"] for g in graines]
    dens = [psi[g]["densite"] for g in graines]
    fracs = [psi[g]["fraction_coordination_6"] for g in graines]
    x = np.arange(len(graines))
    a.bar(x, vals, width=0.45, color=[TURQUOISE, CYAN], edgecolor=BLANC, linewidth=0.7)
    for xi, (v, de, fr) in enumerate(zip(vals, dens, fracs)):
        a.text(xi, v + 0.02, f"ψ₆ = {v:.3f}", ha="center", color=BLANC,
               fontsize=10, fontweight="bold")
        a.text(xi, 0.05, f"ρ = {de:.3f}\n{int(round(fr * 100))} % de liaisons à 6",
               ha="center", color=PANNEAU, fontsize=8.3, fontweight="bold")
    a.axhline(1.0, color=BLANC, lw=1.0, ls="--", alpha=0.7)
    a.text(1.42, 0.975, "cristal parfait", ha="right", va="top", color=BLANC, fontsize=8.3)
    a.set_xticks(x, ["graine 4000\ncristal partiel", "graine 12256\nquasi-cristal"], fontsize=9)
    a.set_ylim(0, 1.12)
    a.set_ylabel("ψ₆ — ordre hexagonal local")
    a.set_title("ψ₆ — l'état atteint n'est pas un verre")

    pied(fig, "resultats/e04_empilement.json + e04_diagnostic.json")
    fig.savefig(ASSETS / "fig_e04_empilement.png", dpi=150)
    plt.close(fig)


def main() -> None:
    style_labo()
    ASSETS.mkdir(exist_ok=True)
    fig_e01()
    fig_e02()
    fig_e03()
    fig_e04()
    print("4 figures écrites dans assets/ :")
    for f in sorted(ASSETS.glob("fig_*.png")):
        print(" -", f.relative_to(RACINE), f"({f.stat().st_size // 1024} Ko)")


if __name__ == "__main__":
    main()
