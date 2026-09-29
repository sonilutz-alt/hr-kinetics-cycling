# =============================================================================
#  REVISION EJAP-D-26-00784 - CELDA 06 - RESULTADOS FINALES + FIGURA 2
#  v1.2  (analisis congelado) - p con F/t y gl de nlme (validado con R); Excel sin redondear
#
#  Decisiones aplicadas:
#   - BBDD revision 1 v1.3; solo se excluye la potencia imposible (41 archivos)
#   - variables del paper con su definicion original; rHRI ELIMINADO
#   - tasa de bajada de FC en VALOR ABSOLUTO (mas alto = mas rapido)
#   - modelo principal: valor ~ cuartil * grupo, intercepto aleatorio por ARCHIVO
#   - grupo Top-10 vs 11-30 con la posicion real
#   - sensibilidad: + carrera (M1), sin Milan-San Remo (M3), posicion continua
#     1-30 (P1), sin cuartiles con cobertura < 80 % (S)
#   - comparaciones de cuartiles: Tukey; grupo dentro de cuartil: Holm (4 Q)
#   - eta2p aproximado con los gl del denominador por nivel (logica de nlme)
#   - analisis de sensibilidad de potencia (no a priori)
#   - antropometria / VO2max (Sitko 2022) solo si existe la plantilla rellena
#
#  Salida en ".../revision 1/":
#     Resultados_finales_v1_<fecha>.xlsx
#     Figura2_resultados_v1.png (300 ppp) y .pdf (vectorial)
# =============================================================================

# ----------------------------------------------------------------------------
#  BLOQUE 1. CONFIGURACION
# ----------------------------------------------------------------------------
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent if "__file__" in globals() else Path.cwd()
REVISION_DIR = BASE / "output"
BBDD = None                     # None = la BBDD_revision1_v1.3_*.xlsx mas reciente
ANTROPOMETRIA = BASE / "data" / "anthropometry.xlsx"   # opcional, no se publica
CORTE_TOP = 10
CP_COLUMNA = "CP_formula_antigua_1_5_12_W"
ALPHA = 0.05

# ----------------------------------------------------------------------------
#  BLOQUE 2. Imports y estilo
# ----------------------------------------------------------------------------
import warnings
from datetime import datetime
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower
from tqdm.auto import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
SELLO = datetime.now().strftime("%Y%m%d_%H%M")
Q = ["Q1", "Q2", "Q3", "Q4"]
G = ["11-30", "Top10"]
CELDAS = [f"{q}_{g}" for q in Q for g in G]

# ----------------------------------------------------------------------------
#  BLOQUE 3. Carga y variables
# ----------------------------------------------------------------------------
if BBDD is None:
    cand = sorted(REVISION_DIR.glob("BBDD_revision1_v1.3_*.xlsx"), key=lambda p: p.stat().st_mtime)
    if not cand:
        raise FileNotFoundError(f"No hay BBDD_revision1_v1.3_*.xlsx en {REVISION_DIR}")
    BBDD = cand[-1]
print(f"BBDD usada: {BBDD.name}")

lar0 = pd.read_excel(BBDD, sheet_name="1_BBDD_larga")
arch0 = pd.read_excel(BBDD, sheet_name="3_Por_archivo_CP_calidad")
flujo = pd.read_excel(BBDD, sheet_name="3b_Flujo_muestra")
arch = arch0[arch0["incluido"] == True].copy()
lar = lar0[lar0["incluido"] == True].merge(arch[["archivo", CP_COLUMNA]], on="archivo", how="left")
for d_ in (lar, arch):
    d_["grupo"] = np.where(d_["posicion"] <= CORTE_TOP, "Top10", "11-30")

lar["rate_HR_increase"] = lar["dFC_subida_bpm_s"]
lar["rate_HR_decrease"] = lar["dFC_bajada_bpm_s"].abs()          # valor absoluto
lar["maximal_HR"] = lar["FCmax_tramos_subida_bpm"]
lar["relative_power"] = lar["pot_tramos_subida_W"] / lar[CP_COLUMNA] * 100
lar["celda"] = pd.Categorical(lar["cuartil"] + "_" + lar["grupo"], categories=CELDAS)

VARS = [("rate_HR_increase", "Rate of HR increase (bpm·s⁻¹)"),
        ("rate_HR_decrease", "Rate of HR decrease (bpm·s⁻¹)"),
        ("maximal_HR", "Maximal HR (bpm)"),
        ("relative_power", "Relative power (%CP)")]

# covariables para sensibilidad (codificacion de suma construida a mano)
CARRERAS = sorted(lar["carrera"].dropna().unique())
for j, c in enumerate(CARRERAS[:-1]):
    lar[f"car_{j}"] = np.where(lar.carrera == c, 1.0, np.where(lar.carrera == CARRERAS[-1], -1.0, 0.0))
CAR_COLS = [f"car_{j}" for j in range(len(CARRERAS) - 1)]
for i, q in enumerate(Q):
    lar[f"qd_{i}"] = (lar.cuartil == q).astype(float)

print(f"Archivos: {lar.archivo.nunique()} | corredores: {lar.corredor.nunique()} | filas: {len(lar)}")

# ----------------------------------------------------------------------------
#  BLOQUE 4. Utilidades del modelo
# ----------------------------------------------------------------------------
def ajusta(d, formula, reml=True):
    ultimo = None
    for metodo in ("lbfgs", "bfgs", "powell", "nm"):
        try:
            fit = smf.mixedlm(formula, d, groups=d["archivo"], re_formula="1").fit(
                reml=reml, method=metodo, maxiter=5000)
            if np.all(np.isfinite(fit.bse.values[:len(fit.fe_params)])):
                return fit
        except Exception as e:
            ultimo = e
    raise RuntimeError(f"No converge: {ultimo}")


def vec(pesos, n):
    v = np.zeros(n)
    for k, w in pesos.items():
        v[CELDAS.index(k)] += w
    return v


def matrices(n):
    LQ = {q: vec({f"{q}_{g}": 0.5 for g in G}, n) for q in Q}
    dif = {q: vec({f"{q}_Top10": 1, f"{q}_11-30": -1}, n) for q in Q}
    return {"cuartil": np.vstack([LQ[q] - LQ["Q4"] for q in Q[:3]]),
            "grupo": vec({**{f"{q}_Top10": .25 for q in Q}, **{f"{q}_11-30": -.25 for q in Q}}, n),
            "cuartil x grupo": np.vstack([dif[q] - dif["Q4"] for q in Q[:3]]),
            "_LQ": LQ, "_dif": dif}


def combo(fit, L):
    b = fit.fe_params.values
    V = np.asarray(fit.cov_params())[:len(b), :len(b)]
    L = np.atleast_2d(L)
    return L @ b, L @ V @ L.T


def wald(fit, L):
    e, c = combo(fit, L)
    chi2 = float(e @ np.linalg.pinv(c) @ e)
    df = np.linalg.matrix_rank(np.atleast_2d(L))
    return chi2, df, float(stats.chi2.sf(chi2, df))


def gl_den(d):
    return {"intra": max(len(d) - d.archivo.nunique() - 6, 1), "entre": max(d.archivo.nunique() - 2, 1)}


def hedges_g(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (na + nb) - 9))


def analiza(d0, var, formula="y ~ 0 + celda"):
    d = d0[d0[var].notna()].copy()
    d["y"] = d[var].astype(float)
    fit = ajusta(d, formula)
    M = matrices(len(fit.fe_params))
    gl = gl_den(d)
    omni = {}
    for ef, niv in [("cuartil", "intra"), ("grupo", "entre"), ("cuartil x grupo", "intra")]:
        chi2, df, _ = wald(fit, M[ef])
        F = chi2 / df
        p = float(stats.f.sf(F, df, gl[niv]))          # como anova() de nlme
        omni[ef] = {"F": F, "gl1": df, "gl2": gl[niv], "p": p, "eta2p": F * df / (F * df + gl[niv])}
    emm = {}
    for q in Q:
        for g_ in G:
            e, c = combo(fit, vec({f"{q}_{g_}": 1}, len(fit.fe_params)))
            emm[(q, g_)] = (float(e[0]), float(np.sqrt(c[0, 0])))
        e, c = combo(fit, M["_LQ"][q])
        emm[(q, "ambos")] = (float(e[0]), float(np.sqrt(c[0, 0])))
    pares = []
    for i in range(4):
        for j in range(i + 1, 4):
            e, c = combo(fit, M["_LQ"][Q[i]] - M["_LQ"][Q[j]])
            se = float(np.sqrt(c[0, 0])); z = float(e[0]) / se
            pares.append({"contraste": f"{Q[i]} vs {Q[j]}", "dif": float(e[0]), "EE": se,
                          "p_Tukey": min(float(stats.studentized_range.sf(abs(z) * np.sqrt(2), 4, gl["intra"])), 1.0)})
    grupo = []
    for q in Q:
        e, c = combo(fit, M["_dif"][q])
        se = float(np.sqrt(c[0, 0]))
        a = d.loc[(d.cuartil == q) & (d.grupo == "Top10"), var]
        b = d.loc[(d.cuartil == q) & (d.grupo == "11-30"), var]
        grupo.append({"cuartil": q, "dif_Top10_menos_11_30": float(e[0]), "EE": se,
                      # t con gl entre-archivos (conservador: el contraste mezcla varianza entre e intra)
                      "IC95_inf": float(e[0]) - stats.t.ppf(0.975, gl["entre"]) * se,
                      "IC95_sup": float(e[0]) + stats.t.ppf(0.975, gl["entre"]) * se,
                      "p": float(2 * stats.t.sf(abs(float(e[0]) / se), gl["entre"])),
                      "g_Hedges": hedges_g(a, b), "n_Top10": int(a.notna().sum()), "n_11_30": int(b.notna().sum())})
    padj = multipletests([r["p"] for r in grupo], method="holm")[1]
    for r, pa in zip(grupo, padj):
        r["p_Holm"] = pa
    return {"fit": fit, "d": d, "omni": omni, "emm": emm, "pares": pares, "grupo": grupo}


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


# ----------------------------------------------------------------------------
#  BLOQUE 5. Modelo principal
# ----------------------------------------------------------------------------
R = {}
for var, _ in tqdm(VARS, desc="Modelo principal"):
    R[var] = analiza(lar, var)

# --- Tabla 1: por cuartil (ambos grupos) + omnibus ---
t1 = []
for var, etiq in VARS:
    r = R[var]
    for q in Q:
        v = r["d"].loc[r["d"].cuartil == q, var]
        m, se = r["emm"][(q, "ambos")]
        t1.append({"Variable": etiq, "Cuartil": q, "n": int(v.notna().sum()),
                   "Media": v.mean(), "DE": v.std(ddof=1),
                   "EMM": m, "IC95_inf": m - 1.96 * se, "IC95_sup": m + 1.96 * se})
t1 = pd.DataFrame(t1)
t1_omni = []
for var, etiq in VARS:
    o = R[var]["omni"]
    sig = [p["contraste"] for p in R[var]["pares"] if p["p_Tukey"] < ALPHA]
    fila = {"Variable": etiq}
    for ef in ["cuartil", "grupo", "cuartil x grupo"]:
        fila[f"{ef} F(gl1,gl2)"] = f"F({o[ef]['gl1']},{o[ef]['gl2']}) = {o[ef]['F']:.2f}"
        fila[f"{ef} p"] = fmt_p(o[ef]["p"])
        fila[f"{ef} eta2p"] = round(o[ef]["eta2p"], 3)
    fila["Contrastes de cuartil significativos (Tukey)"] = ", ".join(sig) if sig else "None"
    t1_omni.append(fila)
t1_omni = pd.DataFrame(t1_omni)
pares = pd.DataFrame([{"Variable": e, **p} for v, e in VARS for p in R[v]["pares"]])

# --- Tabla 2 (suplementaria): grupo x cuartil ---
t2 = []
for var, etiq in VARS:
    r = R[var]
    for gq in r["grupo"]:
        q = gq["cuartil"]
        fila = {"Variable": etiq, "Cuartil": q}
        for g_ in ["Top10", "11-30"]:
            v = r["d"].loc[(r["d"].cuartil == q) & (r["d"].grupo == g_), var]
            m, se = r["emm"][(q, g_)]
            fila[f"{g_} media"] = v.mean(); fila[f"{g_} DE"] = v.std(ddof=1); fila[f"{g_} n"] = int(v.notna().sum())
            fila[f"{g_} EMM IC95"] = f"{m:.2f} [{m - 1.96 * se:.2f}, {m + 1.96 * se:.2f}]"
        fila.update({k: gq[k] for k in ["dif_Top10_menos_11_30", "IC95_inf", "IC95_sup", "g_Hedges", "p", "p_Holm"]})
        t2.append(fila)
t2 = pd.DataFrame(t2)

# ----------------------------------------------------------------------------
#  BLOQUE 6. Sensibilidad (R1.1, R1.2, cobertura)
# ----------------------------------------------------------------------------
sens = []
ESCEN = {"Principal": (lar, "y ~ 0 + celda"),
         "M1 + carrera como covariable": (lar, "y ~ 0 + celda + " + " + ".join(CAR_COLS)),
         "M3 sin Milan-San Remo": (lar[lar.carrera != "MSR"], "y ~ 0 + celda"),
         "Sin cuartiles con cobertura < 80 %": (lar[~lar.cobertura_baja_cuartil.astype(bool)], "y ~ 0 + celda")}
for (esc, (d0, form)), (var, etiq) in tqdm([(e, v) for e in ESCEN.items() for v in VARS], desc="Sensibilidad"):
    r = R[var] if esc == "Principal" else analiza(d0, var, form)
    g1 = [g for g in r["grupo"] if g["cuartil"] == "Q1"][0]
    sens.append({"Escenario": esc, "Variable": etiq, "n_archivos": r["d"].archivo.nunique(),
                 "p cuartil": fmt_p(r["omni"]["cuartil"]["p"]),
                 "p grupo": fmt_p(r["omni"]["grupo"]["p"]),
                 "p cuartil x grupo": fmt_p(r["omni"]["cuartil x grupo"]["p"]),
                 "Q1 dif Top10-(11-30)": round(g1["dif_Top10_menos_11_30"], 3),
                 "Q1 p": fmt_p(g1["p"])})
sens = pd.DataFrame(sens)

# posicion continua 1-30 (R1.2)
pos = []
for var, etiq in tqdm(VARS, desc="Posicion continua"):
    d = lar[lar[var].notna()].copy()
    d["y"] = d[var].astype(float)
    for i in range(4):
        d[f"qx_{i}"] = d[f"qd_{i}"] * (d["posicion"] - d["posicion"].mean())
    fit = ajusta(d, "y ~ 0 + " + " + ".join([f"qd_{i}" for i in range(4)] + [f"qx_{i}" for i in range(4)]))
    nm = list(fit.fe_params.index); n = len(nm)
    ix = [nm.index(f"qx_{i}") for i in range(4)]
    L = np.zeros((3, n))
    for k in range(3):
        L[k, ix[k]] = 1; L[k, ix[3]] = -1
    chi2_i, df_i, _ = wald(fit, L)
    p_int = float(stats.f.sf(chi2_i / df_i, df_i, max(len(d) - d.archivo.nunique() - 6, 1)))
    fila = {"Variable": etiq, "p cuartil x posicion": fmt_p(p_int)}
    for i, q in enumerate(Q):
        b, se = float(fit.fe_params.iloc[ix[i]]), float(fit.bse.iloc[ix[i]])
        fila[f"{q} pendiente por puesto"] = round(b, 4)
        fila[f"{q} p"] = fmt_p(float(2 * stats.t.sf(abs(b / se), d.archivo.nunique() - 2)))
    pos.append(fila)
pos = pd.DataFrame(pos)

# Milan-San Remo descriptiva
msr = (lar[lar.carrera == "MSR"].melt(id_vars=["grupo", "cuartil"], value_vars=[v for v, _ in VARS],
                                      var_name="variable", value_name="valor")
       .dropna().groupby(["variable", "grupo", "cuartil"])["valor"].agg(["size", "mean", "std"]).reset_index())

# ----------------------------------------------------------------------------
#  BLOQUE 7. Muestra, potencia por grupo (R3.3) y potencia estadistica (R3.4)
# ----------------------------------------------------------------------------
pot_vars = [(CP_COLUMNA, "CP (W)"), ("MMP_1min_W", "Best 1-min power (W)"),
            ("MMP_5min_W", "Best 5-min power (W)"), ("MMP_20min_W", "Best 20-min power (W)")]
grp = []
for v, etiq in pot_vars:
    a = arch.loc[arch.grupo == "Top10", v].dropna(); b = arch.loc[arch.grupo == "11-30", v].dropna()
    grp.append({"Variable": etiq, "Top-10": f"{a.mean():.0f} ± {a.std():.0f}", "11th-30th": f"{b.mean():.0f} ± {b.std():.0f}",
                "p (Welch)": fmt_p(stats.ttest_ind(a, b, equal_var=False).pvalue), "g": round(hedges_g(a, b), 2)})

# antropometria y VO2max (Sitko et al. 2022) si la plantilla esta rellena
antro_ok = False
if ANTROPOMETRIA.exists():
    an = pd.read_excel(ANTROPOMETRIA)
    cfn = [c for c in an.columns if c.startswith("fecha_nacimiento")][0]
    an["nac"] = pd.to_datetime(an[cfn], dayfirst=True, errors="coerce")
    an["altura_cm"] = pd.to_numeric(an["altura_cm"], errors="coerce")
    an["peso_kg"] = pd.to_numeric(an["peso_kg"], errors="coerce")
    if an[["nac", "altura_cm", "peso_kg"]].notna().all().all():
        antro_ok = True
        x = arch.merge(an[["corredor", "nac", "altura_cm", "peso_kg"]], on="corredor", how="left")
        # edad a 1 de abril del anyo de la carrera (los monumentos van de marzo a octubre)
        x["edad"] = (pd.to_datetime(x["anyo"].astype(int).astype(str) + "-04-01") - x["nac"]).dt.days / 365.25
        x["VO2max_est"] = 16.6 + 8.87 * (x["MMP_5min_W"] / x["peso_kg"])
        for v, etiq in [("edad", "Age (years)"), ("altura_cm", "Height (cm)"), ("peso_kg", "Body mass (kg)"),
                        ("VO2max_est", "Estimated VO2max (ml·kg⁻¹·min⁻¹)")]:
            a = x.loc[x.grupo == "Top10", v].dropna(); b = x.loc[x.grupo == "11-30", v].dropna()
            grp.append({"Variable": etiq, "Top-10": f"{a.mean():.1f} ± {a.std():.1f}",
                        "11th-30th": f"{b.mean():.1f} ± {b.std():.1f}",
                        "p (Welch)": fmt_p(stats.ttest_ind(a, b, equal_var=False).pvalue), "g": round(hedges_g(a, b), 2)})
        todos = x[["edad", "altura_cm", "peso_kg", "VO2max_est"]].agg(["mean", "std"]).T.round(1)
if not antro_ok:
    print("\nAVISO: antropometria no disponible o incompleta -> sin edad/talla/peso/VO2max.")
grp = pd.DataFrame(grp)

n1 = int((arch.grupo == "Top10").sum()); n2 = int((arch.grupo == "11-30").sum())
pw = TTestIndPower()
d_min = pw.solve_power(nobs1=n1, ratio=n2 / n1, alpha=0.05, power=0.8)
g_q1 = {e: [g for g in R[v]["grupo"] if g["cuartil"] == "Q1"][0]["g_Hedges"] for v, e in VARS[:2]}
potencia = pd.DataFrame([{"n Top-10": n1, "n 11-30": n2, "d minimo detectable (80 %, alfa 0.05)": round(d_min, 2),
                          **{f"potencia para g Q1 {k}": round(pw.power(abs(g), nobs1=n1, ratio=n2 / n1, alpha=0.05), 2)
                             for k, g in g_q1.items()}}])
muestra = pd.DataFrame([{"archivos": len(arch), "corredores": arch.corredor.nunique(),
                         "archivos Top-10": n1, "corredores Top-10": arch.loc[arch.grupo == "Top10", "corredor"].nunique(),
                         "archivos 11-30": n2, "corredores 11-30": arch.loc[arch.grupo == "11-30", "corredor"].nunique(),
                         "corredores en ambos grupos": len(set(arch.loc[arch.grupo == "Top10", "corredor"]) &
                                                          set(arch.loc[arch.grupo == "11-30", "corredor"])),
                         "temporadas": f"{int(arch.anyo.min())}-{int(arch.anyo.max())}",
                         "duracion media (h)": round(arch.n_segundos_rejilla.mean() / 3600, 2)}])
carrera_grupo = pd.crosstab(arch.carrera, arch.grupo, margins=True)

# ----------------------------------------------------------------------------
#  BLOQUE 8. Exportacion
# ----------------------------------------------------------------------------
ruta = REVISION_DIR / f"Resultados_finales_v1.2_{SELLO}.xlsx"
with pd.ExcelWriter(ruta, engine="openpyxl") as w:
    t1.to_excel(w, sheet_name="T1_cuartiles", index=False)
    t1_omni.to_excel(w, sheet_name="T1_omnibus", index=False)
    pares.to_excel(w, sheet_name="T1_pares_Tukey", index=False)
    t2.to_excel(w, sheet_name="S_T2_grupo_cuartil", index=False)
    sens.to_excel(w, sheet_name="S_sensibilidad", index=False)
    pos.to_excel(w, sheet_name="S_posicion_continua", index=False)
    msr.to_excel(w, sheet_name="S_MSR_descriptivos", index=False)
    grp.to_excel(w, sheet_name="Muestra_por_grupo", index=False)
    muestra.to_excel(w, sheet_name="Muestra_resumen", index=False)
    carrera_grupo.to_excel(w, sheet_name="Carrera_x_grupo")
    flujo.to_excel(w, sheet_name="Flujo_muestra", index=False)
    potencia.to_excel(w, sheet_name="Potencia_sensibilidad", index=False)
    if antro_ok:
        todos.to_excel(w, sheet_name="Antropometria_total")
print("Escrito:", ruta)

# ----------------------------------------------------------------------------
#  BLOQUE 9. FIGURA 2 - EMM +/- IC95 por grupo y cuartil
# ----------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "savefig.dpi": 300})
TXT, GRID = "#0b0b0b", "#d9d8d4"
ESTILO = {"Top10": dict(color="#2a78d6", marker="o", ls="-", label="Top-10"),
          "11-30": dict(color="#eb6834", marker="s", ls=(0, (4, 2)), label="11th–30th")}
DESPL = {"Top10": -0.10, "11-30": 0.10}

fig, axes = plt.subplots(2, 2, figsize=(174 / 25.4, 128 / 25.4))
for ax, (var, etiq), letra in zip(axes.flat, VARS, "abcd"):
    r = R[var]
    for g_ in ["Top10", "11-30"]:
        x = np.arange(4) + DESPL[g_]
        m = np.array([r["emm"][(q, g_)][0] for q in Q]); se = np.array([r["emm"][(q, g_)][1] for q in Q])
        st = ESTILO[g_]
        ax.errorbar(x, m, yerr=1.96 * se, color=st["color"], marker=st["marker"], ls=st["ls"],
                    lw=1.2, ms=4.5, capsize=2.5, elinewidth=0.8, mec="white", mew=0.6, label=st["label"])
    top = max(r["emm"][(q, g_)][0] + 1.96 * r["emm"][(q, g_)][1] for q in Q for g_ in G)
    bot = min(r["emm"][(q, g_)][0] - 1.96 * r["emm"][(q, g_)][1] for q in Q for g_ in G)
    rng = top - bot
    for i, gq in enumerate(r["grupo"]):
        if gq["p_Holm"] < ALPHA:
            ax.text(i, top + 0.04 * rng, "*", ha="center", va="bottom", fontsize=11, color=TXT)
    ax.set_ylim(bot - 0.08 * rng, top + 0.18 * rng)
    ax.set_xticks(range(4)); ax.set_xticklabels(Q); ax.set_xlim(-0.4, 3.4)
    ax.set_ylabel(etiq)
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    o = r["omni"]
    pt = lambda x: "p < 0.001" if x < 0.001 else f"p = {x:.3f}"
    ax.set_title(f"Q: {pt(o['cuartil']['p'])}  |  G: {pt(o['grupo']['p'])}  |  "
                 f"Q×G: {pt(o['cuartil x grupo']['p'])}", fontsize=6.5, color="#52514e", pad=4, loc="right")
    ax.text(-0.22, 1.1, letra, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
for ax in axes[1]:
    ax.set_xlabel("Race quartile")
handles = [Line2D([0], [0], color=ESTILO[g_]["color"], marker=ESTILO[g_]["marker"], ls=ESTILO[g_]["ls"],
                  lw=1.2, ms=4.5, mec="white", label=ESTILO[g_]["label"]) for g_ in ["Top10", "11-30"]]
handles.append(Line2D([0], [0], color="none", marker="$*$", ms=5, mfc=TXT, mec=TXT, mew=0.3,
                      label="Top-10 vs 11th–30th within quartile, p < 0.05 (Holm)"))
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.14, hspace=0.45, wspace=0.36)
for ext in ("png", "pdf"):
    fig.savefig(REVISION_DIR / f"Figura2_resultados_v1.1.{ext}", dpi=300)
print("Escrito:", REVISION_DIR / "Figura2_resultados_v1.1.png")

# ----------------------------------------------------------------------------
#  BLOQUE 10. Resumen en pantalla
# ----------------------------------------------------------------------------
pd.set_option("display.width", 220)
print("\n--- Tabla 1 (omnibus) ---"); print(t1_omni.to_string(index=False))
print("\n--- Grupo por cuartil ---")
print(t2[["Variable", "Cuartil", "dif_Top10_menos_11_30", "g_Hedges", "p", "p_Holm"]].round(3).to_string(index=False))
print("\n--- Sensibilidad ---"); print(sens.to_string(index=False))
print("\n--- Posicion continua ---"); print(pos.to_string(index=False))
print("\n--- Muestra ---"); print(muestra.to_string(index=False)); print(grp.to_string(index=False))
print(potencia.to_string(index=False))
