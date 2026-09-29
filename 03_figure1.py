# =============================================================================
#  REVISION EJAP-D-26-00784 - CELDA 05 - FIGURA 1: ESQUEMA DEL METODO
#  v1.2  - leyendas debajo de cada panel (ancho completo); sin solapes
#
#  Figura con datos reales de un archivo que muestra, paso a paso, como se
#  obtienen las variables (R1.11, y aclara R1.7/R1.12):
#     a  potencia de toda la carrera (media movil de 30 s) y cuartiles Q1-Q4
#     b  FC de toda la carrera, mismo eje temporal
#     c  detalle: potencia a 1 s (suavizada 2 s) con los tramos de subida y
#        bajada detectados (umbral +/- 0,3 DE de la derivada de potencia)
#     d  detalle: FC suavizada; tramos conservados (relleno) y descartados por
#        el filtro de signo (rayado)
#  Mismo preprocesado que la celda 02 v1.3 (version corregida).
#  Salida: ".../revision 1/Figura1_metodo_v1.2.png" (300 ppp) y ".pdf" (vectorial)
# =============================================================================

# ----------------------------------------------------------------------------
#  BLOQUE 1. CONFIGURACION
# ----------------------------------------------------------------------------
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent if "__file__" in globals() else Path.cwd()
ARCHIVO = BASE / "data" / "raw" / "example_race.xlsx"   # cualquier archivo de carrera
OUTPUT_DIR = BASE / "output"

VENTANA_S = 60             # duracion del detalle (paneles c y d)
INICIO_VENTANA_S = None    # None = se elige sola (ver BLOQUE 4); o un segundo concreto

POWER_MAX, HR_MIN, HR_MAX = 1800, 30, 220
MAX_INTERP_S = 5
THRESH_SD = 0.3

# ----------------------------------------------------------------------------
#  BLOQUE 2. Imports y estilo
# ----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})
TXT = "#0b0b0b"          # texto y FC
TXT2 = "#52514e"         # potencia y elementos secundarios
GRID = "#d9d8d4"
C_SUB = "#2a78d6"        # tramos de subida de potencia
C_BAJ = "#eb6834"        # tramos de bajada de potencia

# ----------------------------------------------------------------------------
#  BLOQUE 3. Lectura y preprocesado (identico a la celda 02, version corregida)
# ----------------------------------------------------------------------------
def a_segundos(serie):
    if np.issubdtype(serie.dtype, np.datetime64):
        s = pd.to_datetime(serie, errors="coerce")
        return (s - s.min()).dt.total_seconds().values
    num = pd.to_numeric(serie, errors="coerce")
    if num.notna().mean() > 0.95:
        return (num - num.min()).values.astype(float)
    s = pd.to_datetime(serie, errors="coerce", utc=True)
    return (s - s.min()).dt.total_seconds().values


def na_approx(x, max_gap):
    s = pd.Series(np.asarray(x, float))
    isn = s.isna()
    run_len = isn.groupby((isn != isn.shift()).cumsum()).transform("sum")
    f = s.interpolate(limit_area="inside")
    f[isn & (run_len > max_gap)] = np.nan
    return f.values


def rollmean_k2(x):
    x = np.asarray(x, float)
    y = np.empty_like(x)
    y[:-1] = (x[:-1] + x[1:]) / 2
    y[-1] = y[-2]
    return y


def find_sequences(v, thr, direction):
    seqs, start = [], None
    for i, val in enumerate(v):
        if np.isnan(val):
            if start is not None:
                seqs.append((start, i)); start = None
            continue
        cond = val > thr if direction == "increase" else val < -thr
        if cond:
            if start is None:
                start = i
        elif start is not None:
            seqs.append((start, i)); start = None
    if start is not None:
        seqs.append((start, len(v) - 1))
    return seqs


d = pd.read_excel(ARCHIVO)
t_raw = a_segundos(d["timestamp"])
ok = ~np.isnan(t_raw)
df = pd.DataFrame({"s": np.round(t_raw[ok]).astype(int),
                   "p": pd.to_numeric(d["power"], errors="coerce").values[ok],
                   "h": pd.to_numeric(d["heart_rate"], errors="coerce").values[ok]})
df.loc[(df.p > POWER_MAX) | (df.p < 0), "p"] = np.nan
df.loc[(df.h < HR_MIN) | (df.h > HR_MAX), "h"] = np.nan
g = df.groupby("s")[["p", "h"]].mean()
g = g.reindex(np.arange(g.index.min(), g.index.max() + 1))
t = (g.index.values - g.index.min()).astype(float)
p = na_approx(g.p.values, MAX_INTERP_S)
h = na_approx(g.h.values, MAX_INTERP_S)
ps, hs = rollmean_k2(p), rollmean_k2(h)
p_der = np.r_[np.nan, np.diff(ps)]
fc_der = np.r_[np.nan, np.diff(hs)]
thr = THRESH_SD * np.nanstd(p_der, ddof=1)

tramos = []
for tipo in ("increase", "decrease"):
    for s_, e_ in find_sequences(p_der, thr, tipo):
        m = np.nanmean(fc_der[s_:e_ + 1]) if np.any(~np.isnan(fc_der[s_:e_ + 1])) else np.nan
        conservado = (not np.isnan(m)) and ((m >= 0) if tipo == "increase" else (m <= 0))
        tramos.append({"tipo": tipo, "ini": s_, "fin": e_, "fc_media": m, "conservado": conservado})
tramos = pd.DataFrame(tramos)
edges = np.linspace(0, t.max(), 5)

# ----------------------------------------------------------------------------
#  BLOQUE 4. Eleccion de la ventana de detalle
#  Criterio (sin mirar resultados): datos completos en la ventana, dentro de Q2,
#  y la ventana con mas tramos conservados de ambos tipos.
# ----------------------------------------------------------------------------
if INICIO_VENTANA_S is None:
    mejor, puntos = None, -1
    for ini in range(int(edges[1]), int(edges[2]) - VENTANA_S, 10):
        fin = ini + VENTANA_S
        if np.isnan(ps[ini:fin]).any() or np.isnan(hs[ini:fin]).any():
            continue
        w = tramos[(tramos.ini >= ini) & (tramos.fin < fin)]
        n_up = int(((w.tipo == "increase") & w.conservado).sum())
        n_dn = int(((w.tipo == "decrease") & w.conservado).sum())
        sc = min(n_up, n_dn) * 100 + n_up + n_dn
        if sc > puntos:
            mejor, puntos = ini, sc
    INICIO_VENTANA_S = mejor
W0, W1 = INICIO_VENTANA_S, INICIO_VENTANA_S + VENTANA_S
print(f"Ventana de detalle: {W0}-{W1} s ({W0/3600:.2f}-{W1/3600:.2f} h)")

# ----------------------------------------------------------------------------
#  BLOQUE 5. Figura
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(174 / 25.4, 175 / 25.4))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1.25, 1.25], hspace=0.62)
axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
axC, axD = fig.add_subplot(gs[2]), fig.add_subplot(gs[3])

th = t / 3600
p30 = pd.Series(p).rolling(60, min_periods=40, center=True).mean().values
h30 = pd.Series(h).rolling(60, min_periods=40, center=True).mean().values

for ax, y, col, lab in [(axA, p30, TXT2, "Power (W)"), (axB, h30, TXT, "HR (bpm)")]:
    ax.plot(th, y, color=col, lw=0.6)
    ax.set_ylabel(lab)
    ax.set_xlim(0, th.max())
    for e in edges[1:-1] / 3600:
        ax.axvline(e, color=TXT2, lw=0.6, ls=(0, (3, 2)))
    ax.axvspan(W0 / 3600, W1 / 3600, color=GRID, alpha=0.9, lw=0, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
axA.tick_params(labelbottom=False)
axB.set_xlabel("Race time (h)")
for i in range(4):
    axA.text((edges[i] + edges[i + 1]) / 2 / 3600, 1.03, f"Q{i + 1}",
             transform=axA.get_xaxis_transform(), ha="center", va="bottom",
             color=TXT, fontsize=7, clip_on=False)
axA.set_ylim(0, np.nanmax(p30) * 1.05)

# --- detalle ---
sl = slice(W0, W1)
tz = t[sl] - W0
axC.plot(tz, ps[sl], color=TXT2, lw=0.8)
axD.plot(tz, hs[sl], color=TXT, lw=0.8)
w = tramos[(tramos.fin >= W0) & (tramos.ini < W1)]
for _, r in w.iterrows():
    # el tramo propiamente dicho va de ini a fin-1 (fin es la fila que rompe la
    # racha, convencion heredada del R); asi no se solapan tramos contiguos
    x0, x1 = max(r.ini, W0) - W0 - 0.5, min(r.fin - 1, W1 - 1) - W0 + 0.5
    if x1 <= x0:
        continue
    col = C_SUB if r.tipo == "increase" else C_BAJ
    axC.axvspan(x0, x1, color=col, alpha=0.28, lw=0)
    if r.conservado:
        axD.axvspan(x0, x1, color=col, alpha=0.28, lw=0)
    else:
        axD.axvspan(x0, x1, facecolor="none", edgecolor=col, hatch="////", lw=0, alpha=0.6)
for ax, lab in [(axC, "Power (W)"), (axD, "HR (bpm)")]:
    ax.set_xlim(0, VENTANA_S)
    ax.set_ylabel(lab)
    ax.yaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
axC.tick_params(labelbottom=False)
axD.set_xlabel(f"Time within the detail window (s)")

# leyendas fuera de los datos
leg_top = [Line2D([0], [0], color=TXT2, lw=0.8, ls=(0, (3, 2)), label="Quartile boundary"),
           Patch(color=GRID, label="Detail window (c, d)")]
axA.legend(handles=leg_top, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False,
           handlelength=2.2, columnspacing=2.5)
leg_c = [Patch(facecolor=C_SUB, alpha=0.28, label="Power-increase segment"),
         Patch(facecolor=C_BAJ, alpha=0.28, label="Power-decrease segment")]
axC.legend(handles=leg_c, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2, frameon=False,
           handlelength=2.2, columnspacing=2.5)
leg_d = [Patch(facecolor=C_SUB, alpha=0.28, label="Retained (HR rises)"),
         Patch(facecolor=C_BAJ, alpha=0.28, label="Retained (HR falls)"),
         Patch(facecolor="none", edgecolor=TXT2, hatch="////", label="Discarded (opposite HR direction)")]
axD.legend(handles=leg_d, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3, frameon=False,
           handlelength=2.2, columnspacing=2.0)

for ax, letra in zip([axA, axB, axC, axD], "abcd"):
    ax.text(-0.075, 1.08, letra, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left", color=TXT)

fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.10)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUTPUT_DIR / f"Figura1_metodo_v1.2.{ext}", dpi=300)
print("Escrito:", OUTPUT_DIR / "Figura1_metodo_v1.2.png")

# datos de apoyo para el pie de figura
wc = w[(w.ini >= W0) & (w.fin < W1)]
print(f"Tramos en la ventana: subida {int((wc.tipo=='increase').sum())} "
      f"(conservados {int(((wc.tipo=='increase') & wc.conservado).sum())}), "
      f"bajada {int((wc.tipo=='decrease').sum())} "
      f"(conservados {int(((wc.tipo=='decrease') & wc.conservado).sum())}); "
      f"umbral = {thr:.1f} W/s")
