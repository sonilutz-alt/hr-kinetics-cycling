# =============================================================================
#  REVISION EJAP-D-26-00784 - CELDA 02 - EXTRACCION -> NUEVA BBDD
#  v1.3  - anade la opcion B de la decision 2: pendiente de FC con retraso
#          fisiologico tras cada tramo (sin filtro de signo). Solo columnas extra.
#
#  Hace DOS extracciones de cada archivo con el mismo codigo:
#
#   A) "ANTIGUO"   : port fiel a Python de RESPUESTA_CARDIO.R, con sus fallos.
#                    Solo sirve para VALIDAR que el port reproduce el output
#                    antiguo (hoja 4_Validacion_vs_R). No se usa para analizar.
#
#   B) "CORREGIDO" : mismas definiciones operativas del paper (tramos de
#                    subida/bajada de potencia, umbral 0,3*DE, suavizado k=2,
#                    sigmoide por tramo, cuartiles de tiempo), corrigiendo solo
#                    fallos inequivocos:
#                      1. rejilla de 1 s real: los huecos quedan como NaN
#                      2. derivadas y tramos que no cruzan huecos
#                      3. timestamps duplicados promediados (no dt = 1e-5)
#                      4. sin CP por defecto de 350 W (queda NaN + bandera)
#                      5. posicion, carrera, anyo y corredor desde el nombre
#                      6. sin imputar nada
#                    Las opciones DISCUTIBLES no se deciden aqui: se guardan
#                    en columnas paralelas (ver BLOQUE 1 y hoja 0_Diccionario).
#
#  Salida en  ".../paper 1/revision 1/":
#     BBDD_revision1_v1.3_<fecha>.xlsx   y   PARAMS_extraccion_v1.3_<fecha>.json
# =============================================================================

# ----------------------------------------------------------------------------
#  BLOQUE 1. CONFIGURACION - lo unico que hay que editar
# ----------------------------------------------------------------------------
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent if "__file__" in globals() else Path.cwd()
INPUT_DIR = BASE / "data" / "raw"          # un .xlsx por carrera (exportado del .fit)
OUTPUT_DIR = BASE / "output"
RIDER_MAP = BASE / "data" / "rider_map.csv"  # archivo,corredor_id  (NO se publica)

# Salida del R antiguo SIN imputar, para validar el port (hoja 4).
# Es el archivo "resultados_resumen_por_archivo_cp_quartiles_clean.xlsx".
OLD_SUMMARY = None        # p. ej. Path(".../resultados_resumen_por_archivo_cp_quartiles_clean.xlsx")

# --- Parametros del paper original (NO tocar: definen las variables) ---
POWER_MAX, HR_MIN, HR_MAX = 1800, 30, 220
SMOOTH_K = 2              # media movil del R antiguo
THRESH_SD = 0.3           # umbral = 0,3 * DE de la derivada de potencia
MIN_SEQ_FIT = 6           # filas minimas para ajustar la sigmoide
CP_DUR_ANTIGUO = (60, 300, 720)

# --- PENDIENTE DE DECIDIR JUNTOS (solo afecta a la version CORREGIDA) ---
MAX_INTERP_S = 5          # huecos de FC/potencia <= N s se interpolan; mas largos quedan NaN
                          # (el R antiguo interpolaba huecos de CUALQUIER longitud)

# --- CRITERIOS DE INCLUSION (decision 1, acordada 23-09-2026) ---
EXCLUIR_POTENCIA_IMPOSIBLE = True   # cualquier muestra cruda > POWER_MAX W -> archivo fuera
NIVEL_COBERTURA = "ninguno"         # la cobertura NO excluye. Opciones para sensibilidad:
                                    # "ninguno" | "archivo" | "cuartil"
MIN_COBERTURA = 80                  # % usado solo para MARCAR (y para sensibilidad)

# --- DECISION 2, OPCION B (en evaluacion, no decidida) ---
# Pendiente de la FC suavizada entre a y b segundos DESPUES del inicio de cada
# tramo de subida/bajada de potencia: (FC[inicio+b] - FC[inicio+a]) / (b - a).
# Sin filtro de signo. Exige FC valida en toda la ventana.
VENTANAS_RETRASO = [(5, 15), (10, 20), (5, 20)]
MIN_TRAMO_RETRASO = [1, 6]          # todos los tramos / solo tramos >= 6 s

# ----------------------------------------------------------------------------
#  BLOQUE 2. Imports
# ----------------------------------------------------------------------------
import re
import json
import warnings
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from tqdm.auto import tqdm

warnings.filterwarnings("ignore")
SELLO = datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
Q = ["Q1", "Q2", "Q3", "Q4"]

files = sorted([p for p in INPUT_DIR.iterdir()
                if p.suffix.lower() in (".xlsx", ".xls") and not p.name.startswith("~$")])
print(f"Carpeta : {INPUT_DIR}\nSalida  : {OUTPUT_DIR}\nArchivos: {len(files)}")

# ----------------------------------------------------------------------------
#  BLOQUE 3. Identificacion del archivo (corredor, carrera, anyo, posicion)
#  El corredor se asigna con esta tabla (41 archivos del paper). Las grafias
#  distintas del mismo corredor se unifican aqui: REVISAR y CONFIRMAR.
#  Cualquier archivo que no este en la tabla sale marcado "REVISAR".
# ----------------------------------------------------------------------------
# Identificador anonimo del corredor por archivo. Se lee de un CSV local que
# NO se publica (confidencialidad de los deportistas): columnas archivo,corredor_id
CORREDOR = (pd.read_csv(RIDER_MAP).set_index("archivo")["corredor_id"].astype(str).to_dict()
            if RIDER_MAP.exists() else {})


def identifica(nombre):
    base = Path(nombre).stem
    low = base.lower()
    pos = re.search(r"[NP]\s*(\d{1,3})\s*$", base)
    carrera = ("LBL" if re.search(r"lvl|lieja|liege|lbl", low) else
               "MSR" if re.search(r"msr|san ?remo|sanremo", low) else
               "RVV" if re.search(r"ronde|flande|rvv", low) else
               "PR" if re.search(r"roubaix|\bpr\s?\d", low) else
               "LOM" if re.search(r"lombard", low) else None)
    m = re.match(r"^(\d{2})_", base) or re.search(
        r"(?:deronde|pr|lieja|msr|roubaix|san ?remo|flande|lombardia)\s?(\d{2})", low)
    corredor = CORREDOR.get(nombre)
    return {"archivo": nombre,
            "corredor": corredor if corredor else "REVISAR",
            "carrera": carrera, "anyo": 2000 + int(m.group(1)) if m else None,
            "posicion": int(pos.group(1)) if pos else None,
            "id_revisar": corredor is None or carrera is None or pos is None or m is None}


# ----------------------------------------------------------------------------
#  BLOQUE 4. Funciones comunes
# ----------------------------------------------------------------------------
R_PATTERNS = {"timestamp": r"timestamp", "power": r"power",
              "heart_rate": r"heart_rate|hr|bpm|pulse|pulse_rate|fc|heart|rate|beat"}


def elige_columnas(cols, modo):
    """ANTIGUO: primera coincidencia de grepl (como el R).
       CORREGIDO: nombre exacto si existe; si no, como el R."""
    out = {}
    for k, pat in R_PATTERNS.items():
        r_col = next((c for c in cols if re.search(pat, str(c).lower())), None)
        exacta = next((c for c in cols if str(c).strip().lower() == k), None)
        out[k] = r_col if modo == "antiguo" else (exacta or r_col)
    return out


def a_segundos(serie):
    if np.issubdtype(serie.dtype, np.datetime64):
        s = pd.to_datetime(serie, errors="coerce")
        return (s - s.min()).dt.total_seconds().values
    num = pd.to_numeric(serie, errors="coerce")
    if num.notna().mean() > 0.95:
        return (num - num.min()).values.astype(float)
    s = pd.to_datetime(serie, errors="coerce", utc=True)
    return (s - s.min()).dt.total_seconds().values


def rollmean_k2(x):
    """zoo::rollmean(k=2, align='center', fill='extend'): (x[i]+x[i+1])/2,
       el ultimo repite el anterior. NaN se propaga (no puentea huecos)."""
    x = np.asarray(x, float)
    if len(x) < 2:
        return x.copy()
    y = np.empty_like(x)
    y[:-1] = (x[:-1] + x[1:]) / 2
    y[-1] = y[-2]
    return y


def na_approx(x, max_gap=None):
    """zoo::na.approx(na.rm=FALSE): interpola NaN interiores. Con max_gap solo
       rellena rachas de NaN de longitud <= max_gap (en muestras)."""
    s = pd.Series(np.asarray(x, float))
    if max_gap is None:
        return s.interpolate(limit_area="inside").values
    isn = s.isna()
    run_id = (isn != isn.shift()).cumsum()
    run_len = isn.groupby(run_id).transform("sum")
    filled = s.interpolate(limit_area="inside")
    filled[isn & (run_len > max_gap)] = np.nan
    return filled.values


def find_sequences(v, thr, direction):
    """Port exacto de find_sequences() del R (indices 0-based, fin incluido:
       el fin es la fila que ROMPE la racha, igual que c(start, i) en R)."""
    seqs, start = [], None
    n = len(v)
    for i in range(n):
        val = v[i]
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
        seqs.append((start, n - 1))
    return seqs


def sigmoid(x, A, B, C, D):
    return A + B / (1 + np.exp(-C * (x - D)))


def rhri_sigmoide(time, fc):
    """Ajuste con limites como nls(algorithm='port') y rHRI = pendiente maxima
       sobre 1000 puntos. Devuelve NaN si falla o si sale negativa."""
    ok = ~(np.isnan(time) | np.isnan(fc))
    time, fc = time[ok], fc[ok]
    if len(fc) < MIN_SEQ_FIT:
        return np.nan
    lo_A, hi_A = fc.min(), fc.max()
    if hi_A <= lo_A:
        hi_A = lo_A + 1e-9
    rng = fc.max() - fc.min()
    lower = [lo_A, 0, -1, time.min()]
    upper = [hi_A, rng + 1, 1, time.max() if time.max() > time.min() else time.min() + 1e-9]
    p0 = [fc.min(), max(rng, 1e-12), 0.1, time.mean()]
    p0 = [min(max(p, l), u) for p, l, u in zip(p0, lower, upper)]
    try:
        par, _ = curve_fit(sigmoid, time, fc, p0=p0, bounds=(lower, upper), maxfev=2000)
    except Exception:
        return np.nan
    x = np.linspace(time.min(), time.max(), 1000)
    y = sigmoid(x, *par)
    r = np.nanmax(np.diff(y) / np.diff(x))
    return r if r >= 0 else np.nan


def mmp_rejilla(p, w):
    s = pd.Series(p).rolling(w, min_periods=w).mean()
    return float(s.max()) if s.notna().any() else np.nan


def cp_formula_antigua(mmps):
    """lm(1/P ~ 1/t): CP = 1/intercepto (formula del R antiguo)."""
    t = np.array(list(mmps.keys()), float); P = np.array(list(mmps.values()), float)
    ok = np.isfinite(P) & (P > 0)
    if ok.sum() < 2:
        return np.nan
    b, a = np.polyfit(1 / t[ok], 1 / P[ok], 1)
    return 1 / a if a != 0 else np.nan


def cp_lineal(mmps):
    """Modelo lineal P = CP + W'/t. Devuelve (CP, W')."""
    t = np.array(list(mmps.keys()), float); P = np.array(list(mmps.values()), float)
    ok = np.isfinite(P)
    if ok.sum() < 2:
        return np.nan, np.nan
    b, a = np.polyfit(1 / t[ok], P[ok], 1)
    return a, b


# ----------------------------------------------------------------------------
#  BLOQUE 5. Resumen por cuartil a partir de los tramos (comun a ambas)
# ----------------------------------------------------------------------------
def pendiente_retraso(fc_s, s, a, b):
    """(FC[s+b] - FC[s+a]) / (b-a) si toda la ventana tiene FC valida."""
    if s + b >= len(fc_s):
        return np.nan
    w = fc_s[s + a: s + b + 1]
    if np.any(np.isnan(w)):
        return np.nan
    return (fc_s[s + b] - fc_s[s + a]) / (b - a)


def resume_tramos(t, fc_s, p_s, fc_der, qlab, thr_series, con_retraso=False):
    """Devuelve filas por tramo (subida/bajada) con las metricas del R antiguo
       y ademas la derivada media SIN el filtro de signo. Con con_retraso=True
       anade las pendientes de FC con retraso (opcion B)."""
    filas = []
    for tipo in ("increase", "decrease"):
        for s, e in find_sequences(thr_series[0], thr_series[1], tipo):
            sl = slice(s, e + 1)
            fd = np.nanmean(fc_der[sl]) if np.any(~np.isnan(fc_der[sl])) else np.nan
            fd_filtrada = fd
            if not np.isnan(fd):
                if tipo == "increase" and fd < 0:
                    fd_filtrada = np.nan
                if tipo == "decrease" and fd > 0:
                    fd_filtrada = np.nan
            fila = {"tipo": tipo, "cuartil": qlab[s], "n_filas": e - s + 1,
                    "fc_deriv_media": fd_filtrada,
                    "fc_deriv_media_sin_filtro_signo": fd,
                    "pot_media": np.nanmean(p_s[sl]) if np.any(~np.isnan(p_s[sl])) else np.nan,
                    "fc_max": np.nanmax(fc_s[sl]) if np.any(~np.isnan(fc_s[sl])) else np.nan,
                    "rHRI": np.nan}
            if tipo == "increase" and (e - s + 1) >= MIN_SEQ_FIT and np.any(~np.isnan(fc_s[sl])):
                fila["rHRI"] = rhri_sigmoide(t[sl], fc_s[sl])
            if con_retraso:
                for a, b in VENTANAS_RETRASO:
                    fila[f"ret_{a}_{b}"] = pendiente_retraso(fc_s, s, a, b)
            filas.append(fila)
    return pd.DataFrame(filas)


def media_o_nan(x):
    x = pd.Series(x).dropna()
    return x.mean() if len(x) else np.nan


def max_o_nan(x):
    x = pd.Series(x).dropna()
    return x.max() if len(x) else np.nan


# ----------------------------------------------------------------------------
#  BLOQUE 6. Extraccion ANTIGUA (port fiel del R, solo para validar)
# ----------------------------------------------------------------------------
def extrae_antiguo(d, cols):
    t_raw = a_segundos(d[cols["timestamp"]])
    p = pd.to_numeric(d[cols["power"]], errors="coerce").values.astype(float)
    h = pd.to_numeric(d[cols["heart_rate"]], errors="coerce").values.astype(float)
    p[(p > POWER_MAX) | (p < 0)] = np.nan
    h[(h < HR_MIN) | (h > HR_MAX)] = np.nan
    p = na_approx(p); h = na_approx(h)
    if len(d) < 100 or np.sum(~np.isnan(p)) < 100 or np.sum(~np.isnan(h)) < 100:
        return None, "Sparse Data"
    # En el R: approx(as.numeric(timestamp), n = nrow(data)) interpola en los
    # indices 1..n y devuelve los MISMOS timestamps: no cambia nada. Los huecos
    # se conservan y la derivada a traves de un hueco es dFC / duracion_hueco.
    t = t_raw - np.nanmin(t_raw)
    # CP con rollmean por FILAS (como el R) y 350 W por defecto
    mm = {}
    for w in CP_DUR_ANTIGUO:
        mm[w] = (pd.Series(p).rolling(w, center=True).mean().max()
                 if len(p) >= w else 0.0)
    cp = cp_formula_antigua(mm)
    if np.isnan(cp) or cp < 200 or cp > 500:
        cp = 350.0
    ps = rollmean_k2(p); hs = rollmean_k2(h)
    ok = ~(np.isnan(ps) | np.isnan(hs))
    t, ps, hs = t[ok], ps[ok], hs[ok]
    pct = ps / cp * 100
    dt = np.diff(t); dt[dt == 0] = 1e-5
    pct_der = np.r_[np.nan, np.diff(pct) / dt]
    fc_der = np.r_[np.nan, np.diff(hs) / dt]
    sd = np.nanstd(pct_der, ddof=1)
    thr = 0.3 * sd if (sd and not np.isnan(sd)) else 1e-5
    edges = np.linspace(0, t.max(), 5)
    qlab = np.array(Q)[np.clip(np.searchsorted(edges, t, side="left") - 1, 0, 3)]
    tr = resume_tramos(t, hs, pct, fc_der, qlab, (pct_der, thr))
    out = {"CP_antiguo_W": cp}
    for q in Q:
        up = tr[(tr.tipo == "increase") & (tr.cuartil == q)]
        dn = tr[(tr.tipo == "decrease") & (tr.cuartil == q)]
        n = q[-1]
        out[f"avg_rHRI_increase_bpm_per_s_q{n}"] = media_o_nan(up.rHRI)
        out[f"avg_fc_deriv_increase_bpm_per_s_q{n}"] = media_o_nan(up.fc_deriv_media)
        out[f"avg_fc_deriv_decrease_bpm_per_s_q{n}"] = media_o_nan(dn.fc_deriv_media)
        out[f"avg_power_increase_percent_cp_q{n}"] = media_o_nan(up.pot_media)
        out[f"max_fc_increase_bpm_q{n}"] = max_o_nan(up.fc_max)
    return out, "Success"


# ----------------------------------------------------------------------------
#  BLOQUE 7. Extraccion CORREGIDA
# ----------------------------------------------------------------------------
def extrae_corregido(d, cols):
    t_raw = a_segundos(d[cols["timestamp"]])
    p_raw = pd.to_numeric(d[cols["power"]], errors="coerce").values.astype(float)
    h_raw = pd.to_numeric(d[cols["heart_rate"]], errors="coerce").values.astype(float)
    ok_t = ~np.isnan(t_raw)
    seg = np.round(t_raw[ok_t]).astype(int)
    # rejilla de 1 s; duplicados promediados; segundos ausentes = NaN
    df = pd.DataFrame({"s": seg, "p": p_raw[ok_t], "h": h_raw[ok_t]})
    n_dup = int(df.duplicated("s").sum())
    n_pot_imposible = int(((df.p > POWER_MAX) | (df.p < 0)).sum())
    df.loc[(df.p > POWER_MAX) | (df.p < 0), "p"] = np.nan
    df.loc[(df.h < HR_MIN) | (df.h > HR_MAX), "h"] = np.nan
    g = df.groupby("s")[["p", "h"]].mean()
    grid = g.reindex(np.arange(g.index.min(), g.index.max() + 1))
    t = (grid.index.values - grid.index.min()).astype(float)
    p = na_approx(grid.p.values, MAX_INTERP_S)
    h = na_approx(grid.h.values, MAX_INTERP_S)

    info = {"n_segundos_rejilla": len(grid), "n_timestamps_duplicados": n_dup,
            "n_muestras_pot_>1800_o_<0": n_pot_imposible,
            "n_seg_sin_registro": int(len(grid) - len(g)),
            "pct_fc_valida_tras_interp": 100 * np.mean(~np.isnan(h)),
            "pct_pot_valida_tras_interp": 100 * np.mean(~np.isnan(p))}

    # --- CP: varias opciones, SIN valor por defecto ---
    mm = {w: mmp_rejilla(p, w) for w in (60, 180, 300, 600, 720, 1200)}
    info.update({f"MMP_{w // 60}min_W": v for w, v in mm.items()})
    info["CP_formula_antigua_1_5_12_W"] = cp_formula_antigua({w: mm[w] for w in (60, 300, 720)})
    info["CP_lineal_1_5_12_W"], info["Wprima_lineal_1_5_12_J"] = cp_lineal({w: mm[w] for w in (60, 300, 720)})
    info["CP_lineal_3_5_12_W"], info["Wprima_lineal_3_5_12_J"] = cp_lineal({w: mm[w] for w in (180, 300, 720)})
    cpa = info["CP_formula_antigua_1_5_12_W"]
    info["flag_CP_antigua_fuera_200_500"] = bool(np.isnan(cpa) or cpa < 200 or cpa > 500)

    # --- suavizado, derivadas (dt = 1 s; NaN no se puentea) ---
    ps = rollmean_k2(p); hs = rollmean_k2(h)
    p_der = np.r_[np.nan, np.diff(ps)]        # W/s; umbral equivalente al de %CP
    fc_der = np.r_[np.nan, np.diff(hs)]
    sd = np.nanstd(p_der, ddof=1)
    thr = 0.3 * sd if (sd and not np.isnan(sd)) else 1e-5
    edges = np.linspace(0, t.max(), 5)
    qlab = np.array(Q)[np.clip(np.searchsorted(edges, t, side="left") - 1, 0, 3)]
    tr = resume_tramos(t, hs, ps, fc_der, qlab, (p_der, thr), con_retraso=True)

    filas = []
    for q in Q:
        m = qlab == q
        up = tr[(tr.tipo == "increase") & (tr.cuartil == q)]
        dn = tr[(tr.tipo == "decrease") & (tr.cuartil == q)]
        filas.append({
            "cuartil": q,
            "seg_cuartil": int(m.sum()),
            "pct_fc_valida": 100 * np.mean(~np.isnan(h[m])) if m.any() else np.nan,
            "pct_pot_valida": 100 * np.mean(~np.isnan(p[m])) if m.any() else np.nan,
            "n_tramos_subida": len(up), "n_tramos_bajada": len(dn),
            "n_tramos_ajuste_intentado": int((up.n_filas >= MIN_SEQ_FIT).sum()) if len(up) else 0,
            "n_tramos_ajuste_ok": int(up.rHRI.notna().sum()) if len(up) else 0,
            # --- definiciones del paper (con filtro de signo, como el R) ---
            "rHRI_media_bpm_s": media_o_nan(up.rHRI) if len(up) else np.nan,
            "dFC_subida_bpm_s": media_o_nan(up.fc_deriv_media) if len(up) else np.nan,
            "dFC_bajada_bpm_s": media_o_nan(dn.fc_deriv_media) if len(dn) else np.nan,
            "pot_tramos_subida_W": media_o_nan(up.pot_media) if len(up) else np.nan,
            "FCmax_tramos_subida_bpm": max_o_nan(up.fc_max) if len(up) else np.nan,
            # --- alternativas para decidir despues ---
            "dFC_subida_sin_filtro_bpm_s": media_o_nan(up.fc_deriv_media_sin_filtro_signo) if len(up) else np.nan,
            "dFC_bajada_sin_filtro_bpm_s": media_o_nan(dn.fc_deriv_media_sin_filtro_signo) if len(dn) else np.nan,
            "pot_media_cuartil_W": float(np.nanmean(p[m])) if m.any() and np.any(~np.isnan(p[m])) else np.nan,
            "FCmax_cuartil_bpm": float(np.nanmax(hs[m])) if m.any() and np.any(~np.isnan(hs[m])) else np.nan,
            "FCmedia_cuartil_bpm": float(np.nanmean(h[m])) if m.any() and np.any(~np.isnan(h[m])) else np.nan,
        })
        # --- opcion B: pendiente de FC con retraso tras los tramos ---
        for a, b in VENTANAS_RETRASO:
            for lmin in MIN_TRAMO_RETRASO:
                suf = f"{a}_{b}s" + ("" if lmin <= 1 else f"_tramos{lmin}s")
                u = up[up.n_filas >= lmin] if len(up) else up
                dd = dn[dn.n_filas >= lmin] if len(dn) else dn
                vu = media_o_nan(u[f"ret_{a}_{b}"]) if len(u) else np.nan
                vd = media_o_nan(dd[f"ret_{a}_{b}"]) if len(dd) else np.nan
                filas[-1][f"B_dFC_tras_subida_{suf}"] = vu
                filas[-1][f"B_dFC_tras_bajada_{suf}"] = vd
                filas[-1][f"B_dFC_neta_{suf}"] = vu - vd if not (np.isnan(vu) or np.isnan(vd)) else np.nan
                filas[-1][f"B_n_tramos_subida_{suf}"] = int(u[f"ret_{a}_{b}"].notna().sum()) if len(u) else 0
    return info, filas


# ----------------------------------------------------------------------------
#  BLOQUE 8. Recorrido de archivos
# ----------------------------------------------------------------------------
antiguo, archivo_info, larga, errores = [], [], [], []

for f in tqdm(files, desc="Extraccion"):
    ident = identifica(f.name)
    try:
        xl = pd.ExcelFile(f)
        d = xl.parse("Sheet1" if "Sheet1" in xl.sheet_names else xl.sheet_names[0])
        cols_a = elige_columnas(list(d.columns), "antiguo")
        cols_c = elige_columnas(list(d.columns), "corregido")
        if any(v is None for v in cols_c.values()):
            raise ValueError(f"Columnas no encontradas: {cols_c}")
    except Exception as e:
        errores.append({"archivo": f.name, "fase": "lectura", "error": str(e)})
        continue

    # A) antiguo
    try:
        res_a, estado = extrae_antiguo(d, cols_a)
        antiguo.append({"archivo": f.name, "estado": estado,
                        "cols_R": json.dumps(cols_a, ensure_ascii=False), **(res_a or {})})
    except Exception as e:
        errores.append({"archivo": f.name, "fase": "antiguo", "error": str(e)})

    # B) corregido
    try:
        info, filas = extrae_corregido(d, cols_c)
        info.update(ident)
        info["cols_usadas"] = json.dumps(cols_c, ensure_ascii=False)
        info["cols_R_distintas"] = cols_a != cols_c
        archivo_info.append(info)
        for r in filas:
            larga.append({**ident, **r})
    except Exception as e:
        errores.append({"archivo": f.name, "fase": "corregido", "error": str(e)})

ant = pd.DataFrame(antiguo)
arch = pd.DataFrame(archivo_info)
lar = pd.DataFrame(larga)
err = pd.DataFrame(errores)

# orden de columnas de la hoja por archivo
id_cols = ["archivo", "corredor", "carrera", "anyo", "posicion", "id_revisar"]
if len(arch):
    arch = arch[id_cols + [c for c in arch.columns if c not in id_cols]]

# ----------------------------------------------------------------------------
#  BLOQUE 8b. Criterios de inclusion (no se borra nada: se marca)
# ----------------------------------------------------------------------------
if len(arch):
    qmin = lar.groupby("archivo")[["pct_pot_valida", "pct_fc_valida"]].min()
    qmin.columns = ["pct_pot_valida_peor_cuartil", "pct_fc_valida_peor_cuartil"]
    arch = arch.merge(qmin, on="archivo", how="left")
    arch["excl_cobertura_archivo"] = ((arch.pct_pot_valida_tras_interp < MIN_COBERTURA) |
                                      (arch.pct_fc_valida_tras_interp < MIN_COBERTURA))
    arch["excl_cobertura_cuartil"] = ((arch.pct_pot_valida_peor_cuartil < MIN_COBERTURA) |
                                      (arch.pct_fc_valida_peor_cuartil < MIN_COBERTURA))
    arch["excl_potencia_imposible"] = (EXCLUIR_POTENCIA_IMPOSIBLE &
                                       (arch["n_muestras_pot_>1800_o_<0"] > 0))
    if NIVEL_COBERTURA in ("archivo", "cuartil"):
        excl_cob = arch[f"excl_cobertura_{NIVEL_COBERTURA}"].astype(bool)
    else:
        excl_cob = pd.Series(False, index=arch.index)

    def motivo(i):
        m = []
        if arch.loc[i, "excl_potencia_imposible"]:
            m.append("potencia imposible")
        if excl_cob.loc[i]:
            m.append(f"cobertura < {MIN_COBERTURA}% ({NIVEL_COBERTURA})")
        return "; ".join(m)
    arch["motivo_exclusion"] = [motivo(i) for i in arch.index]
    arch["incluido"] = arch.motivo_exclusion == ""
    # marca por cuartil (no excluye): util para analisis de sensibilidad
    lar["cobertura_baja_cuartil"] = ((lar.pct_pot_valida < MIN_COBERTURA) |
                                     (lar.pct_fc_valida < MIN_COBERTURA))
    lar = lar.merge(arch[["archivo", "incluido"]], on="archivo", how="left")

    fin = arch[arch.incluido]
    flujo = [{"paso": "Archivos en la carpeta", "archivos": len(files), "corredores": np.nan},
             {"paso": "Sin columna de FC", "archivos": -int((err.fase == "lectura").sum()) if len(err) else 0, "corredores": np.nan},
             {"paso": "Con FC y procesados", "archivos": len(arch), "corredores": arch.corredor.nunique()},
             {"paso": "Excluidos por potencia imposible", "archivos": -int(arch.excl_potencia_imposible.sum()), "corredores": np.nan}]
    if NIVEL_COBERTURA in ("archivo", "cuartil"):
        flujo.append({"paso": f"Excluidos por cobertura < {MIN_COBERTURA}% ({NIVEL_COBERTURA})",
                      "archivos": -int((excl_cob & ~arch.excl_potencia_imposible).sum()), "corredores": np.nan})
    flujo += [{"paso": "MUESTRA FINAL", "archivos": len(fin), "corredores": fin.corredor.nunique()},
              {"paso": f"  de ellos con cobertura < {MIN_COBERTURA}% en algun cuartil (marcados)",
               "archivos": int(fin.excl_cobertura_cuartil.sum()), "corredores": np.nan},
              {"paso": f"  cuartiles marcados (de {4 * len(fin)})",
               "archivos": int(lar[lar.archivo.isin(fin.archivo)].cobertura_baja_cuartil.sum()), "corredores": np.nan}]
    flujo = pd.DataFrame(flujo)
else:
    flujo = pd.DataFrame()

# ----------------------------------------------------------------------------
#  BLOQUE 9. Formato ancho (una fila por archivo) de la version corregida
# ----------------------------------------------------------------------------
metricas = [c for c in lar.columns if c not in id_cols + ["cuartil", "incluido"]] if len(lar) else []
ancha = (lar.pivot_table(index="archivo", columns="cuartil", values=metricas, aggfunc="first")
         if len(lar) else pd.DataFrame())
if len(ancha):
    ancha.columns = [f"{m}_{q}" for m, q in ancha.columns]
    ancha = arch[id_cols + ["incluido", "motivo_exclusion"]].merge(
        ancha.reset_index(), on="archivo", how="left")

# ----------------------------------------------------------------------------
#  BLOQUE 10. Validacion del port frente al R antiguo
# ----------------------------------------------------------------------------
val, val_det = pd.DataFrame(), pd.DataFrame()
if OLD_SUMMARY is not None and Path(OLD_SUMMARY).exists() and len(ant):
    old = pd.read_excel(OLD_SUMMARY).rename(
        columns={"file": "archivo", "critical_power_W": "CP_antiguo_W"})
    mm = old.merge(ant, on="archivo", how="outer", suffixes=("_R", "_py"), indicator=True)
    filas, det = [], []
    for c in ["CP_antiguo_W"] + [c for c in ant.columns if c.startswith(("avg_", "max_"))]:
        if f"{c}_R" not in mm or f"{c}_py" not in mm:
            continue
        a_, b_ = mm[f"{c}_R"].astype(float), mm[f"{c}_py"].astype(float)
        ok = a_.notna() & b_.notna()
        filas.append({"variable": c, "n_ambos": int(ok.sum()),
                      "NaN_solo_R": int((a_.isna() & b_.notna()).sum()),
                      "NaN_solo_py": int((a_.notna() & b_.isna()).sum()),
                      "r_pearson": a_[ok].corr(b_[ok]) if ok.sum() > 2 else np.nan,
                      "dif_abs_max": (a_[ok] - b_[ok]).abs().max() if ok.any() else np.nan,
                      "dif_abs_media": (a_[ok] - b_[ok]).abs().mean() if ok.any() else np.nan})
        det.append(pd.DataFrame({"archivo": mm.archivo, "variable": c,
                                 "R": a_, "python": b_, "dif": b_ - a_}))
    val = pd.DataFrame(filas)
    val_det = pd.concat(det, ignore_index=True) if det else pd.DataFrame()
    solo = mm.loc[mm["_merge"] != "both", ["archivo", "_merge"]]
    if len(solo):
        print("Archivos que no estan en ambos lados:\n", solo.to_string(index=False))

# ----------------------------------------------------------------------------
#  BLOQUE 11. Diccionario de variables
# ----------------------------------------------------------------------------
dicc = pd.DataFrame([
    ("rHRI_media_bpm_s", "Media, en el cuartil, de la pendiente maxima de una sigmoide ajustada a la FC suavizada en cada tramo de subida de potencia de >= 6 s (definicion del paper)"),
    ("dFC_subida_bpm_s", "Media de las derivadas medias de FC en los tramos de subida de potencia; se descartan tramos con media negativa (filtro de signo del paper)"),
    ("dFC_bajada_bpm_s", "Idem en tramos de bajada de potencia; se descartan tramos con media positiva"),
    ("dFC_subida_sin_filtro_bpm_s", "ALTERNATIVA: como dFC_subida pero sin descartar tramos por signo"),
    ("dFC_bajada_sin_filtro_bpm_s", "ALTERNATIVA: como dFC_bajada pero sin descartar tramos por signo"),
    ("pot_tramos_subida_W", "Potencia media en los tramos de subida (el %CP del paper = esto / CP * 100)"),
    ("pot_media_cuartil_W", "ALTERNATIVA: potencia media de todo el cuartil (lo que dice el texto de Metodos)"),
    ("FCmax_tramos_subida_bpm", "FC maxima suavizada dentro de los tramos de subida (definicion del codigo)"),
    ("FCmax_cuartil_bpm", "ALTERNATIVA: FC maxima suavizada de todo el cuartil (lo que dice Metodos)"),
    ("FCmedia_cuartil_bpm", "FC media del cuartil (para responder al comentario de la deriva)"),
    ("CP_*", "Opciones de CP desde la propia carrera; ninguna se impone todavia. Sin 350 W por defecto"),
    ("pct_fc_valida / pct_pot_valida", "Cobertura del cuartil tras interpolar solo huecos <= MAX_INTERP_S"),
    ("n_tramos_*", "Numero de tramos detectados y de ajustes de sigmoide intentados / correctos"),
    ("B_dFC_tras_subida_a_bs", "OPCION B: media de (FC[inicio+b]-FC[inicio+a])/(b-a) tras cada tramo de SUBIDA de potencia; sin filtro de signo"),
    ("B_dFC_tras_bajada_a_bs", "OPCION B: idem tras cada tramo de BAJADA de potencia"),
    ("B_dFC_neta_a_bs", "OPCION B: subida menos bajada (respuesta neta, descuenta la deriva de la FC)"),
    ("*_tramos6s", "OPCION B restringida a tramos de >= 6 s"),
], columns=["variable", "definicion"])

# ----------------------------------------------------------------------------
#  BLOQUE 12. Exportacion
# ----------------------------------------------------------------------------
ruta = OUTPUT_DIR / f"BBDD_revision1_v1.3_{SELLO}.xlsx"
with pd.ExcelWriter(ruta, engine="openpyxl") as w:
    dicc.to_excel(w, sheet_name="0_Diccionario", index=False)
    lar.to_excel(w, sheet_name="1_BBDD_larga", index=False)
    ancha.to_excel(w, sheet_name="2_BBDD_ancha", index=False)
    arch.to_excel(w, sheet_name="3_Por_archivo_CP_calidad", index=False)
    flujo.to_excel(w, sheet_name="3b_Flujo_muestra", index=False)
    val.to_excel(w, sheet_name="4_Validacion_vs_R", index=False)
    val_det.to_excel(w, sheet_name="4b_Validacion_detalle", index=False)
    ant.to_excel(w, sheet_name="5_Port_antiguo", index=False)
    err.to_excel(w, sheet_name="6_Errores", index=False)

params = {k: (str(v) if isinstance(v, (Path, tuple)) else v) for k, v in {
    "INPUT_DIR": INPUT_DIR, "POWER_MAX": POWER_MAX, "HR_MIN": HR_MIN, "HR_MAX": HR_MAX,
    "SMOOTH_K": SMOOTH_K, "THRESH_SD": THRESH_SD, "MIN_SEQ_FIT": MIN_SEQ_FIT,
    "CP_DUR_ANTIGUO": CP_DUR_ANTIGUO, "MAX_INTERP_S": MAX_INTERP_S,
    "MIN_COBERTURA": MIN_COBERTURA, "NIVEL_COBERTURA": NIVEL_COBERTURA,
    "EXCLUIR_POTENCIA_IMPOSIBLE": EXCLUIR_POTENCIA_IMPOSIBLE,
    "VENTANAS_RETRASO": str(VENTANAS_RETRASO), "MIN_TRAMO_RETRASO": str(MIN_TRAMO_RETRASO),
    "n_archivos": len(files), "n_ok": len(arch), "n_errores": len(err), "fecha": SELLO}.items()}
(OUTPUT_DIR / f"PARAMS_extraccion_v1.3_{SELLO}.json").write_text(
    json.dumps(params, indent=2, ensure_ascii=False))

print(f"\nArchivos procesados : {len(arch)} de {len(files)}  | errores: {len(err)}")
if len(arch):
    print(f"Corredores unicos   : {arch.corredor.nunique()}  | a revisar: {int(arch.id_revisar.sum())}")
    print(f"Carreras            : {arch.carrera.value_counts().to_dict()}")
if len(flujo):
    print("\nFlujo de la muestra:")
    print(flujo.to_string(index=False))
if len(val):
    print("\nValidacion del port frente al R antiguo:")
    print(val.round(3).to_string(index=False))
print(f"\nEscrito: {ruta}")
