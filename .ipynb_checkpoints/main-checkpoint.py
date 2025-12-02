import streamlit as st
import pandas as pd
import os
import re

# ============================================================
#   TITULO
# ============================================================
st.title("🔍 TopKAT Explorer – Comparador de Grupos")


# ============================================================
#   MÓDULO 1 — Cargar datos desde data/raw/
# ============================================================

st.header("📂 1. Carga del archivo principal")

RAW_PATH = os.path.join("data", "raw", "datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv.gz")

if not os.path.exists(RAW_PATH):
    st.error(f"❌ Archivo no encontrado: {RAW_PATH}")
    st.info("Asegúrate de que el archivo comprimido esté en: data/raw/")
    st.stop()

df = pd.read_csv(RAW_PATH, compression="gzip")
st.success("✅ Archivo cargado correctamente")


# ============================================================
#   MÓDULO 2 — Inferir tipo y Fanconi por muestra
# ============================================================

def infer_sample_type(name):
    name = name.lower()
    if "stroma_ad_carcinoma" in name:
        return "stroma_ad_carcinoma"
    if "stroma_ad_dysplasia" in name:
        return "stroma_ad_dysplasia"
    if "carcinoma" in name:
        return "carcinoma"
    if "dysplasia" in name:
        return "dysplasia"
    if "stroma" in name:
        return "stroma"
    return "unknown"

def infer_fanconi(name):
    return "Fanconi" if re.search(r"(?:^|_)f[a-z0-9]", name.lower()) else "No_Fanconi"


df["sample_type"] = df["archivo_base"].apply(infer_sample_type)
df["fanconi_status"] = df["archivo_base"].apply(infer_fanconi)

# Resumen a nivel de PID
sample_summary = (
    df.groupby(["PID", "sample_type", "fanconi_status"])
      .size()
      .reset_index()
      .drop(columns=0)
)

st.header("📊 2. Conteo de muestras por Tipo + Fanconi")

tabla_resumen = (
    sample_summary.groupby(["sample_type", "fanconi_status"])["PID"]
    .nunique()
    .reset_index()
    .rename(columns={"PID": "n_muestras"})
)

st.dataframe(tabla_resumen)


# ============================================================
#   MÓDULO 3 — Comparador de Grupos (reproducible)
# ============================================================

st.header("⚔️ 3. Comparador de Grupos")

tipos_disponibles = sorted(df["sample_type"].unique())
fanconi_opciones = ["Fanconi", "No_Fanconi"]

colA, colB = st.columns(2)

# ---------------------- Grupo A ----------------------
with colA:
    st.subheader("🔵 Grupo A")
    tipo_A = st.selectbox("Tipo", tipos_disponibles, key="tipoA")
    fanconi_A = st.selectbox("Fanconi", fanconi_opciones, key="fanconiA")

# ---------------------- Grupo B ----------------------
with colB:
    st.subheader("🟠 Grupo B")
    tipo_B = st.selectbox("Tipo", tipos_disponibles, key="tipoB")
    fanconi_B = st.selectbox("Fanconi", fanconi_opciones, key="fanconiB")


# Filtrar PIDs reales
grupoA = sample_summary[
    (sample_summary["sample_type"] == tipo_A) &
    (sample_summary["fanconi_status"] == fanconi_A)
]

grupoB = sample_summary[
    (sample_summary["sample_type"] == tipo_B) &
    (sample_summary["fanconi_status"] == fanconi_B)
]

# Cantidades previas
st.write(f"🔵 Grupo A ({tipo_A} - {fanconi_A}): **{len(grupoA)} muestras**")
st.write(f"🟠 Grupo B ({tipo_B} - {fanconi_B}): **{len(grupoB)} muestras**")

if len(grupoA) == 0 or len(grupoB) == 0:
    st.warning("⚠️ Alguno de los grupos no tiene muestras suficientes para comparar.")
    st.stop()


# ============================================================
#   Selección reproducible
# ============================================================

st.header("🎯 4. Selección reproducible de muestras")

nA = st.number_input(
    "¿Cuántas muestras tomar de Grupo A?",
    min_value=1, max_value=len(grupoA), value=min(3, len(grupoA))
)

nB = st.number_input(
    "¿Cuántas muestras tomar de Grupo B?",
    min_value=1, max_value=len(grupoB), value=min(3, len(grupoB))
)

if st.button("Seleccionar muestras"):

    # === REPRODUCIBILIDAD TOTAL ===
    grupoA_sorted = grupoA.sort_values("PID")
    grupoB_sorted = grupoB.sort_values("PID")

    # Semilla fija
    selA = grupoA_sorted["PID"].drop_duplicates().sample(nA, random_state=123).tolist()
    selB = grupoB_sorted["PID"].drop_duplicates().sample(nB, random_state=123).tolist()

    st.success("🎉 Selección completada — Reproducible 100%")

    st.write("🔵 **Grupo A seleccionado:**", selA)
    st.write("🟠 **Grupo B seleccionado:**", selB)

    st.info("Estas muestras pueden enviarse al pipeline de TopKAT, kernels, o comparaciones.")


# ============================================================
#   FIN
# ============================================================
