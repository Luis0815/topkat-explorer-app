import streamlit as st
import pandas as pd
import os
import re

st.title("🔍 TopKAT Explorer – Versión 1")

# Ruta del archivo comprimido dentro del repositorio
csv_path = os.path.join("data", "datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv.gz")

# ====== 1. Cargar el archivo ======
st.header("📂 Carga del archivo")

if not os.path.exists(csv_path):
    st.error(f"❌ Archivo no encontrado: {csv_path}")
    st.info("Asegúrate de tener el archivo comprimido (.csv.gz) en la carpeta /data.")
    st.stop()

else:
    df = pd.read_csv(csv_path, compression="gzip")
    st.success("✅ Archivo cargado correctamente")
    st.write(df.head())


# ====== 2. Clasificación de muestras ======
st.header("🧬 Clasificación automática de muestras")

# Funciones para inferir tipo y Fanconi
def infer_sample_type(name):
    name = name.lower()
    if "carcinoma" in name:
        return "carcinoma"
    elif "dysplasia" in name:
        return "dysplasia"
    elif "stroma_ad_carcinoma" in name:
        return "stroma_ad_carcinoma"
    elif "stroma_ad_dysplasia" in name:
        return "stroma_ad_dysplasia"
    return "unknown"

def infer_fanconi(name):
    return "Fanconi" if re.search(r"f\d|fanconi|_f", name.lower()) else "No_Fanconi"


# ====== 3. Crear tabla resumida ======
df["sample_type"] = df["archivo_base"].apply(infer_sample_type)
df["fanconi_status"] = df["archivo_base"].apply(infer_fanconi)

summary = df.groupby(["sample_type", "fanconi_status"])["PID"].nunique()
summary = summary.reset_index().rename(columns={"PID": "n_muestras"})


st.subheader("📊 Resumen de muestras por categoría")
st.dataframe(summary)


# ====== 4. Selección del usuario ======
st.header("🎛️ Selección de muestras para análisis")

# Tipos disponibles detectados en el CSV
tipos_disponibles = sorted(df["sample_type"].unique())

tipo_seleccionado = st.selectbox("Seleccione el tipo de muestra", tipos_disponibles)

fanconi_opciones = ["Fanconi", "No_Fanconi"]
fanconi_sel = st.selectbox("Seleccione Fanconi / No Fanconi", fanconi_opciones)

# Subconjunto filtrado
subset = df[(df["sample_type"] == tipo_seleccionado) &
            (df["fanconi_status"] == fanconi_sel)]

n_disponibles = subset["PID"].nunique()

st.write(f"🧪 Muestras disponibles para esta selección: **{n_disponibles}**")

n_elegir = st.number_input("¿Cuántas muestras quieres seleccionar aleatoriamente?",
                           min_value=1,
                           max_value=n_disponibles,
                           value=min(5, n_disponibles))

if st.button("Seleccionar muestras"):
    muestras = subset["PID"].drop_duplicates().sample(n_elegir, random_state=42)
    st.success("Muestras seleccionadas:")
    st.write(muestras.tolist())
