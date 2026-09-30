
import os
import numpy as np
import pandas as pd

# 1. CONFIGURACION
#Cambiar a carpeta donde esten las secuencias
base_path = (
    r"C:\Users\Yuyu\Documents\Maestria"
    r"\brca1-variant-embedding-main"
    r"\brca1-variant-embedding-main\Sequences"
)

input_file = os.path.join(
    base_path,
    "Inventario_Tipos_Variantes.xlsx"
)

output_file = os.path.join(
    base_path,
    "BRCA1_SNV_emparejadas_162.xlsx"
)

SEED = 42

# 2. CARGAR Y VALIDAR

df = pd.read_excel(
    input_file,
    sheet_name="Todas_las_variantes"
)

required_columns = {
    "file_name",
    "class",
    "variant_type",
    "start_position",
    "end_position",
    "sequence_length",
    "value_A",
    "value_C",
    "value_G",
    "value_T"
}

missing = required_columns - set(df.columns)

if missing:
    raise ValueError(
        f"Faltan columnas: {sorted(missing)}"
    )

if df["file_name"].str.contains(
    "control",
    case=False,
    na=False
).any():
    raise ValueError(
        "El inventario contiene archivos de control."
    )

# 3. SELECCIONAR SOLO SNV

snv = df[
    df["variant_type"].eq("SNV")
].copy()

print("\n===== SNV DISPONIBLES =====")
print(snv["class"].value_counts())

# 4. EMPAREJAR POR VENTANA GENOMICA

window_columns = [
    "chromosome",
    "start_position",
    "end_position",
    "sequence_length"
]

missing_windows = (
    set(window_columns) - set(snv.columns)
)

if missing_windows:
    raise ValueError(
        "Faltan columnas para identificar ventanas: "
        f"{sorted(missing_windows)}"
    )

selected_groups = []
window_summary = []

for window_key, group in snv.groupby(
    window_columns,
    dropna=False,
    sort=True
):

    benign = group[
        group["class"].eq("benign")
    ]

    pathogenic = group[
        group["class"].eq("pathogenic")
    ]

    n_pairs = min(
        len(benign),
        len(pathogenic)
    )

    if n_pairs == 0:
        continue

    selected_benign = benign.sample(
        n=n_pairs,
        random_state=SEED
    )

    selected_pathogenic = pathogenic.sample(
        n=n_pairs,
        random_state=SEED
    )

    selected_groups.extend([
        selected_benign,
        selected_pathogenic
    ])

    summary = dict(
        zip(window_columns, window_key)
    )

    summary.update({
        "benign_disponibles": len(benign),
        "pathogenic_disponibles": len(pathogenic),
        "pares_seleccionados": n_pairs
    })

    window_summary.append(summary)

if not selected_groups:
    raise ValueError(
        "No se encontraron ventanas con SNV "
        "de ambas clases."
    )

selected = pd.concat(
    selected_groups,
    ignore_index=True
)

summary_df = pd.DataFrame(window_summary)

# 5. VERIFICAR BALANCE

counts = selected["class"].value_counts()

n_benign = counts.get("benign", 0)
n_pathogenic = counts.get("pathogenic", 0)

if n_benign != n_pathogenic:
    raise ValueError(
        "La selección no quedó balanceada."
    )

balance_by_window = pd.crosstab(
    index=[
        selected[col]
        for col in window_columns
    ],
    columns=selected["class"]
)

balance_by_window = balance_by_window.reindex(
    columns=["benign", "pathogenic"],
    fill_value=0
)

if not (
    balance_by_window["benign"]
    == balance_by_window["pathogenic"]
).all():
    raise ValueError(
        "Hay ventanas con cantidades distintas "
        "de benignas y patogenicas."
    )

# 6. PREPARAR DATASET PARA PREPROCESSING

features = [
    "value_A",
    "value_C",
    "value_G",
    "value_T",
    "start_position",
    "sequence_length"
]

dataset_preprocessing = selected[
    features + ["class"]
].copy()

if dataset_preprocessing[features].isna().any().any():
    raise ValueError(
        "Hay valores faltantes en las características."
    )

if not np.isfinite(
    dataset_preprocessing[features].to_numpy(
        dtype=float
    )
).all():
    raise ValueError(
        "Hay valores no finitos en las caracteristicas."
    )

# 7. GUARDAR RESULTADOS

with pd.ExcelWriter(output_file) as writer:

    dataset_preprocessing.to_excel(
        writer,
        sheet_name="Dataset_preprocessing",
        index=False
    )

    selected.to_excel(
        writer,
        sheet_name="Muestras_seleccionadas",
        index=False
    )

    summary_df.to_excel(
        writer,
        sheet_name="Resumen_ventanas",
        index=False
    )

    balance_by_window.reset_index().to_excel(
        writer,
        sheet_name="Verificacion_balance",
        index=False
    )

print("\n===== SELECCION FINALIZADA =====")

print("Benignas seleccionadas:", n_benign)
print("Patogenicas seleccionadas:", n_pathogenic)
print("Total:", len(selected))
print("Ventanas compartidas:", len(summary_df))

print("\nArchivo guardado en:")
print(output_file)

if n_benign != 81:
    print(
        "\nAVISO: Se obtuvieron",
        n_benign,
        "pares, no 81. Revisa el resumen de ventanas."
    )
