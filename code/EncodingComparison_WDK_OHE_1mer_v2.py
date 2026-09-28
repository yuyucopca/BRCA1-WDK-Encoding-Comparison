import os
import re
import ast
import time
import numpy as np
import pandas as pd

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import accuracy_score, f1_score, cohen_kappa_score


#1.Configuration

base_path = (
    r"C:\Users\Yuyu\Documents\Maestria"
    r"\brca1-variant-embedding-main"
    r"\brca1-variant-embedding-main\Sequences"
)

matched_file = os.path.join(
    base_path,
    "BRCA1_SNV_emparejadas_162.xlsx"
)

path_benign = os.path.join(base_path, "SeqBRCA1_Benigns")
path_pathogenic = os.path.join(base_path, "SeqBRCA1_Pathogenics")

output_file = os.path.join(
    base_path,
    "Encoding_Comparison_WDK_OHE_1mer.xlsx"
)

SEED = 42
N_SPLITS = 5


# 2. Sequencer reader
#    Based in EncodingNucleotidsV2.py

def read_fasta(filepath):
    with open(filepath, "r", encoding="utf-8-sig") as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError(f"Archivo vacao: {filepath}")

    sequences = []
    header = None
    content_lines = []

    def save_sequence(current_header, current_lines):
        if current_header is None:
            raise ValueError(
                f"Secuencia sin encabezado FASTA en {filepath}"
            )

        if not current_lines:
            raise ValueError(
                f"Encabezado sin secuencia en {filepath}"
            )

        content = "\n".join(current_lines).strip()

        if content.startswith("["):
            try:
                parsed = ast.literal_eval(content)
            except (ValueError, SyntaxError) as error:
                raise ValueError(
                    f"No se pudo interpretar la lista en {filepath}"
                ) from error

            if not isinstance(parsed, list):
                raise ValueError(
                    f"El contenido no es una lista en {filepath}"
                )

            if not all(
                isinstance(nuc, str) and len(nuc) == 1
                for nuc in parsed
            ):
                raise ValueError(
                    f"Lista de nucleótidos inválida en {filepath}"
                )

            sequence = "".join(parsed).upper()

        else:
            valid_lines = []

            for line in current_lines:
                clean_line = re.sub(
                    r"\s*\(\s*\d+\s*-\s*\d+\s*\)\s*$",
                    "",
                    line
                )

                clean_line = "".join(
                    clean_line.split()
                ).upper()

                if not clean_line:
                    continue

                invalid = set(clean_line) - set("ACGT")

                if invalid:
                    raise ValueError(
                        f"Caracteres invalidos {sorted(invalid)} "
                        f"en {filepath}"
                    )

                valid_lines.append(clean_line)

            sequence = "".join(valid_lines)

        if not sequence:
            raise ValueError(
                f"No se encontro una secuencia valida en {filepath}"
            )

        sequences.append((current_header, sequence))

    for line in lines:

        if line.startswith(">"):

            if header is not None:
                save_sequence(header, content_lines)

            header = line[1:]
            content_lines = []

        else:
            content_lines.append(line)

    if header is not None:
        save_sequence(header, content_lines)

    return sequences


# 3. Coding

NUC_TO_INDEX = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}


def one_hot_encode(sequence):
    """
    OHE de dimensión L x 4.
    """
    X = np.zeros((len(sequence), 4), dtype=np.uint8)

    for i, nuc in enumerate(sequence):
        X[i, NUC_TO_INDEX[nuc]] = 1

    return X


def one_hot_padded_flat(sequence, max_length):
    """
    OHE + padding con [0,0,0,0] hasta max_length
    y aplanado a 4 * max_length caracteristicas.
    """
    X = np.zeros((max_length, 4), dtype=np.uint8)

    for i, nuc in enumerate(sequence):
        X[i, NUC_TO_INDEX[nuc]] = 1

    return X.reshape(-1)


def one_mer_encode(sequence):
    """
    Frecuencias normalizadas de 1-mers:
    [freq_A, freq_C, freq_G, freq_T]
    """
    counts = np.zeros(4, dtype=np.float64)

    for nuc in sequence:
        counts[NUC_TO_INDEX[nuc]] += 1.0

    return counts / len(sequence)


def positional_weight(i):
    """
    w(i) = 2^(-1 / (2(1+i)))
    """
    return 2 ** (-1 / (2 * (1 + i)))


def wdk_encode(sequence):
    """
    WDK de 4 componentes, normalizado por longitud real.
    Mantiene enumerate(..., start=1) para reproducir
    EncodingNucleotidsV2.py.
    """
    values = np.zeros(4, dtype=np.float64)

    for i, nuc in enumerate(sequence, start=1):
        values[NUC_TO_INDEX[nuc]] += positional_weight(i)

    return values / len(sequence)


# 4. Loading New DATASET 

df = pd.read_excel(
    matched_file,
    sheet_name="Muestras_seleccionadas"
)

required_columns = {
    "file_name",
    "class",
    "chromosome",
    "start_position",
    "end_position",
    "sequence_length"
}

missing = required_columns - set(df.columns)

if missing:
    raise ValueError(
        f"Faltan columnas en el dataset emparejado: {sorted(missing)}"
    )

if df["file_name"].str.contains(
    "control",
    case=False,
    na=False
).any():
    raise ValueError(
        "El dataset emparejado contiene controles."
    )

print("\n===== DATASET EMPAREJADO =====")
print("Total:", len(df))
print(df["class"].value_counts())


# 5. Recover original sequences

sequences = []
errors = []

for row_idx, row in df.iterrows():

    label = row["class"]
    filename = str(row["file_name"]).strip()

    if label == "benign":
        folder = path_benign
    elif label == "pathogenic":
        folder = path_pathogenic
    else:
        raise ValueError(
            f"Clase desconocida en fila {row_idx}: {label}"
        )

    filepath = os.path.join(folder, filename)

    if not os.path.isfile(filepath):
        errors.append(
            f"No existe: {filepath}"
        )
        sequences.append(None)
        continue

    try:
        fasta_records = read_fasta(filepath)

        if len(fasta_records) != 1:
            raise ValueError(
                f"{filename} contiene {len(fasta_records)} "
                "secuencias; se esperaba exactamente 1."
            )

        _, sequence = fasta_records[0]

        if len(sequence) != int(row["sequence_length"]):
            raise ValueError(
                f"Longitud distinta para {filename}: "
                f"archivo={len(sequence)}, "
                f"dataset={row['sequence_length']}"
            )

        sequences.append(sequence)

    except Exception as error:
        errors.append(
            f"{filename}: {error}"
        )
        sequences.append(None)


if errors:
    print("\n===== ERRORES AL RECUPERAR SECUENCIAS =====")
    for error in errors:
        print(error)

    raise ValueError(
        "No se pudo reconstruir exactamente el conjunto "
        "de secuencias. Corrige los errores antes de continuar."
    )

df["sequence"] = sequences

print("\nSecuencias recuperadas correctamente:", len(sequences))


# 6. Verify WDK vs saved values

wdk_columns = ["value_A", "value_C", "value_G", "value_T"]

if all(col in df.columns for col in wdk_columns):

    recalculated_wdk = np.vstack([
        wdk_encode(seq)
        for seq in df["sequence"]
    ])

    stored_wdk = df[wdk_columns].to_numpy(dtype=np.float64)

    max_abs_difference = np.max(
        np.abs(recalculated_wdk - stored_wdk)
    )

    print(
        "\nMáxima diferencia absoluta "
        "WDK recalculado vs almacenado:",
        f"{max_abs_difference:.16e}"
    )

    if not np.allclose(
        recalculated_wdk,
        stored_wdk,
        rtol=1e-10,
        atol=1e-12
    ):
        print(
            "AVISO: Los valores WDK recalculados no coinciden "
            "exactamente con los almacenados. "
            "Para la comparacion se usaran los recalculados "
            "desde las secuencias originales."
        )


# 7. Generating the three representations

max_length = int(df["sequence_length"].max())

print("\n===== DIMENSIONALIDAD =====")
print("Longitud maxima:", max_length)

print("\nGenerando WDK...")
X_wdk = np.vstack([
    wdk_encode(seq)
    for seq in df["sequence"]
])

print("Generando 1-mer...")
X_1mer = np.vstack([
    one_mer_encode(seq)
    for seq in df["sequence"]
])

print("Generando OHE + padding...")
X_ohe = np.vstack([
    one_hot_padded_flat(seq, max_length)
    for seq in df["sequence"]
])

print("\nWDK:", X_wdk.shape)
print("1-mer:", X_1mer.shape)
print("OHE + padding:", X_ohe.shape)


# 8. Label and genomic groups

y = df["class"].map({
    "benign": 0,
    "pathogenic": 1
}).to_numpy(dtype=int)

if np.isnan(y.astype(float)).any():
    raise ValueError(
        "Hay etiquetas desconocidas."
    )

window_columns = [
    "chromosome",
    "start_position",
    "end_position",
    "sequence_length"
]

groups = (
    df[window_columns]
    .astype(str)
    .agg("|".join, axis=1)
    .to_numpy()
)

print("\nVentanas genómicas distintas:", len(np.unique(groups)))


# 9. Creating the same Folds for all representations

cv = StratifiedGroupKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=SEED
)

splits = list(
    cv.split(
        X_wdk,
        y,
        groups=groups
    )
)

fold_report = []

for fold, (train_idx, test_idx) in enumerate(
    splits,
    start=1
):

    train_groups = set(groups[train_idx])
    test_groups = set(groups[test_idx])

    if train_groups & test_groups:
        raise ValueError(
            f"Hay ventanas compartidas en fold {fold}."
        )

    train_counts = np.bincount(
        y[train_idx],
        minlength=2
    )

    test_counts = np.bincount(
        y[test_idx],
        minlength=2
    )

    if (
        (train_counts == 0).any()
        or (test_counts == 0).any()
    ):
        raise ValueError(
            f"El fold {fold} no contiene ambas clases."
        )

    fold_report.append({
        "Fold": fold,
        "Train_samples": len(train_idx),
        "Test_samples": len(test_idx),
        "Train_benign": train_counts[0],
        "Train_pathogenic": train_counts[1],
        "Test_benign": test_counts[0],
        "Test_pathogenic": test_counts[1],
        "Train_windows": len(train_groups),
        "Test_windows": len(test_groups)
    })

fold_df = pd.DataFrame(fold_report)

print("\n===== PARTICIONES =====")
print(fold_df.to_string(index=False))


# 10. Computational data about representations

representations = {
    "WDK": X_wdk,
    "OHE_padding": X_ohe,
    "1-mer": X_1mer
}

representation_info = []

for name, X in representations.items():

    representation_info.append({
        "Encoding": name,
        "Samples": X.shape[0],
        "Features": X.shape[1],
        "Dtype": str(X.dtype),
        "Memory_bytes": X.nbytes,
        "Memory_MB": X.nbytes / (1024 ** 2),
        "Bytes_per_sample": X.nbytes / X.shape[0]
    })

representation_df = pd.DataFrame(
    representation_info
)

print("\n===== COSTO DE REPRESENTACION =====")
print(representation_df.to_string(index=False))


# 11. Gradient Boosting evaluation

summary_results = []
fold_results = []

for encoding_name, X in representations.items():

    print(
        f"\n========== {encoding_name} =========="
    )

    total_start = time.perf_counter()

    for fold, (train_idx, test_idx) in enumerate(
        splits,
        start=1
    ):

        X_train = X[train_idx]
        X_test = X[test_idx]

        y_train = y[train_idx]
        y_test = y[test_idx]

        model = GradientBoostingClassifier(
            random_state=SEED
        )

        train_start = time.perf_counter()

        model.fit(
            X_train,
            y_train
        )

        train_time = (
            time.perf_counter() - train_start
        )

        predict_start = time.perf_counter()

        y_pred = model.predict(
            X_test
        )

        predict_time = (
            time.perf_counter() - predict_start
        )

        accuracy = accuracy_score(
            y_test,
            y_pred
        )

        f1 = f1_score(
            y_test,
            y_pred
        )

        kappa = cohen_kappa_score(
            y_test,
            y_pred
        )

        fold_results.append({
            "Encoding": encoding_name,
            "Fold": fold,
            "Accuracy": accuracy,
            "F1": f1,
            "Kappa": kappa,
            "Train_time_s": train_time,
            "Predict_time_s": predict_time
        })

        print(
            f"Fold {fold}: "
            f"Accuracy={accuracy:.4f}, "
            f"F1={f1:.4f}, "
            f"Kappa={kappa:.4f}, "
            f"Train={train_time:.4f}s"
        )

    encoding_time = (
        time.perf_counter() - total_start
    )

    current = pd.DataFrame(
        [
            row for row in fold_results
            if row["Encoding"] == encoding_name
        ]
    )

    summary_results.append({
        "Encoding": encoding_name,
        "Features": X.shape[1],

        "Accuracy_mean": current["Accuracy"].mean(),
        "Accuracy_std": current["Accuracy"].std(ddof=1),

        "F1_mean": current["F1"].mean(),
        "F1_std": current["F1"].std(ddof=1),

        "Kappa_mean": current["Kappa"].mean(),
        "Kappa_std": current["Kappa"].std(ddof=1),

        "Train_time_mean_s":
            current["Train_time_s"].mean(),

        "Predict_time_mean_s":
            current["Predict_time_s"].mean(),

        "Total_CV_time_s": encoding_time
    })


# 12. Results

summary_df = pd.DataFrame(
    summary_results
)

fold_results_df = pd.DataFrame(
    fold_results
)

print("\n===== RESUMEN FINAL =====")
print(summary_df.to_string(index=False))


# 13. Saving

with pd.ExcelWriter(output_file) as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Resumen",
        index=False
    )

    fold_results_df.to_excel(
        writer,
        sheet_name="Resultados_por_fold",
        index=False
    )

    representation_df.to_excel(
        writer,
        sheet_name="Representaciones",
        index=False
    )

    fold_df.to_excel(
        writer,
        sheet_name="Particiones",
        index=False
    )

print("\n===== COMPARACION FINALIZADA =====")
print("Resultados guardados en:")
print(output_file)
