
import os
import ast
import numpy as np
import pandas as pd

# 1. CONFIGURATION

base_path = r"C:\Users\Yuyu\Documents\Maestria\brca1-variant-embedding-main\brca1-variant-embedding-main\Sequences"

path_benign = os.path.join(base_path, "SeqBRCA1_Benigns")
path_pathogenic = os.path.join(base_path, "SeqBRCA1_Pathogenics")

output_file = os.path.join(
    base_path,
    "BRCA1_dataset_encoded_v2.xlsx"
)

# 2. READING SEQUENCES

def read_fasta(filepath):
    import ast
    import re

    sequences = []

    with open(filepath, "r", encoding="utf-8-sig") as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError("El archivo esta vacio.")

    header = None
    content_lines = []

    def save_sequence(current_header, current_lines):

        if current_header is None:
            raise ValueError(
                "Se encontro una secuencia sin encabezado FASTA."
            )

        if not current_lines:
            raise ValueError(
                "El encabezado no tiene una secuencia."
            )

        content = "\n".join(current_lines).strip()

        if content.startswith("["):

            try:
                parsed = ast.literal_eval(content)

            except (ValueError, SyntaxError) as error:
                raise ValueError(
                    "No se pudo interpretar la lista de nucleótidos."
                ) from error

            if not isinstance(parsed, list):
                raise ValueError(
                    "El contenido entre corchetes no es una lista."
                )

            if not all(
                isinstance(nuc, str) and len(nuc) == 1
                for nuc in parsed
            ):
                raise ValueError(
                    "La lista contiene elementos inválidos."
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
                        f"Linea no reconocida: {line!r}. "
                        f"Caracteres invalidos: {sorted(invalid)}"
                    )

                valid_lines.append(clean_line)

            sequence = "".join(valid_lines)

        if not sequence:
            raise ValueError(
                "No se encontro una secuencia valida."
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
# 3. header coordinates

def parse_header(header):

    for part in header.split():

        if part.startswith("chromosome:"):

            chrom_data = part.split(":")

            if len(chrom_data) >= 6:

                genome_version = chrom_data[1]
                chromosome = chrom_data[2]

                start = int(chrom_data[3])
                end = int(chrom_data[4])

                strand = chrom_data[5]

                return (
                    genome_version,
                    chromosome,
                    start,
                    end,
                    strand
                )

    raise ValueError(
        "No se encontraron coordenadas genómicas "
        "válidas en el encabezado."
    )


# 4. ONE-HOT ENCODING

def one_hot_encode(sequence):

    mapping = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1]
    }

    return np.array(
        [mapping[nuc] for nuc in sequence],
        dtype=np.float64
    )

# 5. WDK

def positional_weight(i):

    return 2 ** (-1 / (2 * (1 + i)))

# 6. CODING SEQUENCE

def encrypt_sequence(one_hot_matrix):

    value_A = 0.0
    value_C = 0.0
    value_G = 0.0
    value_T = 0.0
    for i, vec in enumerate(one_hot_matrix, start=1):

        w = positional_weight(i)

        value_A += vec[0] * w
        value_C += vec[1] * w
        value_G += vec[2] * w
        value_T += vec[3] * w

    return value_A, value_C, value_G, value_T


# 7. PROCESSING FOLDER

def process_folder(folder_path, label):

    if not os.path.isdir(folder_path):
        raise FileNotFoundError(
            f"No existe la carpeta: {folder_path}"
        )

    data = []
    excluded_controls = 0
    errors = []

    for filename in sorted(os.listdir(folder_path)):

        if not filename.lower().endswith(
            (".fasta", ".fa", ".txt")
        ):
            continue

        if "control" in filename.lower():

            excluded_controls += 1
            print(f"Control excluido: {filename}")

            continue

        filepath = os.path.join(folder_path, filename)

        try:

            sequences = read_fasta(filepath)

            for header, seq in sequences:

                (
                    genome_version,
                    chromosome,
                    start,
                    end,
                    strand
                ) = parse_header(header)

                if start <= 0 or end < start:

                    raise ValueError(
                        f"Coordenadas incorrectas: "
                        f"start={start}, end={end}"
                    )

                reference_length = end - start + 1

                sequence_length = len(seq)
                one_hot = one_hot_encode(seq)

                (
                    val_A,
                    val_C,
                    val_G,
                    val_T
                ) = encrypt_sequence(one_hot)

                data.append({

                    "file_name": filename,
                    "header": header,

                    "genome_version": genome_version,
                    "chromosome": chromosome,

                    "start_position": start,
                    "end_position": end,

                    "genomic_length": reference_length,
                    "sequence_length": sequence_length,

                    "strand": strand,

                    "value_A": val_A / sequence_length,
                    "value_C": val_C / sequence_length,
                    "value_G": val_G / sequence_length,
                    "value_T": val_T / sequence_length,

                    "class": label
                })

        except Exception as error:

            errors.append({
                "file_name": filename,
                "class": label,
                "error": str(error)
            })

            print(
                f"ERROR en {filename}: {error}"
            )

    return data, excluded_controls, errors


# 8. BENIGNS PATHOGENICS PROCESSING

print("\n===== PROCESANDO BENIGNAS =====")

(
    data_benign,
    controls_benign,
    errors_benign
) = process_folder(
    path_benign,
    "benign"
)

print("\n===== PROCESANDO PATOGÉNICAS =====")

(
    data_pathogenic,
    controls_pathogenic,
    errors_pathogenic
) = process_folder(
    path_pathogenic,
    "pathogenic"
)

# 9. VALIDATION

df = pd.DataFrame(
    data_benign + data_pathogenic
)

errors = errors_benign + errors_pathogenic

print("\n===== RESUMEN =====")

print("Benignas codificadas:", len(data_benign))
print("Patogénicas codificadas:", len(data_pathogenic))

print(
    "Controles excluidos:",
    controls_benign + controls_pathogenic
)

print("Archivos con errores:", len(errors))
print("Total de registros:", len(df))

if errors:

    error_file = os.path.join(
        base_path,
        "Encoding_Errors_v2.xlsx"
    )

    pd.DataFrame(errors).to_excel(
        error_file,
        index=False
    )

    print("\nSe guardó el reporte de errores en:")
    print(error_file)

    raise ValueError(
        "Se encontraron archivos con errores. "
        "Revisa Encoding_Errors_v2.xlsx antes "
        "de generar el dataset definitivo."
    )


if df.empty:

    raise ValueError(
        "No se codificó ninguna secuencia."
    )


controls_in_dataset = df[
    df["file_name"].str.contains(
        "control",
        case=False,
        na=False
    )
]

if not controls_in_dataset.empty:

    raise ValueError(
        "Se encontraron controles en el dataset."
    )

features = [
    "value_A",
    "value_C",
    "value_G",
    "value_T",
    "start_position",
    "sequence_length"
]

if df[features].isna().any().any():

    raise ValueError(
        "El dataset contiene valores faltantes."
    )


if not np.isfinite(
    df[features].to_numpy(dtype=float)
).all():

    raise ValueError(
        "El dataset contiene valores infinitos."
    )


# 10. EXPORT NEW DATASET 

df.to_excel(
    output_file,
    index=False
)

print("\n===== CODIFICACIÓN FINALIZADA =====")

print("Controles incluidos: 0")
print("Total de muestras:", len(df))

print("\nClases:")
print(df["class"].value_counts())

print("\nArchivo guardado en:")
print(output_file)
