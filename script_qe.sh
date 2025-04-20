#!/bin/bash


# .
# ├── input/                   # input files
# │   ├── scf.in
# │   ├── nscf.in
# │   ├── bands.in
# │   ├── dos.in
# ├── run_qe.sh                # << THIS SCRIPT
# └── tmp/                     # temporary output dir



# === CONFIGURATION ===
BIN_DIR="/opt/qe741/bin"                      # Path to QE binaries (pw.x, bands.x, etc.)
PSEUDO_DIR="$HOME/pseudo"                     # Path to pseudopotentials
INPUT_DIR="./input"
TMP_DIR="./tmp"
OUTPUT_DIR="./output"
NPROC=1                                      # Number of MPI processes. Keep it =1 for GPU computation

# === ENV SETUP ===
export ESPRESSO_PSEUDO=$PSEUDO_DIR
mkdir -p "$TMP_DIR" "$OUTPUT_DIR"

# === FUNCTION TO RUN QE ===
run_qe() {
    stage=$1
    infile="${INPUT_DIR}/${stage}.in"
    outfile="${OUTPUT_DIR}/${stage}.out"
    echo ">>> Running $stage ..."
    mpirun -np $NPROC "$BIN_DIR/pw.x" -in "$infile" > "$outfile"
}

# === WORKFLOW ===
run_qe scf
run_qe nscf

# === POST-PROCESSING ===

echo ">>> Running bands.x ..."
"$BIN_DIR/bands.x" < "$INPUT_DIR/bands.in" > "$OUTPUT_DIR/bands_post.out"

"$BIN_DIR/bands.x" < "$INPUT_DIR/bands.pp.in" > "$OUTPUT_DIR/bands_post.out"

echo ">>> Running dos.x ..."
"$BIN_DIR/dos.x" < "$INPUT_DIR/dos.in" > "$OUTPUT_DIR/dos_post.out"

echo ">>> All calculations complete!"
