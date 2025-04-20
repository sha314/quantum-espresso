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
run_qe bands

# === POST-PROCESSING ===

# After running the bands calculation (bands.in with pw.x), 
# you use bands.pp.in with bands.x to post-process the data and extract eigenvalues, suitable for plotting.
echo ">>> Running bands.x ..."
mpirun -np $NPROC "$BIN_DIR/bands.x" -in "$INPUT_DIR/bands.pp.in" > "$OUTPUT_DIR/bands_post.out"

echo ">>> Running dos.x ..."
mpirun -np $NPROC "$BIN_DIR/dos.x" -in "$INPUT_DIR/dos.in" > "$OUTPUT_DIR/dos_post.out"

echo ">>> All calculations complete!"
