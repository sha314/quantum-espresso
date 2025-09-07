package require pwtk

# First doing the vc-relax computation
# then using optimized parameters,
# Doing SCF, NSCF, Bands and DOS with one pwtk

# --- user section ---------------------------------
set name "si"
set prefix   silicon      ;# change to your prefix
set outdir   ./tmp        ;# change to your outdir
set bandfile  bands.dat    ;# bands.x output file
set dosfile  tdos.dat     ;# dos.x output file
# --------------------------------------------------

# "load_fromPWI" command loads the pw.x input file
load_fromPWI initial-scf.in


CONTROL { 
   calculation = 'vc-relax'
   etot_conv_thr = 1e-5
   forc_conv_thr = 1e-4
   }

IONS {}

CELL {}

# runPW $name.relax

set output_relax $name.relax.out


# last (relaxed) CELL_PARAMETERS block
set cell [exec sh -c "awk '
    /Begin final coordinates/  {final=1}
    final && /CELL_PARAMETERS/ {flag=1; next}
    flag && /^ *$/             {exit}
    flag                       {print}' $output_relax"]

# last (relaxed) ATOMIC_POSITIONS (crystal) block
set pos [exec sh -c "awk '
    /Begin final coordinates/     {final=1}
    final && /ATOMIC_POSITIONS/   {flag=1; next}
    flag && (/^ *$/ || /End final/) {exit}
    flag                          {print}' $output_relax"]


puts $cell
puts $pos


# 2.  SCF CALCULATION

SYSTEM {
   ibrav = 0
}

CELL_PARAMETERS alat "
   $cell
"

ATOMIC_POSITIONS alat "
   $pos
"

runPW $name.scf

set total_energy [pwo_totene $name.scf.out]
puts "Total energy fron NSCF = $total_energy eV"


# 3. NSCF CALCULATION (uniform grid for DOS)

CONTROL { calculation = 'nscf' }

K_POINTS automatic { 16 16 16 0 0 0 }

SYSTEM {
   nosym   = .true.
   occupations = 'tetrahedra'
}

runPW $name.nscf

# 4. BANDS CALCULATION (high-symmetry path automatically generated)

CONTROL { calculation = 'bands' }

SYSTEM {
   occupations = 'fixed'
}

K_POINTS crystal_b {
   10                        !   Number   of   high-symmetry   points 
   0.0    0.0    0.0    40   !   Γ 
   0.5    0.0    0.5    30   !   X 
   0.5    0.25   0.75   25   !   W 
   0.375  0.375  0.75   40   !   K 
   0.0    0.0    0.0    40   !   Γ 
   0.5    0.5    0.5    40   !   L 
   0.625  0.25   0.625  40   !   U 
   0.5    0.25   0.75   40   !   W 
   0.5    0.5    0.5    40   !   L 
   0.375  0.375  0.75   40   !   K
}  

runPW  $name.bands

# 5. Running bands post PROCESSING

BANDS "
    filband  = '$bandfile'
    lsym     = .true.     
"

runBANDS $name.bands.pp

# 6. DOS POST-PROCESSING

DOS "
    fildos   = '$dosfile'
    Emin     = -30.0
    Emax     = 30.0
    DeltaE   = 0.005
    degauss  = 0.005
"

runDOS $name.dos.pp

set efermi [::pwtk::pwo::efermi  $name.nscf.out]
puts "Fermi energy = $efermi eV"


