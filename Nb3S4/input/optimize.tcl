package require pwtk

# First doing the vc-relax computation
# then optimizing parameters


# --- user section ---------------------------------
set name "nb3s4"
set prefix   nb3s4      ;# change to your prefix
set outdir   ./tmp        ;# change to your outdir
set bandfile  bands.dat    ;# bands.x output file
set dosfile  tdos.dat     ;# dos.x output file
# --------------------------------------------------

# "load_fromPWI" command loads the pw.x input file
load_fromPWI scf.in


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

CELL_PARAMETERS angstrom "
   $cell
"

ATOMIC_POSITIONS  (crystal) "
   $pos
"

# runPW $name.scf

set total_energy [pwo_totene $name.scf.out]
puts "Total energy fron NSCF = $total_energy eV"


#------------------------
# scan the energy-cuttofs
#------------------------
# the "scanpar" command loops over a parameter; 
# it is similar to a one-parameter Tcl "foreach" command

scanpar eVar {20.0 30.0 40.0 50.0 60.0 80.0 100.0} {
   # load the new energy cuttofs
   SYSTEM "ecutwfc = $eVar, ecutrho = 4*$eVar"

   # perform the pw.x calculation (I/O files: scf.Si_e$eVar.in & scf.Si_e$eVar.out)
   runPW scf.energy$eVar
}


foreach eVar {20.0 30.0 40.0 50.0 60.0 80.0 100.0} {
   # store $eVar and total energy in the 'ecut.dat' file 

   # Get total energy from output file
   set Etot [pwo_totene scf.energy$eVar.out]

   # Append both ecut and energy to the file
   write ecut.dat "$eVar $Etot\n"

   # write ecut.dat [pwo_totene scf.energy$eVar.out]
}

# # plot the result
# plot -xl "ecutwfc (Ry)" -yl "Total energy (Ry)" ecut.dat





# ---------------------------------------------------
# LOOP OVER DIFFERENT K-POINT MESHES
# ---------------------------------------------------
# Format: {nx ny nz}
foreach kmesh {
    {3 3 6}
    {6 6 12}
    {8 8 16}
    {10 10 20}
    {12 12 24}
} {
    puts $kmesh
    lassign $kmesh nx ny nz

    # Set k-points
    K_POINTS  automatic "
      $nx $ny $nz 0 0 0
   "

    # Run SCF
    runPW scf.kmesh_${nx}x${ny}x${nz}

    # Collect total energy
    set Etot [pwo_totene scf.kmesh_${nx}x${ny}x${nz}.out]

    # Write to convergence data file
    write kpoints.dat "$nx $ny $nz $Etot\n"
}





###  ecutrho optimization
load_fromPWI scf.in
puts "Running ecutrho convergence test"
set eVar 60

scanpar mult {4 6 8 10 12} {
   # load the new energy cuttofs
   set ecutrho [expr {$eVar * $mult}]
   SYSTEM "ecutwfc = $eVar, ecutrho = $ecutrho"
   
   K_POINTS  automatic "
      12 12 24 0 0 0
   "

   # perform the pw.x calculation (I/O files: scf.Si_e$eVar.in & scf.Si_e$eVar.out)
   # runPW scf.energy_ecutrho$ecutrho
}


scanpar mult {4 6 8 10 12} {
   set ecutrho [expr {$eVar * $mult}]

   # Get total energy from output file
   set Etot [pwo_totene scf.energy_ecutrho$ecutrho.out]

   # Append both ecut and energy to the file
   write ecutrho.dat "$ecutrho $Etot\n"

   # write ecut.dat [pwo_totene scf.energy$eVar.out]
}



##########################
## degauss optimization
load_fromPWI scf.in
puts "Running degauss convergence test"

SYSTEM "ecutwfc = $eVar, ecutrho = $ecutrho"
scanpar degauss {.0001 .0005 .005 .05} {
   # load the new energy cuttofs
   SYSTEM "degauss = $degauss"
   
   K_POINTS  automatic "
      12 12 24 0 0 0
   "

   # perform the pw.x calculation (I/O files: scf.Si_e$eVar.in & scf.Si_e$eVar.out)
   # runPW scf.degauss$degauss
}


scanpar degauss {.0001 .0005 .001 .005 .01 .02 .04 .05} {
   

   # Get total energy from output file
   set Etot [pwo_totene scf.degauss$degauss.out]

   # Append both ecut and energy to the file
   write degauss.dat "$degauss $Etot\n"

}


