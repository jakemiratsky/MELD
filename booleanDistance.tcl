# Load the PDB file in vmd
resetpsf
# Select only C-alpha atoms
set sel [atomselect top "protein and name CA"]

# Create a list of residue names and initialize the distance matrix
set resnames [$sel get resname]
puts $resnames
set num_res [llength $resnames]
puts $num_res
set distmat {}
for {set i 0} {$i < $num_res} {incr i} {
  set row {}
  for {set j 0} {$j < $num_res} {incr j} {
    lappend row 0.0
  }
  lappend distmat $row
}

# Measure distances between all pairs of C-alpha atoms
set num_atoms [$sel num]
puts $num_atoms
set residues [$sel get residue]
puts $residues
set indices [$sel get index]
puts $indices
for {set i 0} {$i < $num_atoms} {incr i 1} {
  for {set j [expr {$i+1}]} {$j < $num_atoms} {incr j} {
    set res_i [lindex $residues $i]
    set res_j [lindex $residues $j]
    set atom_i [lindex $indices $i]
    set atom_j [lindex $indices $j]
    set dist [measure bond "$atom_i $atom_j"]
    puts $dist
    if {$dist < 10.0} {
        set dist 1
    }
    if {$dist > 10.0} {
        set dist 0
    }
    lset distmat $res_i $res_j $dist
    lset distmat $res_j $res_i $dist
    }
}

# Output the distance matrix to a file
set outfile "distance.txt"
set fptr [open $outfile "w"]
puts $fptr "Residue Distances (Angstroms): "
puts $fptr [join $resnames "\t"]
foreach row $distmat {
    set rowstr [join $row "\t"]
    puts $fptr "$rowstr"
}
close $fptr
puts "Distance matrix saved to $outfile"

# Cleanup
$sel delete
mol delete top

