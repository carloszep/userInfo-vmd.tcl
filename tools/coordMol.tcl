#|-coordMol006.tcl :|{condText} ================================================
#|  -tcl script for VMD-v.1.9.1 to be sourced in the TkConsole .
#|  -version :-0.0.6 ;
#|  -version information :
#|    -the resulting side length of the cubic solvent box is reported ;
#|  -date :-2016-06-09.Thu ;
#|  -proc coordMol {IDmol {selTxt "all"} {minPad "none"}} :
#|    -report for an atom selection in a specified molecule :
#|      -cartesian coordinates for the geometric center .
#|      -cartesian coordinates for the center of mass .
#|      -dimensions of the cube containing the atoms .
#|      -dimensions (lengths in cartesian coordinates) of a cubic solvent box
#|       _ with a specified minimum solvent pad around the solute atoms ;;
puts "center and box-coordinates dimensions of a molecule v.0.0.6"
puts "usage: coordMol <id> \[<selTxt>\] \[<minPad>\]"

proc coordMol {IDmol {selTxt "all"} {minPad "none"}} {
  puts "\n== Center and box size for molecule: [molinfo $IDmol get name] =="
  puts "  molecule Id: $IDmol"
  puts "  atom selection text: $selTxt"
  puts "  min solvent pad: $minPad"
  set selMol [atomselect $IDmol $selTxt]
  puts "-- Coordinates of the center: --"
  foreach c [measure center $selMol] {
    puts [format "%2.3f" $c]
    }
  puts "-- Center of mass: --"
  foreach c [measure center $selMol weight mass] {
    puts [format "%2.3f" $c]
    }
  set cuboCoord [measure minmax $selMol]
  set iniCoord [lindex $cuboCoord 0]
  set finCoord [lindex $cuboCoord 1]
  puts "-- Size of the coordinates box: --"
  set max 0.0
  set dimBox {}
  for {set c 0} {$c <= 2} {incr c} {
    set dim [expr {[lindex $finCoord $c] - [lindex $iniCoord $c]}]
    lappend dimBox $dim
    if {$dim >= $max} {set max $dim}
    puts [format "%2.3f" $dim]
    }
  $selMol delete
  if {$minPad != "none"} {
    set d [expr {$max + $minPad*2.0}]
    puts "-- Solvent pads for a cubic box of length [format "%2.3f" $d] (Angstroms): --"
    for {set c 0} {$c <= 2} {incr c} {
      puts [format "%2.3f" \
                [expr {(($max + ($minPad*2.0)) - [lindex $dimBox $c])/2.0}]]
      }
    }
  }

