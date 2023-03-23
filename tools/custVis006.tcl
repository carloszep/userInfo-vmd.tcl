#|-custVis006.tcl :
#|  -tcl script library for VMD to manage custom visualizations .
#|  -this library is intended to apply pre-established visualization styles
#|   _ while VMD is running .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :-2017-06-06.Thu ;
#|  -version :-0.0.6 ;
#|  -version information :
#|    -added features :
#|      -procedure rep_newRibbonsGlass, rep_cpkEdgyShiny added
#|       _ and commonRep updated .
#|      -optional argument added to rep_cpk1 ;
#|    -pending features :
#|      -commonRep should be declared with a variable number of arguments and
#|       _ a lot of default values .
#|      -it may be convenient to declare an array to store representations .
#|      -a proc for turning on and off specific representations in molecules .
#|      -to declare the proc defVis to restore to the VMD's default
#|       _ visualization settings .
#|      -a procedure to replace a representation type of all the loaded
#|       _ molecules to an specific representation (as CPK for example) .
#|      -show an introduction an quick explanation about usage ;
#|    -finished version ;
#|  -global variables and parameters :
#|    -custVis_version :-script version in format #.#.# ;
#|    -custVis_list :-list of custom visualization names ;
#|    -custVis_loSt :-global log output stream .
#|      -by default, all output log and messages will be sent to the console ;;
#|  -procedures and functions :
#|    -whiteVis :-apply custom global visualization options ;
#|    -blueBeta :-global visualization option to highlight beta-sheets
#|     _ in blue ;
#|    -commonRep :-creates a predefined representation to a list of molecules ;
#|    -delLastRep :-deletes the last representation of a list of molecules ;
#|    -rep_cpk1 :-set up the properties of a representation in CPK ;
#|    -rep_newCartoon1 :-set up the properties of a rep in newCartoon ;
#|    -rep_newRibbons1 :-set up the properties of a rep in newRibbons ;
#|    -rep_licorice1 :-set up the properties of a representation in Licorice ;
#|    -rep_vdw1 :-set up the properties of a representation in VDW ;;
#|  -code details :

# constants, global and default parameters
global custVis_version custVis_list custVis_loSt
set custVis_version 0.0.6
set custVis_list [list whiteVis blueBeta]   ;# initialization
set custVis_loSt stdout

# +++ showing introduction +++ pending....

# +++ procedures to change global visualization +++

#|-proc whiteVis {} :
#|  -sets white display, gray carbons and silver hydrogens and correct the
#|   _ color of the labels ;
proc whiteVis {} {
  global custVis_loSt
  puts $custVis_loSt "\nsetting custom global visualization to: whiteVis..."
  color Display Background white
  color Name C gray
  color Name H white
  color Restype Nonpolar silver
  color Labels Atoms cyan3
  color Labels Bonds purple
  color Labels Angles violet2
  color Labels Dihedrals orange2
  }

#|-proc blueBeta {} :
#|  -change structure motifs, i.e. beta-sheet: blue,
#|   _ turn: green, random coil: silver, etc. .
#|  -starts from whiteVis ;
proc blueBeta {} {
  global custVis_loSt
  whiteVis
  puts $custVis_loSt "\nsetting custom visualization: blueBeta..."
  color change rgb 24 0.140000 0.180000 0.680000
  color change rgb 20 0.180000 0.690000 0.180000
  color change rgb 22 0.000000 0.760000 1.000000
  color Structure Extended_Beta blue3
  color Structure Bridge_Beta cyan3
  color Structure Turn green3
  color Structure Coil silver
  color Structure 3_10_Helix orange2
  }

#|-proc commonRep {molList repOp {selTxt "all"} {color "Name"}} :
#|  -creates a new representation for a list of loaded molecules .
#|  -arguments :
#|    -molList :-tcl list with VMD molIDs aready loaded ;
#|    -repOp :
#|      -acceptable values :
#|        -cpk1, newCartoon1, licorice1, newRibbons1, newRibbonsBGlass ;
#|    -selTxt :-atom selection considered for the new representation ;
#|    -color :-coloring method string, ex. 'ColorID 0' ;;;
#|  -notes :
#|    -the top flag most be used in the repScript .
#|    -the original top molecule is reset at the end of the proc ;;
proc commonRep {molList repOp {selTxt "all"} {color "Name"}} {
  switch $repOp {
    "cpk1" {rep_cpk1 $selTxt}
    "cpkEdgyShiny" {rep_cpkEdgyShiny $selTxt $color}
    "newCartoon1" {rep_newCartoon1 $selTxt}
    "newRibbons1" {rep_newRibbons1 $selTxt}
    "newRibbonsBGlass" {rep_newRibbonsBGlass $selTxt $color}
    "licorice1" {rep_licorice1 $selTxt $color}
    "vdw1" {rep_vdw1 $selTxt $color}
    default {return "Option not found!"}
    }
  foreach id $molList {
    mol addrep $id
    }
  }   ;# commonRep

#|-proc delLastRep {molList} :
#|  -deletes the last representation in each of a set of molecules in molList ;
proc delLastRep {molList} {
  foreach id $molList {
    mol delrep [expr [molinfo $id get numreps] - 1] $id
    }
  }   ;# delLastRep

#|-rep_cpk1 {{selTxt "all"}} :
#|  -set up the properties of a representation to be created .
#|  -it may include changes in color, rep type, material, etc. .
#|  -the rep is actually not created here, the command mol addrep id should
#|   _ be used after .
#|  -provisionally the selection is all ;
proc rep_cpk1 {{selTxt "all"} {color Name}} {
  mol color $color
  mol representation CPK 0.600000 0.300000 13.000000 13.000000
  mol selection $selTxt
  mol material Diffuse
  }   ;# rep_cpk1

proc rep_cpkEdgyShiny {{selTxt "all"} {color Name}} {
  mol color $color
  mol representation CPK 0.600000 0.300000 13.000000 13.000000
  mol selection $selTxt
  mol material EdgyShiny
  }   ;# rep_cpkEdgyShiny

#|-rep_newCartoon1 {{selTxt "all"}} :
#|  -set up the properties of a representation to be created .
#|  -the rep is actually not created here, the command mol addrep id should
#|   _ be used after ;
proc rep_newCartoon1 {{selTxt "all"}} {
  mol color Structure
  mol representation NewCartoon 0.300000 13.000000 4.100000 0
  mol selection $selTxt
  mol material Diffuse
  }   ;# rep_newCartoon1

#|-proc rep_newRibbons1 {{selTxt "all"}} :
#|  -set up the properties of a representation to be created .
#|  -the rep is actually not created here, the command mol addrep id should
#|   _ be used after ;
proc rep_newRibbons1 {{selTxt "all"}} {
  mol color Structure
  mol representation NewRibbons 0.300000 13.000000 3.000000 0
  mol selection $selTxt
  mol material Diffuse
  }   ;# rep_newRibbons1

#|-proc rep_newRibbonsBGlass {{selTxt "all"} {color "Molecule"}} :
#|  -set up newRibbons representation transparent to BlownGlass material ;
proc rep_newRibbonsBGlass {{selTxt "all"} {color "Name"}} {
  mol color $color
  mol representation NewRibbons 0.300000 13.000000 3.000000 0
  mol selection $selTxt
  mol material BlownGlass
  }   ;# rep_newRibbonsBGlass

#|-proc rep_licorice1 {{selTxt "all"} {color "Name"}} :
#|  -set up the properties of a representation to be created .
#|  -the rep is actually not created here, the command mol addrep id should
#|   _ be used after .
#|  -arguments :
#|    -selTxt :-atomselection to be considered in the rep ;
#|    -color :-coloring method ;;;
proc rep_licorice1 {{selTxt "all"} {color "Name"}} {
  mol color $color
  mol representation Licorice 0.100000 13.000000 13.000000
  mol selection $selTxt
  mol material Diffuse
  }   ;# rep_licorice1

#|-proc rep_vdw1 {{selTxt "all"} {color "Name"}} :
#|  -set up the properties of a representation to be created .
#|  -the rep is actually not created here, the command mol addrep id should
#|   _ be used after .
#|  -arguments :
#|    -selTxt :-atomselection to be considered in the rep ;
#|    -color :-coloring method ;;;
proc rep_vdw1 {{selTxt "all"} {color "Name"}} {
  mol color $color
  mol representation VDW 0.500000 13.000000
  mol selection $selTxt
  mol material Diffuse
  }   ;# rep_vdw1

# +++ applying default custom visualization +++ ;# for the moment none

#|    - ;;
#axes location Off

