#|-anMD.tcl :| {condtext}
#|  -Structural analysis of MD trajectories .
#|  -Scripts library in tcl for VMD v.1.9.2 .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :-2020-08-10.Mon ;
#|  -version :-0.6.0 ;
#|  -version information :
#|    -changes in this version :
#|      -the call to fibVsecStructCA.tcl is updated for the GitHub version :
#|        -the script is read from the anMD/tools/ directory ;;
#|    -previous changes :
#|      -updating to use fibVsecStructCA.tcl by secStructTraj and
#|       _ secStructAv .
#|      -an error processing the ref fragId was corrected in rmsdTrajDat ;
#|    -finished version ;
#|    -to do :
#|      -when loading a traj with trajLoad, the global loadStep has to be
#|       _ considered .
#|      -it is necessary to optimize secStructTraj and secStructAv .
#|      -several procs have to be reincorporated and updated from anMD055.tcl .
#|      -to update several proc to use the trajFragSpec procedure .
#|      -to include in setSelId a keyword to explicitly set default values .
#|      -some analysis procedures need to detect common molId between selIds
#|       _ to avoid unnecessary multiple trajectory loads .
#|      -to implement a progress bar or indicator in all the analysis .
#|      -to implement optional arguments to the procedures needing it .
#|      -to apply the trajSelId namespace to all analysis procedures .
#|      -to write a general procedure to create data labels and output file
#|       _ names based on the selId/trajInfo specifications .
#|      -to include proc putsLog to manage log information output .
#|      -to correct the number of arguments in all procs calling proc
#|       _ getSelInfo .
#|      -update all procedures to work with the new tragFrag scheme .
#|      -the equilibreation fragments could be arranged to display negative
#|       _ values to be plotted together with the prod simulation .
#|      -structural analysis :
#|        -procs anMD_ncIntTraj, anMD_ncIntSels .
#|        -calculation and managing of volumetric data (with volmap) .
#|        -to quantify HBs between the protein and the solvent .
#|        -to calculate force-field energies ;
#|      -plotting procedures :
#|        -to incorporate capability to create gnuplot graphs .
#|        -to generate xmgrace .agr plots from previous .dat files .
#|        -to write a proc to smooth data .
#|        -to include more than one graph, i.e. a matrix of graphs .
#|        -to include optional input in agrXY for specifying colors of series
#|         _ in plots ;
#|      -general library issues :
#|        -to document in detail the selInfo and trajInfo arrays .
#|        -to consider the use of the molInfo array to make anMol library
#|         _ compatible ;;;
#|  -List of global variables and parameter :
#|    -anMD_version :-script version in format #.#.# ;
#|    -anMD_logFileName :-default file name for log output ;
#|    -selInfo... :| (this has to be documented in detail)
#|      -array storing information related to atom selections ;
#|    -trajInfo... :- ;| (this has to be documented in detail)
#|      -array storing information about the trajectories loaded in VMD ;
#|    -anMD_graphicsOn :-stablishes whether graphics are active ;
#|  -List of procedures and functions :
#|    -structural analysis procedures :
#|      -alignTraj :-align each frame of a structure in a trajectory ;
#|      -rmsdTrajDat :-calculates the RMSD for a trajectory in a selection ;
#|      -rmsdTraj :-RMSD for a set of selections ;
#|      -distTrajDat :-calculates distances along a traj ;
#|      -angleTrajDat :-calculates geometric angles along a trajectory ;
#1      -dihedTrajDat :-calculate dihedral angles along a trajectory ;
#|      -geomTraj :-plots structural geometric measurements for selIds ;
#|      -attribTrajDat :-monitors properties of a selection over the time ;
#|      -attribTraj :-monitor a property for a set of selections ;
#|      -anMD_gofrTraj :-calculates the g(r) function for a trajectory ;
#|      -gofrTraj :-Calculation of g(r) function for trajectories ;
#|      -radGyrTrajDat :-radious of gyration for a trajectory ;
#|      -radGyrTraj :-radious of gyration for a set of trajectories ;
#|      -anMD_hbCount :-Number of hydrogen bonds for a single frame ;
#|      -anMD_hbCountTraj :-Number of hydrogen bonds for a trajectory ;
#|      -anMD_hbCountSels :-Plots the number of hydrogen bonds vs time for a
#|       _ set of trajectories ;
#|      -anMD_ncIntCount :-Number of non-covalent interactions ;
#|      -sasaTrajDat :-Solvent-accesible surface area for a trajectory ;
#|      -sasaTraj :-Solvent-accesible surface area for a set of trajs ;
#|      -secStructTraj :-secundary structure incidence for a trajectory ;
#|      -secStructAv :-average secundary structure for categories ;;
#|    -trajInfo- and selInfo-independent procedures :
#|      -nI_rmsdFrame :-calculates the RMSD for a single frame ;
#|      -nI_gofr :-runs the VMD function gofr (g(r)) and returns the result ;
#|      -nI_hbCount :-count the number of hydrogen bonds ;
#|      -nI_sasa :-calculates the solvent-accesible surface area ;
#|      -frame2time :-time in ns of a specific frame in a trajectory ;
#|      -time2frame :-frame number correspondig a specific time in a traj ;;
#|    -utility procedures :
#|      -agrXY :-creates an XY-data agr file for ploting in xmgrace ;
#|      -strNames :-return a comma-separated string made by several tokens ;
#|      -infoAtom :-return all the information associated to an atom ;
#|      -aaCode_1to3 :-return a 3-letter amino acid code from 1-letter ;
#|      -aaCode_3to1 :-return a 1-letter amino acid code from 3-letter ;
#|      -aaCode_capToCam :-return a 3-letter amino acid code in cammel style ;;
#|    -testing procedures :
#|      - ;
#|  -code details :
#

# constants, global and default parameters
global anMD_version anMD_logFileName anMD_maxTrajSizeGB anMD_currTrajSize
global anMD_graphicsOn ind
set anMD_version 0.6.0
set anMD_logFileName "log_anMD-v.$anMD_version.txt"
set anMD_graphicsOn 1



# +++ showing introduction +++

puts "\nanMD-v$anMD_version: Structural analysis of MD trajectories for VMD."
puts "Carlos Z. Gómez Castro, 2010-2017.\n"



# +++ structural analysis procedures +++


#|-proc alignTraj {selId args} :
#|  -alignment (to a reference structure/frame) of a selection of atoms along
#|   _ a trajectory .
#|  -alignment implies the rotation and traslation of the coordinates of the
#|   _ whole molecule so that the rmsd of the selection is minimized .
#|  -traslational and rotational degrees of freedom are deleted from a traj .
#|  -arguments :
#|    -selId :
#|      -selection identifier to access the selInfo array .
#|      -defines the atom-selection text, the molId and range of frames
#|       _ to be aligned :
#|        -selInfo keys considered :
#|          -"selTxt", "molId", "fragId", "first", and "last" .
#|          -"frame" if argument selIdRef has value "selId" ;;
#|    -args (variable arguments) :
#|      -"selIdRef", "refSelId", "ref" :
#|        -selection identifier of the reference structure for alignment .
#|        -selInfo keys considered :-"selTxt", "molId", "frame" ;
#|        -default frame specification is "now" or the current frame .
#|        -acceptable values :
#|          -any key idetifier (selId) for the selInfo array .
#|          -"selId" :
#|            -the same as the selId argument is used ;;
#|        -default value :-"selId" ;;
#|      -"weight", "weights", "wt" :
#|        -perform the alignment considering a specified weighting
#|         _ squeme of the atom's coordinates .
#|        -acceptable values ... :
#|          -any atom-selection property (data field) or a list of weights .
#|          -"none" :-no weights are considered ;
#|          -"mass" :-mass of each atom is used for weighting ;
#|          -"charge" :-atom charges are used for weighting the coordinates ;;
#|        -default value :-"mass" ;;
#|      -"order", "indOrderList", "reorder" :
#|        -list of 0-based indices indicating how to reorder the atoms in
#|         _ selection 2 .
#|        -order optional flag of the VMD measure fit command .
#|        -default value :-{} ;;
#|      -"exclude", "except", "exFragId" - "excl" :
#|        -fragId exclude list for the trajFragSpec proc .
#|        -default value :-{} ;;
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;;;
#|  -notes :
#|    -the "order" option of the VMD's measure fit proc is available but not
#|     _ tested .
#|    -the traj fragments specified by the selId's must be already loaded for
#|     _ this proc to work (no traj is loaded) .
#|    -the step key of the selInfo array is ignored .
#|    -the proc has to be updated to manage the new 'selLbl' and 'des'
#|     _ variables from proc 'getSelInfo' ;;
proc alignTraj {selId args} {
# invoke global variables and check default values for arguments
  global trajInfo selInfo
  namespace import trajSelId::getDatX trajSelId::getFrame
# default values for arguments
  set loSt stdout; set procName [lindex [info level 0] 0]
  set selIdRef "selId"; set weight "mass"; set exclude {}; set fragIdOver ""
  set order {}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "selIdRef" - "refSelId" - "ref" {set selIdRef $val}
          "weight" - "weights" - "wt" {set weight $val}
          "order" - "indOrderList" - "reorder" {set order {}}
          "exclude" - "except" - "exFragId" - "excl" {set exclude $val}
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {puts $loSt "$procName: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "$procName: Odd number of variable arguments! args: $args"
      return ""
      }
  if {$selIdRef == "selId"} {set selIdRef $selId}
# extracting and reporting input information
  puts $loSt "\nalignTraj: Aligning structure in a trajectory. selId: $selId"
  puts $loSt "  Reference structure selId: $selIdRef"
  if {[llength $weight] > 1} {
    puts $loSt "  Using user-provided weights"
  } else {
    puts $loSt "  Weighting scheme: $weight"}
  showSelInfo $selId $loSt
  if {$selId != $selIdRef} {showSelInfo $selIdRef $loSt}
# information for the reference selection
  getSelInfo $selIdRef refTxt refId frag frmR frst lst stp upd ldstp selLbl des
  if {![getTrajProp "loaded" $refId $frag]} {
    puts $loSt "alignTraj: Error: Reference structure not loaded."
    return
    }
  set frag [trajFragSpec $frag $refId exclude $exclude loSt $loSt]
  set frag [simFrag $frmR $refId $frag] 
  puts $loSt "\n  Reference atom-selection text: $refTxt"
  puts $loSt "  Reference molId: $refId"
  puts $loSt "  Reference time/frame: $frmR"
# information for molecule to be aligned
  getSelInfo $selId selTxt id frag frm frst lst stp upd ldstp selLbl des
  if {![getTrajProp "loaded" $id $frag]} {
    puts $loSt "alignTraj: Error: Structure not loaded."
    return
    }
  puts $loSt "  Molecule atom-selection text: $selTxt"
  puts $loSt "  Molecule to be aligned: $id"
  puts $loSt "  Trajectory frames to be aligned: $frst - $lst"
# atom selections for alignment  
  set refSel [atomselect $refId $refTxt frame $frmR]
  set molSel [atomselect $id $selTxt]
  set totSel [atomselect $id "all"]
#  puts -nonewline $loSt "\naligning frame:"
  while {[trajSelId::iterate $selId keep $frag excl $exclude log $loSt]} {
#    puts -nonewline $loSt " [getFrame]"
    $molSel frame [getFrame]
    $totSel frame [getFrame]
    if {$order == {}} {
      $totSel move [measure fit $molSel $refSel weight $weight]
    } else {
      $totSel move [measure fit $molSel $refSel weight $weight order $order]
      }
    }
  puts $loSt ""
  $refSel delete
  $molSel delete
  $totSel delete
  puts $loSt "\nalignTraj: Done."
  }   ;# alignTraj

#|-proc rmsdTrajDat {selId {refSelId "selId"} {alignSelTxt ""} {lbl "selId"}
#|           _ {weight "mass"} {seqTime 1} {excl {}} {ord {}} {loSt stdout}} :
#|  -calculates the RMSD for an atom selection along a trajectory .
#|  -returns a data list with the format {l_timeDat l_rmsdDat setLbl} used by
#|   _ the agrXY proc .
#|  -arguments :
#|    -selId :
#|      -selection Id (selInfo array) of the structure which the RMSD will
#|       _ be calculated for .
#|      -selInfo array keys considered :
#|        -molId, selTxt, fragId, first, last, step, and loadStep .
#|        -if the trajectory is already loaded the loadStep is ignored ;;
#|    -"refSelId", "selIdRef", "ref" :
#|      -selId to access the selInfo array for reference structure .
#|      -selInfo array keys to be considered :
#|        -molId, selTxt, fragId, and frame ;
#|      -acceptable values :
#|        -"selId" :-the same molecules as that indicated by the selId
#|         _ argument is used as reference ;;
#|      -default value :-"selId" ;;
#|    -alignSelTxt :
#|      -VMD string for an atom selection to perform the struct alignment only .
#|      -the atom selection must be compatible for both selId and refSelId .
#|      -the rmsd is calculated for a different atom selection, thus it is
#|       _ not minimized (different from 'the best fit') .
#|      -only overrides the "selTxt" from the selId/refSelId for alignment,
#|       _ the rest of the selInfo keys are taken the same .
#|      -default value :-"" ;
#|      -notes :
#|        -to consider to accept a selId instead of a selection text ;;
#|    -lbl :
#|      -data set label to identify the data in plots .
#|      -short string (preferently with no blank spaces) .
#|      -included at index position 2 of the returned output data list .
#|      -acceptable values :
#|        -"selId" :-the specified selId name is used as label ;
#|        -"selInfo" :-the label in the selInfo array is used ;
#|        -"refSelId", "selIdRef", "ref" :
#|          -the specified refSelId is used as label ;
#|        -"fragId", "frag" :-the fragId specified in the selInfo is used ;
#|        -"name" :-the name stored in trajInfo(id,name) is used ;
#|        -"full", "auto" :-a label is generated with format:
#|         _ idName(selId)_refIdName(refSelId) ;;
#|      -default value :-"selId" ;;
#|    -weight :
#|      -weighting scheme considered to calculate the RMSD .
#|      -the coordinates of each atom are weighted (scaled) by a property .
#|      -acceptable values... :
#|        -any atom-selection property (data field) or a list of weights .
#|        -"none" :-no weights are considered ;
#|        -"mass" :-mass of each atom is used for weighting ;
#|        -"charge" :-atom charges are used for weighting the coordinates ;
#|        -other values are possible... ;
#|      -default value :-"mass" ;;
#|    -seqTime :
#|      -indicate whether the trajectory x data should be taken as time (ns)
#|       _ or as frame number .
#|      -"0" value can be used to ignore time data of a simulation .
#|      -default value :-1 ;;
#|    -excl :-fragId exclude list for the trajFragSpec proc .
#|      -default value :-{} ;;
#|    -ord :
#|      -list of 0-based indices indicating how to reorder the atoms in
#|       _ selection 2 .
#|      -"order" optional flag of the VMD measure fit command .
#|      -default value :-{} ;;
#|    -loSt :-output stream for log messages .-default value :-stdout ;;;
#|  -notes :
#|    -if the fragId traj is not loaded it will be loaded for the
#|     _ calculation, and then deleted .
#|    -if the trajFrags considered in the selId are all loaded in a consecutive
#|     _ order but have a non-consecutive time course, the frame number instead
#|     _ of the time will be considered in the X-axis .
#|    -all the trajectory frames considered in selId will be aligned with the
#|     _ the reference frame .
#|    -the proc has to be updated to manage the new variables 'selLbl' and 'des'
#|     _ from proc 'getSelInfo' ;;
proc rmsdTrajDat {selId {refSelId "selId"} {alignSelTxt ""} {lbl "selId"} \
               {weight "mass"} {seqTime 1} {excl {}} {ord {}} {loSt stdout}} {
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
  set prevMolId [molinfo top]
  set prevFrame [molinfo $prevMolId get frame]
  puts $loSt "\nrmsdTrajDat: Calculating RMSD data for a trajectory..."
# extracting information for the reference structure
  if {$refSelId == "selId"} {set refSelId $selId}
  getSelInfo $refSelId refTxt refId fragR tmR frst lst stp upd ldstp selLbl des
  if {$tmR == "now"} {set tmR [getTrajProp "iniTime" $refId $fragR excl $excl]}
# limit the ref frag specification to the single frag containing the ref time
  set fragR [trajFragSpec $fragR $refId loSt $loSt]
  set fragR [simFrag $tmR $refId $fragR]
  if {![getTrajProp "loaded" $refId $fragR]} {set loadR 1} else {set loadR 0}
  if {$loadR} {trajLoad $refId $fragR loadStep $ldstp loSt $loSt}
  set frmR [simFrame $tmR $refId $fragR]
  puts $loSt " Ref. structure:"
  puts $loSt "  selId: $refSelId  id: $refId  name: $trajInfo($refId,name)"
  puts $loSt "  trajFrag: $fragR  ref time: $tmR  frame $frmR"
  puts $loSt "  atom selection: $refTxt"
  if {$alignSelTxt != ""} {
    puts $loSt "Atom selection for structural alignment: $alignSelTxt"}
# extracting information for the trajectory
  getSelInfo $selId selTxt id frag frm frst lst stp upd ldstp selLbl des
# report default values
  puts $loSt " Trajectory:"
  puts $loSt "  selId: $selId  id: $id  name: $trajInfo($id,name)"
  puts $loSt "  trajFrag: $frag  fragNames: \
                  [trajFragSpec $frag $id excl $excl loSt $loSt]"
  puts $loSt "  time Interval: $frst - $lst  step: $stp  load step: $ldstp"
  puts $loSt "  atom selection: $selTxt"
  puts $loSt " Weighting scheme: $weight"
# setting up the output label for files and legends
  switch $lbl {
    "auto" - "full" {
      if {$selId == $refSelId} {
        set lbl "$trajInfo($id,name)($selId)_ref.t$tmR"
      } else {
        set lbl \
          "$trajInfo($id,name)($selId)_$trajInfo($refId,name)($refSelId)"}
      }
    "selId" {set lbl $selId}
    "selInfo" {set lbl $selLbl}
    "refSelId" - "selIdRef" - "ref" {set lbl $refSelId}
    "fragId" - "frag" {set lbl [join $frag ","]}
    "name" {set lbl $trajInfo($id,name)}
    }
  puts $loSt " Data set label: $lbl"
# setting up the atom selections
  set refSel [atomselect $refId $refTxt frame $frmR]   ;# ref atom selection
  set molSel [atomselect $id $selTxt]   ;# traj atom selection
  set totSel [atomselect $id all]   ;# total traj atom selection
  if {$alignSelTxt != ""} {
    set refSelAlign [atomselect $refId $alignSelTxt frame $frmR]
    set molSelAlign [atomselect $id $alignSelTxt]
  } else {
    set refSelAlign [atomselect $refId $refTxt frame $frmR]
    set molSelAlign [atomselect $id $selTxt]
    }
# calculating rmsd
  set l_xDat {}
  set l_yDat {}
  trajSelId::setSeqTime $seqTime
  while {[trajSelId::iterate $selId keep $fragR excl $excl log $loSt]} {
    $molSel frame [getFrame]
    $totSel frame [getFrame]
    $molSelAlign frame [getFrame]
    if {$ord == {}} {
      $totSel move [measure fit $molSelAlign $refSelAlign weight $weight]
    } else {
      $totSel move [measure fit $molSelAlign $refSelAlign weight $weight \
                                                              order $ord]
      }
    lappend l_xDat [getDatX]
    lappend l_yDat [measure rmsd $molSel $refSel weight $weight]
    }
  $refSel delete
  $molSel delete
  $totSel delete
  $refSelAlign delete
  $molSelAlign delete
  if {$loadR} {trajDelete $refId $fragR loSt $loSt}
  puts $loSt "\nrmsdTrajDat: done."
# return to the original molId and frame
  mol top $prevMolId
  if {$prevFrame != -1} {animate goto $prevFrame}
  return [list $l_xDat $l_yDat $lbl]
  }   ;# rmsdTrajDat

#|-proc rmsdTraj {selIds args} :
#|  -Generates a RMSD plot for atom selection(s) along a trajectory .
#|  -there is a common reference for all the trajectories specified .
#|  -uses rmsdTrajDat to calculate the RMSD, see that proc for the options .
#|  -arguments :
#|    -selIds :-list of atom selection Ids for the rmsdTrajDat proc ;
#|    -args (variable argument names) :
#|      -argument names for rmsd calculation :
#|        -"refSelId", "selIdRef", "ref", "reference" :
#|          -equivalent to the "refSelId" argument for the rmsdTrajDat proc .
#|          -acceptable values :
#|            -a single selId that will be considered as a reference for all
#|             _ selIds specfied :
#|              -all atom selections must be compatible ;
#|            -"selId" :
#|              -same as "selId" value in the rmsdTrajDat proc .
#|              -will apply to each of the specified selIds to calculate ;
#|            -a list of selIds that will be considered as reference for each
#|             _ of the specified selIds to calculate, respectively :
#|              -this requires the list to be the same length as the selIds
#|               _ list (and optionally to the label list) ;;
#|          -default value :-"selId" ;;
#|        -"alignSelTxt", "align", "alignSel" :
#|          -"alignSelTxt" argument for proc "rmsdTrajDat" .
#|          -a list of values corresponding to each of the selIds specified .
#|          -note :-for a single sel text value formed by several words (i.e.
#|           _ "protein and name CA") it is necessary to enclose in curly
#|           _ brackets ({"protein and name CA"}) to avoid the text to be
#|           _ splitted as a list ;
#|          -default value :-"" ;;
#|        -"weight", "weights", "wt" :
#|          -weighting scheme for rmsd calculation and structure alignment .
#|          -acceptable values :
#|            -the same as the "weight" argument for the rmsdTrajDat proc .
#|            -a list of values corresponding to each of the selIds specified ;
#|          -default value :-"mass" ;;
#|        -"order", "ord", "indOrderList", "reorder" :
#|          -"order" flag of the VMD's measure fit command .
#|          -acceptable values :
#|            -the same as the "ord" argument of the rmsdTrajDat proc .
#|            -list of lists of reorder indices corresponding to each of the
#|             _ selIds specified ;
#|          -default value :-{} ;;
#|        -"exclude", "except", "exFragId", "excl" :
#|          -"exclude" frgId list for the trajFragSpec proc .
#|          -default value :-{} ;;
#|        -"seqTime", "nsScale", "nsTime" :
#|          -see "seqTime" argument in rmsdTrajDat .
#|          -default value :-"1" ;;;
#|      -argument names for plotting set up :
#|        -"label", "labels", "lbl", "legend", "dataLabel", "datLabel" :
#|          -equivalent to the "lbl" argument for the rmsdTrajDat proc .
#|          -acceptable values :
#|            -the same as the "lbl" argument of the rmsdTrajDat proc .
#|            -a list of values corresponding to each of the selIds specified ;
#|          -default value :-"auto" ;;
#|        -"prefix", "pref", "outPrefix", "agrName" :
#|          -prefix for .agr and .dat output file names .
#|          -the rest of the file names are built from other parameters .
#|          -acceptable values :
#|            -the same as the "pref" argument of the agrXY proc ;
#|          -default value :-"rmsd" ;;;
#|        -"title", "titl", "plotTitle", "graphTitle" :
#|          -title for the xmgrace plot .
#|          -default value :-"RMSD for selection Id(s): <sels> " ;;
#|        -"subtitle", "stitl", "plotSubtitle", "graphSubtitle" :
#|          -subtitle for the xmgrace plot .
#|          -default value :-"reference: <refSelId>" ;;
#|        -"xLabel", "labelX", "xlbl" :-label for the x axis (abscissa) .
#|          -default value :-"Time (ns)" ;;
#|        -"yLabel", "labelY", "ylbl" :-label for the y axis (ordinate) .
#|          -default value :-"RMSD (Angstroms)" ;;;
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;;;
#|  -notes :
#|    -xLabel optional variable unable to be defined by the user .
#|    -seqTime option is not working .
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'des' from proc 'getSelInfo' ;;
proc rmsdTraj {selIds args} {
  global trajInfo
# default values for arguments
  set loSt stdout; set procName [lindex [info level 0] 0]
  set refSelId "selId"; set weight "mass"; set excl {}; set seqTime 1
  set order {}; set label "auto"; set pref "rmsd"; set align ""
  set titl "RMSD for selection Id(s): [join $selIds ","]"
  set stitl "reference: $refSelId"; set xlbl "Time (ns)"
  set ylbl "RMSD (Angstroms)"
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "selIdRef" - "refSelId" - "ref" - "reference" {
            set refSelId $val
            set stitl "reference: [join $refSelId ","]"}
          "alignSelTxt" - "align" - "alignSel" {set align $val}
          "weight" - "weights" - "wt" {set weight $val}
          "order" - "ord" - "indOrderList" - "reorder" {set order $val}
          "exclude" - "except" - "exFragId" - "excl" {set exclude $val}
          "seqTime" - "nsScale" - "nsTime" {set seqTime $val}
          "label" - "labels" - "lbl" - "legend" - "dataLabel" - "datLabel" {
            set label $val}
          "prefix" - "pref" - "outPrefix" - "agrName" {set pref $val}
          "title" - "titl" - "plotTitle" - "graphTitle" {set titl $val}
          "subtitle" - "stitl" - "plotSubtitle" - "graphSubtitle" {
            set stitl $val}
          "xLabel" - "labelX" - "xlbl" {set xlbl $val}
          "yLabel" - "labelY" - "ylbl" {set ylbl $val}
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {puts $loSt "$procName: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "$procName: Odd number of variable arguments! args: $args"
      return ""
      }
# show info
  puts $loSt "\n$procName: RMSD calculation for a set of atom selections..."
  puts $loSt "  List of atom selection Ids: $selIds"
  puts $loSt "  Reference selection Id: $refSelId"
# analyze input arguments
  set nSel [llength $selIds]
  foreach var {refSelId align weight order label} {
    set val [set $var]
    set nArg [llength $val]
    if {$nArg == 1} {
      for {set i 1} {$i < $nSel} {incr i} {lappend $var [lindex $val 0]}
    } elseif {(($nArg > 1) && ($nArg < $nSel)) || ($nArg > $nSel)} {
      puts $loSt "$procName: Incorrect number of values for argument $var"
      return ""
      }
    }
  set l_plots {}
  set simNames {}
  foreach selId $selIds {
    getSelInfo $selId selTxt id frag frm frst lst stp upd ldstp selLbl des
    if {$seqTime} {
      set seqTime [getTrajProp "seqTime" $id $frag excl $excl loSt $loSt]}
    if {$trajInfo($id,name) != $simNames} {lappend simNames $trajInfo($id,name)}
    }
  for {set i 0} {$i < $nSel} {incr i} {
    set sId [lindex $selIds $i]
    set ref [lindex $refSelId $i]
    set al [lindex $align $i]
    set wt [lindex $weight $i]
    set ord [lindex $order $i]
    set lbl [lindex $label $i]
    lappend l_plots \
      [rmsdTrajDat $sId $ref $al $lbl $wt $seqTime $excl $ord $loSt]
    }
  if {$seqTime} {set xlbl "Time (ns)"} else {set xlbl "Frame number"}
  set sels [join $selIds ","]
  puts $loSt "\n$procName: Creating RMSD plot for selections $sels..."
  agrXY $l_plots "${pref}_[join $simNames ","]" \
        title $titl subtitle $stitl \
        xLabel $xlbl yLabel $ylbl loSt $loSt
  puts $loSt "$procName: Done."
  }   ;# rmsdTraj


#|-proc distTrajDat {selId1 selId2 {weight "mass"} {lbl "selId"} {seqTime 1}
#|                                              _ {fixSel2 0} {loSt stdout}} :
#|  -calculates the distance between two atom selections along a trajectory .
#|  -returns a list of lists with the data for plotting with proc 'agrXY' .
#|  -arguments :
#|    -selId1, selId2 :
#|      -selId (keys) to access the selInfo array storing atom select info .
#|      -selInfo keys from selId1 used to define trajectory parameters :
#|        -'molId', 'fragId', 'first', 'last', 'step', 'updateSel',
#|         _ and 'loadStep' ;
#|      -selInfo keys taken independently from each selId :
#|        -'selTxt', 'label', and desc .
#|        -'molId' when 'fixSel2' equals '1' ;
#|      -selInfo keys taken from selId2 :-'frame' when 'fixSel2' equals '1' ;;
#|    -weight :
#|      -weighting scheme to calculate the geometric center of each selection
#|       _ of atoms .
#|      -used when a selId comprises more than one atom in the selection .
#|      -acceptable values :
#|        -"none" :-no weights are considered ;
#|        -a list of lists with values (weights) per atom in the selTxt for
#|         _ each selId :-NOTE: option not implemented yet ;
#|        -an atom-selection property (data field), i.e. "mass", "charge",... ;
#|      -default value :-"mass" ;;
#|    -lbl :
#|      -data label (legend) to be included in the returned list .
#|      -this label may be included in plots .
#|      -acceptable values :
#|        -any short string (preferently with no blank spaces) to label data .
#|        -"selId" :-the specified selId will be used ;
#|        -"selInfo" :-the label stored in the selInfo array is used ;
#|        -"fragId" :-the fragId specified in the selInfo array is used ;
#|        -"name" :-the name specified in the trajInfo array is used ;
#|        -"full", "auto" :-a label is generated with format:
#|         _ 'geomTerm_selId1-selId2' ;;
#|      -default value :-"auto" ;;
#|    -seqTime :
#|      -flag for the type of time data used by namespace 'trajSelId' .
#|      -acceptable values :
#|        -"0" :-frame number is used in the 'x' axis ;
#|        -"1" :-the time (in nanoseconds) is used in the 'x' axis ;;
#|      -default value :-"0" ;;
#|    -fixSel2 :-flag indicating that selId2 should consider a single
#|     _ structure/frame .
#|      -acceptable values :
#|        -'0' :
#|          -atom selection 2 (selId2) is expected to change with time .
#|          -will consider the same molId and traj parameters as selId1 ;
#|        -'1' :
#|          -atom selection 2 (selId2) is fixed, a single structure is used .
#|          -allows different molIds for each selId ;;
#|      -default value :-'0' ;;
#|    -loSt :-output stream for log messages .
#|      -default value :-stdout ;;;
#|  -notes :
#|    -if a selection comprises more than one atom, the weighted
#|     _ geometric center will be considered to make the calculation .
#|    -'exclude' lists independent from fragId specs are not implemented .
#|    -for the moment the calculation is done only for a single molecule unless 
#|     _ the selId2 was fixed (single strucutre):
#|      -to consider more than one molecule a sort of 'time synchronization'
#|       _ would have to be defined ;
#|    -consider to change the scheme to specify either sel1 or 2 as fixxed ;;
proc distTrajDat {selId1 selId2 {weight "mass"} {lbl "selId"} {seqTime 1} \
                                                {fixSel2 0} {loSt stdout}} {
# namespace and global variables invocation
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
# storing initial 'top' id and frame number
  set prevMolId [molinfo top]
  set prevFrame [molinfo $prevMolId get frame]
# showing input information
  set procName [lindex [info level 0] 0]
  puts $loSt "\n$procName: Distance between atom selections for trajs."
  puts $loSt "  selIds: $selId1, $selId2"
  puts $loSt "  opts: weight: $weight, label: $lbl, seqTime: $seqTime"
# extracting selInfo from selIds
  getSelInfo $selId2 selTxt2 id2 frag2 frm2 fst lst stp upd ldStp2 selLbl2 desc
  getSelInfo $selId1 selTxt1 id frag frm fst lst stp upd ldStp selLbl1 desc
  puts $loSt "  molId: $id   frag: $frag"
  puts $loSt "  selTxt1 (label: $selLbl1): $selTxt1"
  puts $loSt "  selTxt2 (label: $selLbl2): $selTxt2"
  puts $loSt "  time interval: $fst - $lst   step: $stp   loadStep: $ldStp"
# setting up data label
  switch $lbl {
    "selId" {set lbl "$selId1-$selId2"}
    "selInfo" {set lbl "$selLbl1-$selLbl2"}
    "fragId" {set lbl "$frag-$frag2"}
    "name" {set lbl $trajInfo($id,name)}
    "full" - "auto" {set lbl "dist_$selId1-$selId2"}
    }
# computing distance from atom selections
  if {$fixSel2} {
# calculate the fixed center of atom selection 2
    puts $loSt \
      "  steady selection (selId2): molId $id2  frag: $frag2  time/frame: $frm2"
    set frag2 [simFrag $frm2 $id2 [trajFragSpec $frag2 $id2 loSt $loSt]]
    if {[getTrajProp "loaded" $id2 $frag2 loSt $loSt]} {
      set loadFrag2 0
    } else {
      set loadFrag2 1
      trajLoad $id2 $frag2 loadStep $ldStp2 loSt $loSt
      }
    set sel2 [atomselect $id2 $selTxt2 frame [simFrame $frm2 $frag2]]
    set center2 [measure center $sel2 weight $weight]
    } else {set sel2 [atomselect $id2 $selTxt2]}
  set sel1 [atomselect $id $selTxt1]
  set l_xDat {}
  set l_yDat {}
  trajSelId::setSeqTime $seqTime
  while {[trajSelId::iterate $selId1 loSt $loSt]} {
    $sel1 frame [getFrame]
    set center1 [measure center $sel1 weight $weight]
    if {!$fixSel2} {
      $sel2 frame [getFrame]
      set center2 [measure center $sel2 weight $weight]
      }
    lassign $center1 x1 y1 z1
    lassign $center2 x2 y2 z2
    lappend l_xDat [getDatX]
    lappend l_yDat [expr {sqrt(($x2-$x1)**2 + ($y2-$y1)**2 + ($z2-$z1)**2)}]
    }
  $sel1 delete
  $sel2 delete
  if {$fixSel2 && $loadFrag2} {trajDelete $id2 $frag2 loSt $loSt}
  puts $loSt "\n$procName: done."
# return to the original molId and frame
  mol top $prevMolId
  if {$prevFrame != -1} {animate goto $prevFrame}
  return [list $l_xDat $l_yDat $lbl]
  }   ;# distTrajDat


#|-proc angleTrajDat {selId1 selId2 selId3 {weight "mass"} {lbl "selId"}
#|                _ {seqTime 1} {fixSel23 0} {units "deg"} {loSt stdout}} :
#|  -calculation of the angle between three atom selections along a traj .
#|  -returns a lis of lists with data neccesary to plot the angle vs time .
#|  -intended to be called by proc 'geomTraj' .
#|  -arguments :
#|    -selId1, selId2, selId3 :
#|      -array keys (atom selection identifiers) to access the selInfo array .
#|      -the center of each selection is considered to calculate the angle
#|       _ for each frame in the trajectory .
#|      -selId2 is considered as the middle point to claculate the angle .
#|      -'selInfo' array keys to be considered independently from each selId :
#|        -'selTxt', 'fragId', 'frame', 'loadStep', 'label', 'desc' ;
#|      -'selInfo' array keys to be considered only from 'selId1' :
#|        -if 'fixSel23' equals 0 :-'molId', 'fragId' ;
#|        -'first', 'last', 'step' ;
#|      -'selInfo' array keys to be considered from 'selId2' and 'selId3'
#|       (if 'fixSel23' equals 1) :
#|        -'molId', 'fragId', 'frame' ;;
#|    -weight :
#|      -weighting scheme to calculate the center of each atom selection .
#|      -used when a selId comprises more than one atom in the selection .
#|      -acceptable values :
#|        -"none" :-no weights are considered ;
#|        -a list of lists with values (weights) per atom in the selTxt for
#|         _ each selId .
#|        -an atom-selection property (data field), i.e. "mass", "charge",... ;
#|      -default value :-"mass" ;;
#|    -lbl :-data label (legend) to be included in the returned list .
#|      -this label may be included in plots .
#|      -acceptable values :
#|        -any short string (preferently with no blank spaces) to label data .
#|        -"selId" :-the specified selId will be used ;
#|        -"selInfo" :-the label stored in the selInfo array is used ;
#|        -"name" :-the name specified in the trajInfo array is used ;
#|        -"full", "auto" :-a label is generated with format:
#|         _ 'geomTerm_selId1-selId2-selId3' ;;
#|      -default value :-"auto" ;;
#|    -seqTime :
#|      -flag for the type of time data used by namespace 'trajSelId' .
#|      -acceptable values :
#|        -"0" :-frame number is used in the 'x' axis ;
#|        -"1" :-the time (in nanoseconds) is used in the 'x' axis ;;
#|      -default value :-"1" ;;
#|    -fixSel23 :-flag indicating that selId2 and selId3 should consider a
#|     _ single structure/frame .
#|      -acceptable values :
#|        -'0' :
#|          -selId2 and selId3 are expected to change with time .
#|          -will consider the same molId and traj parameters as selId1 ;
#|        -'1' :
#|          -selId2 and selId3 are fixed, a single structure is used .
#|          -allows different molIds for each selId ;;
#|      -default value :-'0' ;;
#|    -units :-units for the angle data returned .
#|      -acceptable values :
#|        -"rad", "radians" :-angles in radians ;
#|        -"deg", "degrees" :-angles in degrees ;;
#|      -default value :-"deg" ;;
#|    -loSt :-output stream for log messages .
#|      -default value :-stdout ;;;
#|  - ;;
proc angleTrajDat {selId1 selId2 selId3 {weight "mass"} {lbl "selId"} \
                 {seqTime 1} {fixSel23 0} {units "degrees"} {loSt stdout}} {
# namespace and global variables invocation
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
# storing initial 'top' id and frame number
  set prevMolId [molinfo top]
  set prevFrame [molinfo $prevMolId get frame]
# showing input information
  set procName [lindex [info level 0] 0]
  puts $loSt "\n$procName: Angle between atom selections for trajs."
  puts $loSt "  selIds: $selId1, $selId2, $selId3"
  puts $loSt \
    "  opts: weight: $weight, label: $lbl, seqTime: $seqTime, units: $units"
# extracting selInfo from selIds
  getSelInfo $selId3 selTxt3 id3 frag3 frm3 fst lst stp upd ldStp3 selLbl3 desc
  getSelInfo $selId2 selTxt2 id2 frag2 frm2 fst lst stp upd ldStp2 selLbl2 desc
  getSelInfo $selId1 selTxt1 id1 frag1 frm1 fst lst stp upd ldStp1 selLbl1 desc
  puts $loSt "  molId: $id1   frag: $frag1"
  puts $loSt "  selTxt1 (label: $selLbl1): $selTxt1"
  puts $loSt "  selTxt2 (label: $selLbl2): $selTxt2"
  puts $loSt "  selTxt3 (label: $selLbl3): $selTxt3"
  puts $loSt "  time interval: $fst - $lst   step: $stp   loadStep: $ldStp1"
# setting up data label
  switch $lbl {
    "selId" {set lbl "$selId1-$selId2-$selId3"}
    "selInfo" {set lbl "$selLbl1-$selLbl2-$selLbl3"}
    "name" {set lbl $trajInfo($id,name)}
    "full" - "auto" {set lbl "angle_$selId1-$selId2-$selId3"}
    }
# computing angle from atom selections
  if {$fixSel23} {
# calculate the fixed center of atom selection 2
    puts $loSt \
      "  steady selection (selId2): molId $id2  frag: $frag2  time/frame: $frm2"
    puts $loSt \
      "  steady selection (selId3): molId $id3  frag: $frag3  time/frame: $frm3"
    set frag2 [simFrag $frm2 $id2 [trajFragSpec $frag2 $id2 loSt $loSt]]
    if {[getTrajProp "loaded" $id2 $frag2 loSt $loSt]} {
      set loadFrag2 0
    } else {
      set loadFrag2 1
      trajLoad $id2 $frag2 loadStep $ldStp2 loSt $loSt
      }
    set sel2 [atomselect $id2 $selTxt2 frame [simFrame $frm2 $frag2]]
    set center2 [measure center $sel2 weight $weight]
# calculate the fixed center of atom selection 3
    set frag3 [simFrag $frm3 $id3 [trajFragSpec $frag3 $id3 loSt $loSt]]
    if {[getTrajProp "loaded" $id3 $frag3 loSt $loSt]} {
      set loadFrag3 0
    } else {
      set loadFrag3 1
      trajLoad $id3 $frag3 loadStep $ldStp3 loSt $loSt
      }
    set sel3 [atomselect $id3 $selTxt3 frame [simFrame $frm3 $frag3]]
    set center3 [measure center $sel3 weight $weight]
  } else {
    set sel2 [atomselect $id2 $selTxt2]
    set sel3 [atomselect $id3 $selTxt3]
    }
  set sel1 [atomselect $id1 $selTxt1]
  set l_xDat {}
  set l_yDat {}
  trajSelId::setSeqTime $seqTime
  while {[trajSelId::iterate $selId1 loSt $loSt]} {
    $sel1 frame [getFrame]
    set center1 [measure center $sel1 weight $weight]
    if {!$fixSel23} {
      $sel2 frame [getFrame]
      set center2 [measure center $sel2 weight $weight]
      $sel3 frame [getFrame]
      set center3 [measure center $sel3 weight $weight]
      }
    set vec21 [vecsub $center1 $center2]
    set vec23 [vecsub $center3 $center2]
    lappend l_xDat [getDatX]
    set angleRad [expr \
      {acos([vecdot $vec21 $vec23]/([veclength $vec21]*[veclength $vec23]))}]
    if {($units == "rad") || ($units == "radians")} {
      lappend l_yDat $angleRad
    } else {
      lappend l_yDat [expr {$angleRad*180/3.141592}]
      }
    }
  $sel1 delete
  $sel2 delete
  $sel3 delete
  if {$fixSel23 && $loadFrag2} {trajDelete $id2 $frag2 loSt $loSt}
  if {$fixSel23 && $loadFrag3} {trajDelete $id3 $frag3 loSt $loSt}
  puts $loSt "\n$procName: done."
# return to the original molId and frame
  mol top $prevMolId
  if {$prevFrame != -1} {animate goto $prevFrame}
  return [list $l_xDat $l_yDat $lbl]
  }   ;# angleTrajDat


#|-proc dihedTrajDat {selId1 selId2 selId3 selId4 {weight "mass"} {lbl "selId"}
#|             _ {seqTime 1} {fixSel234 0} {units "degrees"} {loSt stdout}} { :
#|  -calculates the dihedral angle between four atom selections along a traj .
#|  -returns a list of lists wtih the calculated data to plot with agrXY .
#|  -arguments :
#|    -selId1, selId2, selId3, selId4 :
#|      -atom selection identifers (array keys) to access the selInfo array .
#|      -define the four geometric points necessary to calculate a dihedral .
#|      -selId1 :
#|        -is the leading atom selection where most of the trajectory info
#|         _ is taken from .
#|        -array keys (fields) to be considered :
#|          -'selTxt', 'molId', 'fragId', 'first', 'last', 'step',
#|           _ 'loadStep', 'label', and 'desc' ;;
#|      -selId2, selId3, selId4 :
#|        -these atom selection may be chosen to be fixed ('fixSel234' arg) .
#|        -may refer to molIds different from that specified by selId1 if
#|         _ a fixed (single strucutre) is considered through arg 'fixSel234' .
#|        -array keys (fields) to be considered :
#|          -if arg 'fixSel234' is "0" :
#|            -'selTxt', 'label' ;
#|          -if arg 'fixSel234' is "1" :
#|            -'selTxt', 'molId', 'fragId', 'frame', 'loadStep',
#|             _ 'label' ;;;;
#|    -weight :
#|      -weighting scheme to calculate the center of each atom selection .
#|      -used when a selId comprises more than one atom in the selection .
#|      -acceptable values :
#|        -"none" :-no weights are considered ;
#|        -a list of lists with values (weights) per atom in the selTxt for
#|         _ each selId .
#|        -an atom-selection property (data field), i.e. "mass", "charge",... ;
#|      -default value :-"mass" ;;
#|    -lbl :
#|      -data label (legend) to be included in the returned data list .
#|      -this label might be included in plots .
#|      -acceptable values :
#|        -any short string (preferently with no blank spaces) to label data .
#|        -"selId" :-the specified selId will be used ;
#|        -"selInfo" :-the label stored in the selInfo array is used ;
#|        -"name" :-the name specified in the trajInfo array is used ;
#|        -"full", "auto" :-a label is generated with format:
#|         _ 'geomTerm_selId1-selId2-selId3-selId4' ;;
#|      -default value :-"auto" ;;
#|    -seqTime :
#|      -flag for the type of time data used by namespace 'trajSelId' .
#|      -acceptable values :
#|        -"0" :-frame number is used in the 'x' axis ;
#|        -"1" :-the time (in nanoseconds) is used in the 'x' axis ;;
#|      -default value :-"1" ;;
#|    -fixSel234 :
#|      -flag indicating that selId2, selId3, and selId4 should consider a
#|     _ single structure/frame .
#|      -acceptable values :
#|        -"0" :
#|          -selId2, selId3, and selId4 are expected to change with time .
#|          -will consider the same molId and traj parameters as selId1 ;
#|        -"1" :
#|          -selId2, selId3, and selId4 are fixed, a single structure is used .
#|          -allows different molIds for each selId ;;
#|      -default value :-"0" ;;
#|    -units :
#|      -units used for the angle data returned .
#|      -acceptable values :
#|        -"rad", "radians" :-angles in radians ;
#|        -"deg", "degrees" :-angles in degrees ;;
#|      -default value :-"degrees" ;;
#|    -loSt :-output stream for log messages .
#|      -default value :-stdout ;;;
#|  -references :
#|    -Rahul (https://math.stackexchange.com/users/856/rahul),
#|     _ How do I calculate a dihedral angle given Cartesian coordinates?,
#|     _ URL (version: 2011-06-23): https://math.stackexchange.com/q/47084 ;;
proc dihedTrajDat {selId1 selId2 selId3 selId4 {weight "mass"} {lbl "selId"} \
                   {seqTime 1} {fixSel234 0} {units "degrees"} {loSt stdout}} {
# namespace and global variables invocation
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
# storing initial 'top' id and frame number
  set prevMolId [molinfo top]
  set prevFrame [molinfo $prevMolId get frame]
# showing input information
  set procName [lindex [info level 0] 0]
  puts $loSt "\n$procName: Dihedral angle between atom selections for trajs."
  puts $loSt "  selIds: $selId1, $selId2, $selId3, $selId4"
  puts $loSt \
    "  opts: weight: $weight, label: $lbl, seqTime: $seqTime, units: $units"
# extracting selInfo from selIds
  getSelInfo $selId4 selTxt4 id4 frag4 frm4 fst lst stp upd ldStp4 selLbl4 desc
  getSelInfo $selId3 selTxt3 id3 frag3 frm3 fst lst stp upd ldStp3 selLbl3 desc
  getSelInfo $selId2 selTxt2 id2 frag2 frm2 fst lst stp upd ldStp2 selLbl2 desc
  getSelInfo $selId1 selTxt1 id1 frag1 frm1 fst lst stp upd ldStp1 selLbl1 desc
  puts $loSt "  molId: $id1   frag: $frag1"
  puts $loSt "  selTxt1 (label: $selLbl1): $selTxt1"
  puts $loSt "  selTxt2 (label: $selLbl2): $selTxt2"
  puts $loSt "  selTxt3 (label: $selLbl3): $selTxt3"
  puts $loSt "  selTxt3 (label: $selLbl4): $selTxt4"
  puts $loSt "  time interval: $fst - $lst   step: $stp   loadStep: $ldStp1"
# setting up data label
  switch $lbl {
    "selId" {set lbl "$selId1-$selId2-$selId3-$selId4"}
    "selInfo" {set lbl "$selLbl1-$selLbl2-$selLbl3-$selId4"}
    "name" {set lbl $trajInfo($id,name)}
    "full" - "auto" {set lbl "angle_$selId1-$selId2-$selId3-$selId4"}
    }
# computing angle from atom selections
  if {$fixSel234} {
    puts $loSt \
      "  steady selection (selId2): molId $id2  frag: $frag2  time/frame: $frm2"
    puts $loSt \
      "  steady selection (selId3): molId $id3  frag: $frag3  time/frame: $frm3"
    puts $loSt \
      "  steady selection (selId4): molId $id4  frag: $frag4  time/frame: $frm4"
# calculate the fixed center of atom selection 2
    set frag2 [simFrag $frm2 $id2 [trajFragSpec $frag2 $id2 loSt $loSt]]
    if {[getTrajProp "loaded" $id2 $frag2 loSt $loSt]} {
      set loadFrag2 0
    } else {
      set loadFrag2 1
      trajLoad $id2 $frag2 loadStep $ldStp2 loSt $loSt
      }
    set sel2 [atomselect $id2 $selTxt2 frame [simFrame $frm2 $frag2]]
    set center2 [measure center $sel2 weight $weight]
# calculate the fixed center of atom selection 3
    set frag3 [simFrag $frm3 $id3 [trajFragSpec $frag3 $id3 loSt $loSt]]
    if {[getTrajProp "loaded" $id3 $frag3 loSt $loSt]} {
      set loadFrag3 0
    } else {
      set loadFrag3 1
      trajLoad $id3 $frag3 loadStep $ldStp3 loSt $loSt
      }
    set sel3 [atomselect $id3 $selTxt3 frame [simFrame $frm3 $frag3]]
    set center3 [measure center $sel3 weight $weight]
# calculate the fixed center of atom selection 4
    set frag4 [simFrag $frm4 $id4 [trajFragSpec $frag4 $id4 loSt $loSt]]
    if {[getTrajProp "loaded" $id4 $frag4 loSt $loSt]} {
      set loadFrag4 0
    } else {
      set loadFrag4 1
      trajLoad $id4 $frag4 loadStep $ldStp4 loSt $loSt
      }
    set sel4 [atomselect $id4 $selTxt4 frame [simFrame $frm4 $frag4]]
    set center4 [measure center $sel4 weight $weight]
  } else {
    set sel2 [atomselect $id2 $selTxt2]
    set sel3 [atomselect $id3 $selTxt3]
    set sel4 [atomselect $id4 $selTxt4]
    }
  set sel1 [atomselect $id1 $selTxt1]
# calculating data
  set l_xDat {}
  set l_yDat {}
  trajSelId::setSeqTime $seqTime
  while {[trajSelId::iterate $selId1 loSt $loSt]} {
    $sel1 frame [getFrame]
    set center1 [measure center $sel1 weight $weight]
    if {!$fixSel234} {
      $sel2 frame [getFrame]
      set center2 [measure center $sel2 weight $weight]
      $sel3 frame [getFrame]
      set center3 [measure center $sel3 weight $weight]
      $sel4 frame [getFrame]
      set center4 [measure center $sel4 weight $weight]
      }
# calculating dihedral between centers 1, 2, 3, 4
#    set vec12 [vecsub $center2 $center1]
#    set vec23 [vecsub $center3 $center2]
#    set vec34 [vecsub $center4 $center3]
#    set norm13 [veccross $vec12 $vec23]
#    set norm24 [veccross $vec23 $vec34]
#    set angleRad [expr {acos([vecdot $norm13 $norm24]/([veclength $norm13]*[veclength $norm24]))}]
# calculating dihedral between centers 1, 2, 3, 4
    set b1 [vecnorm [vecsub $center2 $center1]]
    set b2 [vecnorm [vecsub $center3 $center2]]
    set b3 [vecnorm [vecsub $center4 $center3]]
    set n1 [veccross $b1 $b2]
    set n2 [veccross $b2 $b3]
    set m1 [veccross $n1 $b2]
    set x [vecdot $n1 $n2]
    set y [vecdot $m1 $n2]
    set angleRad [expr {-atan2($y,$x)}]
# NOTE: the minus sign was added to coincide in sign with VMD plots
    if {($units == "rad") || ($units == "radians")} {
      lappend l_yDat $angleRad
    } else {
      lappend l_yDat [expr {$angleRad*180/3.141592}]
      }
    lappend l_xDat [getDatX]
    }
  $sel1 delete
  $sel2 delete
  $sel3 delete
  $sel4 delete
  if {$fixSel234 && $loadFrag2} {trajDelete $id2 $frag2 loSt $loSt}
  if {$fixSel234 && $loadFrag3} {trajDelete $id3 $frag3 loSt $loSt}
  if {$fixSel234 && $loadFrag4} {trajDelete $id4 $frag4 loSt $loSt}
  puts $loSt "\n$procName: done."
# return to the original molId and frame
  mol top $prevMolId
  if {$prevFrame != -1} {animate goto $prevFrame}
  return [list $l_xDat $l_yDat $lbl]
  }   ;# dihedTrajDat


#|-proc geomTraj {geomPar args} :
#|  -calculates geometric parameters for atom selections along a trajectory .
#|  -creates an xmgrace .agr plot and .dat files of the calculated data .
#|  -arguments :
#|    -geomPar :-type of geometric parameter to calculate .
#|      -acceptable values :
#|        -"center", "COM", "com", "weightedCenter" :
#|          -to be implemented... .
#|          -will write a .dat file with the following data on each column :
#|            -time (ns), x, y, and z cartesian coordinates .
#|          -note :-a 3d plot with a time-traj scale may be generated ;;;
#|        -"bond", "distance", "dist", "diffCOM", "comDiff", "COMDiff" :
#|          -the distance between two atom selections ;
#|        -"angle", "comAngle", "angleCOM" :
#|          -the (bond) angle between three atom selections ;
#|        -"dihedral", "dihed", "comDihed", "dihedCOM" :
#|          -the dihedral angle between four atom selections ;
#|        -"improper"*... ;;
#|    -args (variable arguments) :
#|      -"selId<N>", "sel<N>"; <N> may be 1 to 4 :
#|        -selection Ids necessary to calculate different geometric
#|         _ parameters .
#|        -each geomParam requires a specific number of atom selections
#|         _ (selIds) from 1 to 4 .
#|        -see procs "distTrajDat", ..., to specific details about each selId .
#|        -acceptable values :
#|          -a single or a list of selId keys to access the selInfo array :
#|            -list of selIds specified :
#|              -the number of values must be the same for each variable ;
#|            -single selId specified :
#|              -if other variables have more than one selId specified, the
#|               _ selId will be reproduced to match the number of values ;;
#|          -"none" :-indicates that the variable is not in use ;;
#|        -default value :-"none" ;;
#|      -"weight", "weigths", "wt" :
#|        -weighting scheme for geometric center of an atom selection .
#|        -used when a selId comprises more than one atom in the selection .
#|        -acceptable values :
#|          -"none" :-no weights are considered ;
#|          -a list of lists with values (weights) per atom in the selTxt for
#|           _ each selId .
#|          -an atom-selection property (data field), ("mass", "charge",...) ;
#|        -default value :-"mass" ;;
#|      -"loSt", "channelId", "log" :
#|        -output stream (channelId) for log messages -
#|        -acceptable values :
#|          -stdout :-default standard channel (output to the console) .
#|          -an identifier of an open channel (see tcl command 'open') ;
#|        -default value :-stdout ;;
#|      -"fixSel", "fixSelIds" :
#|        -flag indicating that some atom selections are fixed (steady
#|         _ structure or frame) according to pre-specifyed conventions .
#|        -acceptable values :
#|          -"0" :-no selId is considered fixed .
#|            -in this case all selId's must refer to the same molecule ;
#|          -"1" :-one or more atom selections are fixed ;
#|        -default value :-"0" ;;
#|      -"seqTime", "nsScale", "nsTime" :
#|        -see 'seqTime' argument in proc 'disTrajDat' .
#|        -default value :-1 ;;
#|      -"units", "angleUnits" :
#|        -units used for the angles measured .
#|        -acceptable values :
#|          -"deg", "degrees", "rad", "radians" ;
#|        -default value :-"degrees" ;;;
#|    -args (variable arguments) affecting the generated plots :
#|      -"label", "labels", "lbl", "legend", "dataLabel", "datLabel" :
#|        -equivalent to lbl argument for procs: 'distTrajDat',  .
#|        -acceptable values :
#|          -see 'lbl' argument for the just mentioned procedures .
#|          -a single value or a list of values for multiple data sets ;
#|        -default value :-"auto" ;;
#|      -"prefix", "pref" "fileName", "agrName", "agrPrefix" :
#|        -name prefix for data and .agr file created by proc 'agrXY' .
#|        -acceptable values :
#|          -string suitable to name a file .
#|          -"auto", "default" :
#|            -generate an automatic name prefix in the form:
#|             _ "<geomParam>_<simNames>"  ;;
#|        -default value :-"auto" ;;;;
#|  -notes :
#|    -exclude argument is disabled outside the fragId specification .
#|    -pending to fix some small details such as the name of the .agr file ;;
proc geomTraj {geomPar args} {
# global variables and namespaces invocation
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
  set procName [lindex [info level 0] 0]
# default values for variable arguments
  foreach N {1 2 3 4} {set selId$N "none"}
  set loSt stdout; set weight "mass"; set fixSel "0"; set loSt stdout
  set pref "auto"; set seqTime 1; set lbl "auto"; set units "degrees"
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "selId1" - "sel1" {set selId1 $val}
          "selId2" - "sel2" {set selId2 $val}
          "selId3" - "sel3" {set selId3 $val}
          "selId4" - "sel4" {set selId4 $val}
          "weight" - "weights" - "wt" {set weight $val}
          "loSt" - "channelId" - "log" {set loSt $val}
          "fixSel" - "fixSelId" {set fixSel $val}
          "seqTime" - "nsScale" - "nsTime" {set seqTime $val}
          "units" - "angleUnits" {set units $val}
          "label" - "labels" - "lbl" - "legend" - "dataLabel" - "datLabel" {
            set lbl $val}
          "prefix" - "pref" - "fileName" - "agrName" - "agrPrefix" {
            set pref $val}
          default {
            puts $loSt "$procName: argument unkown: $arg"
# if the argument is unknown is redirected to other procs through 'argl' list
            lappend argl $arg
            lappend argl $val
            }
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "$procName: Odd number of variable arguments! args: $args"
      return ""
      }
  switch $geomPar {
    "bond" - "distance" - "dist" - "diffCOM" - "comDiff" {
      set nSelId 2
      set l_selId [list selId1 selId2]
      set yLabel "Distance (Angstroms)"
      }
    "angle" - "comAngle" - "angleCOM" {
      set nSelId 3
      set l_selId [list selId1 selId2 selId3]
      set yLabel "Angle ($units)"
      }
    "dihedral" - "dihed" - "comDihed" - "dihedCOM" {
      set nSelId 4
      set l_selId [list selId1 selId2 selId3 selId4]
      set yLabel "Dihedral Angle ($units)"
      }
    default {
      puts $loSt "$procName: unknown geometric parameter: $geomPar"
      return ""
      }
    }
# show info
  puts $loSt "\n$procName: Calculation of geometric parameters in a trajectory"
  puts $loSt "  procedure invocation: geomTraj $geomPar $args"
# analize input arguments
  set nSel 0
  foreach selId $l_selId {
    if {[llength [set $selId]] > $nSel} {set nSel [llength [set $selId]]}
    if {[set $selId] == "none"} {
      puts $loSt "$procName: $nSelId selIds need to be specified"
      return ""
      }
    }
# check for arguments that must contain a list of values
  foreach varn $l_selId {
    foreach var $varn {
      set val [set $var]
      set nArg [llength $val]
      if {$nArg == 1} {
        for {set i 1} {$i < $nSel} {incr i} {lappend $var [lindex $val 0]}
      } elseif {(($nArg > 1) && ($nArg < $nSel)) || ($nArg > $nSel)} {
        puts $loSt "$procName: Incorrect number of values for argument $var"
        return ""
        }
      }
    }
  foreach var [list weight lbl fixSel] {
    set val [set $var]
    set nArg [llength $val]
    if {$nArg == 1} {
      for {set i 1} {$i < $nSel} {incr i} {lappend $var [lindex $val 0]}
    } elseif {(($nArg > 1) && ($nArg < $nSel)) || ($nArg > $nSel)} {
      puts $loSt "$procName: Incorrect number of values for argument $var"
      return ""
      }
    }
# check for sequential trajectories (in simulation time)
# NOTE: only selId1 is checked for the moment
  set simNames {}
  for {set i 0} {$i < $nSel} {incr i} {
    for {set N 1} {$N <= 1} {incr N} {
      set sId [lindex [set selId$N] $i]
      getSelInfo $sId selT id frag frm frst lst stp upd ldstp sLbl des
      if {$seqTime} {set seqTime [getTrajProp "seqTime" $id $frag loSt $loSt]}
      if {[lsearch $simNames $trajInfo($id,name)] == -1} {
        lappend simNames $trajInfo($id,name)}
      }
    }
  if {$seqTime} {set xLabel "Time (ns)"} else {set xLabel "Frame number"}
# creating data sets
  set plots {}
  for {set i 0} {$i < $nSel} {incr i} {
    for {set N 1} {$N <= $nSelId} {incr N} {
      set sId$N [lindex [set selId$N] $i]
      }
    set wt [lindex $weight $i]
    set lbli [lindex $lbl $i]
    set fix [lindex $fixSel $i]
    if {($fix == "") || ($fix == "none") || ($fix == 0)} {
      set fix 0} else {set fix 1}
    switch $nSelId {
      2 {lappend plots [distTrajDat $sId1 $sId2 $wt $lbli $seqTime $fix $loSt]}
      3 {lappend plots \
         [angleTrajDat $sId1 $sId2 $sId3 $wt $lbli $seqTime $fix $units $loSt]}
      4 {lappend plots [dihedTrajDat $sId1 $sId2 $sId3 $sId4 $wt $lbli $seqTime $fix $units $loSt]}
      }
    }
  if {($pref == "auto") || ($pref == "default")} {
    set pref "${geomPar}_[join $simNames ","]"
    }
  agrXY $plots $pref \
        title "Geometric parameters for trajectories" \
        xLabel $xLabel yLabel $yLabel loSt $loSt
  puts $loSt "$procName: Done."
  }   ;# geomTraj


#|-proc attribTrajDat {selId {attrib "num"} {lbl "auto"} {seqTime 1}
#|                                                    _ {loSt stdout}} :
#|  -monitors an attribute of an atom selection along a trajectory .
#|  -arguments :
#|    -selId :-selection ID in the selInfo array which is changing in time ;
#|    -attrib :-attribute relative to atom selections changing in time .
#|      -if it is set to "num" (dafault) it follows the number of atoms .
#|      -any other attribute must be considered for the get option of
#|       _ atomselect .
#|      -for "num" attribute the selection is always updated every frame .
#|      -for any attribute but "num" the selection will be updated depending
#|       _ on the "updateSel" key in selInfo array ;
#|    -lbl :-label used for the output data .
#|      -if lbl is set to "auto" an automatic label will be created ;
#|    -seqTime :-use "0" to override time-consecutive trajectories and use
#|     _ frame numbers as X-data instead ;
#|    -loSt :-optional .-output stream for log messages .
#|      -default value :-stdout ;;;
#|  -notes :
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;;
proc attribTrajDat {selId {attrib "num"} {lbl "auto"} {seqTime 1} \
                                                    {loSt stdout}} {
  global trajInfo selInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
  set exclude {}
  puts $loSt "\nMonitoring $attrib for an atom sel along a traj..."
# extracting selection's information
  getSelInfo $selId selTxt molId frag frm firstFrm lastFrm stpFrm upd loadStep \
             selLbl desc
  puts $loSt "trajectory: $trajInfo($molId,name) ($desc)"
  puts $loSt "trajectory range (fisrt:last:step): $firstFrm:$lastFrm:$stpFrm"
  puts $loSt "atom selection ($selId): $selTxt"
  puts $loSt "traj fragments: $frag"
  if {$lbl == "auto"} {
    set lbl "$trajInfo($molId,name)($selId)_$attrib"
    }
  puts $loSt "data set label: $lbl"
  set l_xDat {}
  set l_yDat {}
  trajSelId::setSeqTime $seqTime
  set molSel [atomselect $molId $selTxt]
  if {$attrib == "num"} {
    while {[trajSelId::iterate $selId loSt $loSt]} {
      $molSel frame [getFrame]
      $molSel update
      lappend l_xDat [getDatX]
      lappend l_yDat [$molSel num]
      }
  } else {
    while {[trajSelId::iterate $selId loSt $loSt]} {
      $molSel frame $frm
      if {$upd} {$molSel update}
      lappend l_xDat [getDatX]
      lappend l_yDat [$molSel get $attrib]
      }
    }
  $molSel delete
  return [list $l_xDat $l_yDat $lbl]
  }   ;# attribTrajDat

#|-proc attribTraj {selIds {attrib "num"} {lbl "auto"} {loSt stdout}} :
#|  -monitors an attribute or property of an atom selection for a set of
#|   _ trajectories .
#|  -creates a agr plot including all the data sets .
#|  -arguments :
#|    -selIds :-list of selId's for the selInfo array ;
#|    -attrib :-atom selection's property changing in time .
#|      -if it is set to "num" (default) it follows the number of atoms .
#|      -any other attribute must be considered for the get option of
#|       _ atomselect .
#|      -for the "num" attribute the selection is always updated every frame .
#|      -for any attribute but "num" the selection will be updated depending
#|       _ on the "updateSel" key in selInfo array ;
#|    -lbl :-string used as prefix of output files .
#|      -"auto" can be used to generate an automatic label ;
#|    -loSt :-optional argument .-output stream for log messages .
#|      -default value :-stdout ;;;
#|  -notes :
#|    -the proc has to be updated to manage the new variable 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;;
proc attribTraj {selIds {attrib "num"} {lbl "auto"} {loSt stdout}} {
  set selNames [join $selIds ","]
  puts $loSt "\nMonitoring attribute $attrib for selections: $selNames"
  set seqTime 1
  set l_plots {}
  foreach selId $selIds {
    getSelInfo $selId selTxt id frag frm frst lst stp upd ldstp selLbl desc
    if {![getTrajProp "seqTime" $id $frag]} {set seqTime 0}
    }
  foreach selId $selIds {
    lappend l_plots [attribTrajDat $selId $attrib $selId $seqTime $loSt]
    }
  if {$lbl == "auto"} {set lbl attrib-${attrib}_$selNames}
  if {$seqTime} {set xLbl "Time (ns)"} else {set xLbl "Frame number"}
  agrXY $l_plots $lbl title "Monitoring attribute $attrib for selections" \
        subtitle "Sels: $selNames" xLabel $xLbl yLabel "Attribute $attrib" \
        loSt $loSt
  }   ;# attribTraj


#|-proc gofrTraj {ll_selId1-selId2-lbl args} :
#|  -Calculate atomic radial pair distribution function (g(r)) between two
#|   _ groups of atoms along a trajectory .
#|  -output plots created (for xmgrace .agr and .dat files) :
#|    -g(r) function .
#|    -integral of g(r) .
#|    -unnormalized histogram ;
#|  -vmd gofr function options values defaulted here as :
#|    -delta r :-0.02 ;-r max :-19.98 ;-usepbc :-1 ;;
#|  -arguments :
#|    -required arguments :
#|      -ll_selId1-selId2-lbl :
#|        -list of lists with the format: {{selId1 selId2 label} ...} .
#|        -each set of {selId1 selId2 label} elements represent as a single
#|         _ data set or series in each output graph .
#|        -selId1 :
#|          -selId of the selInfo array specifying a first group of atoms .
#|          -selInfo array keys specification considered :
#|            -"molId", "first", "last", "step", "fragId", "loadStep",
#|             _ and "updateSel" ;
#|            -these specifications are considered for both atom selections ;
#|        -selId2 :
#|          -selId of the selInfo array specifying the second group of atoms .
#|          -selInfo array keys specification considered :
#|            -"selTxt" ;
#|          -acceptable values :
#|            -"", "selId" :-indicates that selId1 will be used for both groups
#|             _ of atoms or selections ;;;
#|        -label :
#|          -data set label to be included in the output graphs .
#|          -acceptable values :
#|            -"", "auto", "full" :
#|              -traj name and selId are included in the label ;
#|            -name :-the name of the traj is used as label ;
#|            -selId, selId1, selId2 :-selId1 and/or selId2 used as label ;
#|            -fragId :-the fragId specification is used as label ;
#|            -selInfo, selInfo1, selInfo2 :
#|              -labels in selInfo array for the selIds used ;;
#|          -default value :-"auto" ;;;;
#|    -args (variable arguments) :
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -dr, delta :-delta r value (resolution) (Angstroms) ;
#|        -default value :-"0.02" ;;
#|      -rmax :-maximum r value (Angstroms) .
#|        -default value :-"19.98" ;;
#|      -usepbc :-flag to turn on periodic boundary conditions (boolean) .
#|        -default value :-"0" ;;;
#|    -args (variable arguments) for plotting set-up :
#|      -"prefix", "pref", "outPrefix", "agrName" :
#|        -prefix for .agr and .dat output file names .
#|        -the rest of the file names are built from other parameters .
#|        -acceptable values :
#|          -the same as the "pref" argument of the agrXY proc ;
#|        -default value :-"gofrTraj" ;;;
#|      -"title", "titl", "plotTitle", "graphTitle" :
#|        -title for the xmgrace plot .
#|        -default value :-"Atomic radial pair distribution function, g(r)" ;;
#|      -"subtitle", "stitl", "plotSubtitle", "graphSubtitle" :
#|        -subtitle for the xmgrace plot .
#|        -default value :-"reference: <refSelId>" ;;
#|      -"xLabel", "labelX", "xlbl", "xLbl" :-label for the x axis (abscissa) .
#|        -default value :-"r (Angstroms)" ;;;;
#|  -notes :
#|    -optional arguments must be specified as pairs of "arg value" keywords .
#|    -the VMD measure gofr selupdate option are specified by the selInfo
#|     _ array key updateSel associated to selId1 .
#|    -global variables :-trajInfo ;
#|    -the label for y-axis in plots cannot be specified as argument ;;
proc gofrTraj {ll_selId1-selId2-lbl args} {
  global trajInfo
# default parameters and proc arguments
  set dr 0.02; set rm 19.98; set pbc 0; set loSt stdout; set outPref "gofrTraj"
  set title "Atomic radial pair distribution function, g(r)"; set subtitle ""
  set xlbl "r (Angstroms)"
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "loSt" -  "channelId" - "log" {set loSt $val}
          "dr" - "delta" {set dr $val}
          "rmax" {set rm $val}
          "usepbc" - "usePbc" - "pbc" {set pbc $val}
          "prefix" - "pref" - "outPrefix" - "agrName" {set outPref $val}
          "title" - "titl" - "plotTitle" - "graphTitle" {set title $val}
          "subtitle" - "stitl" - "plotSubtitle" - "graphSubtitle" {
            set subtitle $val}
          "xLabel" - "labelX" - "xlbl" - "xLbl" {set xlbl $val}
          default {puts $loSt "gofrTraj: $arg argument unkown"; return ""}
          }
        }
      }
  } else {   ;# odd number of arguments
    puts $loSt "gofrTraj: Odd number of optional arguments! args: $args"
    return ""
    }
# log info header
  puts $loSt "\ngofrTraj: Atomic radial pair distribution function, g(r),"
  puts $loSt "          between two sets of atoms along a trajectory."
  puts $loSt "requested selIds: ${ll_selId1-selId2-lbl}"
# initializing data sets
  set gofr_plotDat {}; set gofr-int_plotDat {}; set gofr-hist_plotDat {}
  foreach l_set ${ll_selId1-selId2-lbl} {   ;# for each data series
    set selId1 [lindex $l_set 0]
    set selId2 [lindex $l_set 1]
    set lblSet [lindex $l_set 2]
    if {($selId2 == "") || ($selId2 == "selId")} {set selId2 $selId1}
    if {$lblSet == ""} {set lblSet "auto"}
    puts $loSt "selId1: $selId1   selId2: $selId2   lblSet: $lblSet"
# extracting information from selId2 (only selTxt2 is considered in this case)
    getSelInfo $selId2 selTxt2 id frag frm fst lst stp upd ldstp selLbl2 desc
# extracting information from selId1 (selTxt1, id, fragId, fst, lst, stp, upd)
    getSelInfo $selId1 selTxt1 id frag frm fst lst stp upd ldstp selLbl1 desc
# load traj
    if {[getTrajProp "loaded" $id $frag]} {set loadT 0} else {set loadT 1}
    if {$loadT} {trajLoad $id $frag loadStep $ldstp loSt $loSt}
# convert sim time (ns) into frame number
    set frag [trajFragSpec $frag $id loSt $loSt]
    set fst [simFrame $fst $id $frag]
    set lst [simFrame $lst $id $frag]   ;# note that simFrame is not exact
    if {($lst == -1) || ($lst >= [molinfo $id get numframes])} {
      set lst [getTrajProp "finFrame" $id $frag]}   ;# maybe useless
    switch $lblSet {
      "name" {set lblSet $trajInfo($id,name)}
      "selId" {set lblSet "$selId1-$selId2"}
      "selId1" {set lblSet $selId1}
      "selId2" {set lblSet $selId2}
      "fragId" {set lblSet $frag}
      "selInfo" {set lblSet "$selLbl1-$selLbl2"}
      "selInfo1" {set lblSet $selLbl1}
      "selInfo2" {set lblSet $selLbl2}
      "" - "auto" - "full" {set lblSet "$trajInfo($id,name)_$selId1-$selId2"}
      }
# log information
    puts $loSt "\nselection 1: $selTxt1 \nselection 2: $selTxt2"
    puts $loSt "mol Id: $id; name: $trajInfo($id,name) - $trajInfo($id,desc)"
    puts $loSt "delta: $dr; rmax: $rm; usepbc: $pbc; selupdate: $upd;"
    puts $loSt "first: $fst; last: $lst; step: $stp; label: $lblSet"
# calculates g(r)
    set res [nI_gofr $id $selTxt1 $selTxt2 $dr $rm $pbc $upd $fst $lst $stp]
    puts $loSt "Writing output: g(r), integral and unnormalized histogram..."
    set rDat [lindex $res 0]
    set gofrDat [lindex $res 1]
    set intDat [lindex $res 2]
    set unnormHistDat [lindex $res 3]
    puts $loSt "Number of frames processed and skipped: [lindex $res 4]"
    lappend gofr_plotDat [list $rDat $gofrDat $lblSet]
    lappend gofr-int_plotDat [list $rDat $intDat $lblSet]
    lappend gofr-hist_plotDat [list $rDat $unnormHistDat $lblSet]
    if {$loadT} {trajDelete $id $frag loSt $loSt}
    }
# setting up the output label for files and legends
  set tSubt $subtitle
  agrXY ${gofr_plotDat} $outPref \
        title $title subtitle $tSubt xlbl $xlbl \
        ylbl "g(r)" loSt $loSt
  if {$subtitle == ""} {set tSubt "g(r) - Integral"}
  agrXY ${gofr-int_plotDat} "${outPref}.int" \
        title $title subtitle $tSubt xlbl $xlbl \
        ylbl "g(r) Integral" loSt $loSt
  if {$subtitle == ""} {set tSubt "g(r) - Unnormalized histogram"}
  agrXY ${gofr-hist_plotDat} "${outPref}.hist" \
        title $title subtitle $tSubt xlbl $xlbl \
        ylbl "g(r) Histogram" loSt $loSt
  }


#|-radGyrTrajDat {selId {lbl "auto"} {loSt stdout}} :
#|  -calculates of the radius of gyration of a selection of atoms along a
#|   _ trajectory .
#|  -returns a list of lists containing the data for plotting with agrXY .
#|  -arguments :
#|    -selId :
#|      -atom selection name of the atoms to be considered in the calculation .
#|      -must be included in the global array selInfo ;
#|    -lbl :
#|      -text label included in the output data .
#|      -"auto" (default) may be specified to generate an automatic label ;
#|    -loSt :-optional output stream for log messages .
#|      -default value :-stdout ;;
#|  -notes :
#|    -in the case that the trajectory was not already loaded it will be loaded
#|     _ with trajLoad and at the end will be deleted using trajDelete .
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;;
proc radGyrTrajDat {selId {lbl "auto"} {weight "mass"} {loSt stdout}} {
  global selInfo trajInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
  set exclude {}
  puts $loSt "\n+++ Calculating the radius of gyration for selection $selId..."
# extracting information from the selInfo array
  getSelInfo $selId selTxt molId frag frm firstFrm lastFrm stepFrm upd \
             loadStep selLbl desc
  if {$lbl == "auto"} {
    set lbl "$trajInfo($molId,name)($selId)"
    }
# reporting info and returning results
  puts $loSt "trajectory: $trajInfo($molId,name) ($trajInfo($molId,desc))"
  puts $loSt "trajectory range (fisrt:last:step): $firstFrm:$lastFrm:$stepFrm"
  puts $loSt "atom selection ($selId): $selTxt"
  puts $loSt "weighting scheme : $weight"
  puts $loSt "data set label: $lbl"
# calculating the radius of gyration along the trajectory
  set sel [atomselect $molId $selTxt]
  while {[trajSelId::iterate $selId excl $exclude loSt $loSt]} {
    $sel frame [getFrame]
    lappend l_xDat [getDatX]
    lappend l_yDat [measure rgyr $sel weight [$sel get $weight]]
    }
  $sel delete
  return [list $l_xDat $l_yDat $lbl]
  }

#|-proc radGyrTraj {selIds {lbl "auto"} {weight "mass"} {loSt stdout}} :
#|  -calculates the radius of gyration for a set of trajectories .
#|  -each trajectory may consider different molecules or atom selections .
#|  -creates a xmgrace agr plot file considering all the trajectories .
#|  -arguments :
#|    -l_selId :
#|      -list for selection IDs corresponding to each trajectory .
#|      -each selId must be included in the selInfo array  ;
#|    -lbl :
#|      -label included in output file name .
#|      -if set to "auto", an automatic label will be created ;
#|    -loSt :-optional output stream for log messages .
#|      -default value :-stdout ;;;;
proc radGyrTraj {selIds {lbl "auto"} {weight "mass"} {loSt stdout}} {
  puts $loSt "\nradGyrTraj: Radius of gyration for traj: $selIds"
  foreach selId $selIds {
    lappend l_plots [radGyrTrajDat $selId $selId $weight $loSt]
    }
  set sels [join $selIds ","]
  if {$lbl == "auto"} {
    set lbl "radGyrTraj_$sels"
    }
  puts $loSt "Creating radius of gyration plot for selections $sels..."
  agrXY $l_plots $lbl \
        title "Radius of gyration" subtitle "selections: $sels" \
        xLabel "Time (ns)" yLabel "Radius of gyration (Angstroms)" loSt $loSt
  puts $loSt "\nradGyrTraj: Done."
  }


#|-proc sasaTrajDat {selId pSelId {lbl "auto"} {srad 1.4} {loSt stdout}} :
#|  -calculates the solvent-accesible surface area along a trajectory .
#|  -return a list of lists for plotting the data with agrXY .
#|  -arguments :
#|    -selId :-name id of a selection stored in the selInfo array .
#|      -this selection is applied in the "-restrict" option in the sasa proc
#|       _ of VMD, as in ref r1 .
#|      -is supposed to be the selection for which the calculation is done .
#|      -the parameters first, last, update, etc. used in the calculation
#|       _ are taken from this selection ;
#|    -pSelId :-id of the bigger selection (typically the whole protein) .
#|      -the mol id (for the moment) must be the same as that in selId ;
#|    -lbl :-label for the data set returned .
#|      -"auto" (dafault) may be specified to generate an automatic label ;
#|    -srad :-sphere radius for the measure sasa VMD procedure (see r1) .
#|      -default value :-1.4 ;;
#|    -loSt :-output stream for log messages ;;
#|  -notes :
#|    -is pending to consider the "points" and "samples" options
#|     _ in the sasa procedure of VMD .
#|    -it is not clear if the selection and the restrict selection may be from
#|     _ different molecules .
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;
#|  -references :
#|    -r1 :-script on internet, http://www.ks.uiuc.edu/Research/vmd/
#|     _mailing_list/vmd-l/att-18670/sasa.tcl ;;;
proc sasaTrajDat {selId pSelId {lbl "auto"} {srad 1.4} {loSt stdout}} {
  global selInfo trajInfo
  namespace import trajSelId::getFrame trajSelId::getDatX
  set exclude {}
  puts $loSt "\n*** sasaTraj: Calculating SASA for a trajectory."
# extracting info
  getSelInfo $pSelId pSelTxt pMolId pFrag frm fst lst stp pUpd ldStp selLbl desc
  getSelInfo $selId selTxt molId frag frm fst lst stp upd ldStp selLbl desc
  puts $loSt "Restricted selection:"
  puts $loSt "selId: $selId ($selTxt)"
  puts $loSt "molId: $molId: $trajInfo($molId,name) - $trajInfo($molId,desc)"
  showSelInfo $selId $loSt
  puts $loSt "Protein selection:"
  puts $loSt "selId: $pSelId ($pSelTxt)"
  if {$lbl == "auto"} {set lbl "$trajInfo($pMolId,name)_sel-${selId}($pSelId)"}
# calculating for each frame
  set selRest [atomselect $molId $selTxt]
  set selProt [atomselect $pMolId $pSelTxt]
  set yDat {}
  set xDat {}
  while {[trajSelId::iterate $selId excl $exclude loSt $loSt]} {
    $selRest frame [getFrame]
    $selProt frame [getFrame]
    if {$upd} {$selRest update}
    if {$pUpd} {$selProt update}
    lappend yDat [measure sasa $srad $selProt -restrict $selRest]
    lappend xDat [getDatX]
    }
  $selRest delete
  $selProt delete
  puts $loSt "\nsasaTrajDat: Done."
  return [list $xDat $yDat $lbl]
}   ;# sasaTrajDat

#|-proc sasaTraj {l_selId-pSelId {lbl "auto"} {srad 1.4} {loSt stdout}} :
#|  -calculate the solvent-accesible surface area for a set of trajectories .
#|  -makes an xmgrace graph with all the data sets .
#|  -see sasaTrajDat or options and configuration .
#|  -arguments :
#|    -ll_selId-pSelId :-list of lists with pairs of name ids of selections
#|     _ stored in the selInfo array (as applied in the sasaTrajDat proc) .
#|      -format :-{{selId pSelId} ...} ;
#|      -each selId is applied in the "-restrict" option in the sasa proc
#|       _ of VMD, as in ref r1 .
#|      -it is supposed to be the selection for which the calculation is done .
#|      -the parameters first, last, update, etc. used in the calculation
#|       _ are taken from this selection .
#|      -pSelId ids for bigger selections (typically the whole protein) ;
#|    -lbl :-label prefix for the output files .
#|      -"auto" (default) may be specified to use automatic name ;
#|    -srad :-sphere radius (see sasaTrajDat) ;
#|    -loSt :-output stream for log messages ;;
#|  -notes :
#|    -pending to improve the output/log info and the automatic file name ;;
proc sasaTraj {ll_selId-pSelId {lbl "auto"} {srad 1.4} {loSt stdout}} {
  global selInfo trajInfo
  puts $loSt "\n*** sasaTraj: Calculating SASA for a set of trajectories."
  set selNames ""
  foreach selPair ${ll_selId-pSelId} {
    lassign $selPair sel pSel
    append selNames "${sel}($pSel)" ","
    lappend pDat [sasaTrajDat $sel $pSel $sel $srad $loSt]
    }
  if {$lbl == "auto"} {set lbl "$selNames"}
  agrXY $pDat $lbl title "Solvent-accesible surface area" \
              subtitle "selections: ${selNames}" \
              xLabel "Time (ns)" yLabel "SASA (Amstrongs^2)" loSt $loSt

  puts $loSt "\nsasaTraj: Done."
  }

#|-proc secStructTraj {selIdl {l_ssCode {T E B H G I C}} {loSt stdout}} :
#|  -calculates the incidence of STRIDE secundary structures for all residues
#|   _ in a selection along a trajectory .
#|  -the evolution of each sec Struct is plotted in an .agr file .
#|  -secondary structure codes :
#|    -T - turn .
#|    -E - extended conformation (beta sheets) .
#|    -B - isolated bridge .
#|    -H - alpha helix .
#|    -G - 3-10 helix .
#|    -I - pi helix .
#|    -C - coil ;
#|  -arguments :
#|    -selIdl :-list of atom selection ids for the selInfo array including the
#|     _ portion of the protein to be analyzed .
#|      -only the atoms with 'name CA' will be considered .
#|      -the first, last and step selInfo keys are considered ;
#|    -l_ssCode :-list of sec struct 1-letter codes to be considered .
#|      -default value :-{T E B H G I C} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;;;
#|  -notes :
#|    -if the trajectory is not already loaded it will temporarily be loaded .
#|    -if several trajs are considered, the coord would have to be not
#|     _ already loaded :
#|      -do not trajLoad or use trajDelete for each traj before calling
#|       _ this proc ;
#|    -for each selId in selIdl when considering more than one traj, the coords
#|     _ will be loaded and deleted every time :
#|      -this will be avoided when multiple selIds in the selInfo array is
#|       _ supported ;
#|    -there is still a problem with the codeX incidence to solve .
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;;
proc secStructTraj {selIdl {l_ssCode {T E B H G I C}} {loSt stdout}} {
  global trajInfo selInfo fibVMolId tclScriptPath fibVSaveTrimer
  namespace import trajSelId::getFrame trajSelId::getDatX
  set exclude {}
  puts $loSt "\n+++ secStructTraj: Incidence of STRIDE secundary structure motifs of proteins along a trajectory. +++"
  puts $loSt "\nlist of selection Ids: $selIdl"
  set prevTop [molinfo top]   ;# stores the original top molecule
  foreach selId $selIdl {
    getSelInfo $selId selTxt molId frag frm fst lst stp upd ldStp selLbl desc
    showSelInfo $selId $loSt
    puts $loSt "Only CA atoms will be considered."
    showTrajInfo $molId $frag loSt $loSt
    puts $loSt "\ntrajectory time interval: $fst - $lst load step: $ldStp"
    set sel [atomselect $molId "$selTxt and name CA"]
    set nRes [$sel num]
    mol top $molId
    foreach ssCode $l_ssCode {
      set dataY($selId-$ssCode) {}   ;# data lists for plotting
      }
    set dataY($selId-X) {}
    set dataX($selId) {}   ;# x data for time axis
    set useCodeX 0   ;# to indicate whether 'X' ss code was used during the traj
# consider the special cas of the periodic fibril model fibV
    if {$trajInfo($molId,name) == "bao2i-s1"} {
      set fibVMolId $molId
      set fibVSaveTrimer 1
      set calcSecStructV 1
    } else {
      set calcSecStructV 0
      }
    puts $loSt "\nCollecting secondary structure data for motifs: $l_ssCode "
#    puts -nonewline $loSt "frame (out of $lst):"
    while {[trajSelId::iterate $selId excl $exclude loSt $loSt]} {
#      puts -nonewline $loSt " [getFrame]"
      foreach ssCode $l_ssCode {set inc($selId-$ssCode) 0.0}
      set inc($selId-X) 0.0
      animate goto [getFrame]
      $sel frame [getFrame]
      if {$calcSecStructV} {
        source "${tclScriptPath}anMD/tools/fibVsecStructCA.tcl"
      } else {
        mol ssrecalc $molId
        }
      foreach ssCodeRes [$sel get structure] {
        set codeFound 0
        foreach ssCode $l_ssCode {
          if {$ssCodeRes == $ssCode} {
            set inc($selId-$ssCode) [expr $inc($selId-$ssCode) + 1.0]
            set codeFound 1
            }
          }
        if {$codeFound == 0} {
          set inc($selId-X) [expr $inc($selId-X) + 1.0]
          set useCodeX 1
          }
        }
      foreach ssCode $l_ssCode {
        lappend dataY($selId-$ssCode) [expr $inc($selId-$ssCode)/$nRes]
        }
      if {$useCodeX} {lappend dataY($selId-X) [expr $inc($selId-X)/$nRes]}
      lappend dataX($selId) [getDatX]
      }
    $sel delete
    }
#  puts $loSt "\file prefix: ssTraj-$trajInfo($molId,name)-[join $selIdl ","]"
  set ll_plots {}
  foreach selId $selIdl {
    foreach ssCode $l_ssCode {
      lappend ll_plots [list $dataX($selId) $dataY($selId-$ssCode) "$selId-ss-$ssCode"]
      }
    if {$useCodeX} {
      lappend ll_plots [list $dataX($selId) $dataY($selId-X) "$selId-ss-others"]
      }
    }
  agrXY $ll_plots \
        "ssTraj-$trajInfo($molId,name)-[join $l_ssCode ""]-[join $selIdl ","]" \
        title "Incidence of secondary structure along a trajectory" \
        subtitle "selIds: [join $selIdl ","]" \
        xLabel "Time (ns)" \
        yLabel "Incidence" \
        loSt $loSt
  mol top $prevTop
  puts $loSt "\nsecStructTraj: Done!"
  }   ;# *** secStructTraj ***

#|-proc secStructAv {selId {attrib "resid"} {loSt stdout}} :
#|  -calculates the time average of the incidence of secondary structure in a
#|   _ trajectory for different categories of protein structure in a selection .
#|  -the categories may refer to resids, resnames, segnames, fragments, etc. .
#|  -the results will be output as a matrix .
#|  -secondary structure codes :
#|    -T - turn .
#|    -E - extended conformation (beta sheets) .
#|    -B - isolated bridge .
#|    -H - alpha helix .
#|    -G - 3-10 helix .
#|    -I - pi helix .
#|    -C - coil ;
#|  -arguments :
#|    -selId :-selection id including the residues to be analyzed .
#|      -only CA atoms would be considered ;
#|    -attrib :-attribute to be analyzed (see atomselect keywords) .
#|      -acceptable values :
#|        -resid, resname, segid, segname, chain, ... ;
#|      -default value :-resid ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;;;
#|  -notes :
#|    -if the trajectory is not already loaded, it will be so temporarily .
#|    -each resid will be considered regardless the different segids where
#|     _ they may appear .
#|    -the list of sec struct code may be set as an argument, including the
#|     _ use of the "others" secundary structure result like in secStructTraj .
#|    -the proc has to be updated to manage the new variables 'selLbl' and
#|     _ 'desc' from proc 'getSelInfo' ;;
proc secStructAv {selId {attrib "resid"} {loSt stdout}} {
  global trajInfo selInfo fibVMolId tclScriptPath fibVSaveTrimer
  namespace import trajSelId::getFrame trajSelId::getDatX
  set exclude {}
  puts $loSt "\n+++ anMD_secStructResAv: Average incidence of STRIDE secondary structure per residue along a trajectory +++"
  set prevTop [molinfo top]   ;# stores the original top molecule
  getSelInfo $selId selTxt molId frag frm fst lst stp upd ldStp selLbl desc
  puts $loSt "\nSelection id: $selId"
  puts $loSt "atomselect attribute to be analyzed: $attrib"
  showSelInfo $selId $loSt
  puts $loSt "Only CA atoms will be considered."
  showTrajInfo $molId $frag loSt $loSt
  puts $loSt "\ntrajectory time interval: $fst - $lst  load step: $ldStp"
  if {$trajInfo($molId,name) == "bao2i-s1"} {
    set fibVMolId $molId
    set fibVSaveTrimer 1
    set calcSecStructV 1
  } else {
    set calcSecStructV 0
    }
  set sel [atomselect $molId "$selTxt and name CA"]
  set resIndL [$sel get index]
  set attribIdL {}   ;# init list with resids included (not repeated)
  foreach ind $resIndL {
    set ca [atomselect $molId "index $ind"]
    set atId [$ca get $attrib]
    $ca delete
    if {[lsearch $attribIdL $atId] == -1} {
      lappend attribIdL $atId
      set freq($atId) 1
    } else {
      incr freq($atId)
      }
    }   ;# list of residues (not repeated), and the frequency of each each resid
  puts $loSt "\nlist of elements found (in attribute: $attrib): $attribIdL"
  foreach atId $attribIdL {
    foreach ssCode {T E B H G I C} {
      set av($atId,$ssCode) 0.0   ;# init sum of time average incidences
      }
    }
  set nFrm 0
  while {[trajSelId::iterate $selId excl $exclude loSt $loSt]} {
    animate goto [getFrame]
    $sel frame [getFrame]
    if {$calcSecStructV} {
      source "${tclScriptPath}anMD/tools/fibVsecStructCA.tcl"
    } else {mol ssrecalc $molId}
    foreach atId $attribIdL {
      foreach ssCode {T E B H G I C} {
        set inc($atId,$ssCode) 0.0   ;# init incidence of sec struct per attrib
        }
      }
    foreach ind $resIndL {
      set ca [atomselect $molId "index $ind" frame [getFrame]]
      set atId [$ca get $attrib]
      set secStruct [$ca get structure]
      $ca delete
      set inc($atId,$secStruct) [expr $inc($atId,$secStruct) + 1.0]
      }   ;# incidence of each sec struct code for each attribute and each frame
    foreach ssCode {T E B H G I C} {
      foreach atId $attribIdL {
        set av($atId,$ssCode) [expr $av($atId,$ssCode) + \
                                    $inc($atId,$ssCode)/$freq($atId)]
        }
      }   ;# calculate time average for each sec struct incidence per attribute
    incr nFrm
    }
  $sel delete
  puts -nonewline $loSt "\n\n$attrib"
  foreach ssCode {T E B H G I C} {puts -nonewline $loSt "\tss-$ssCode"}
  puts $loSt ""
  foreach atId $attribIdL {
    puts -nonewline $loSt "$atId"
    foreach ssCode {T E B H G I C} {
      puts -nonewline $loSt [format "\t%.4f" [expr $av($atId,$ssCode)/$nFrm]]
      }
    puts $loSt ""
    }
  mol top $prevTop
  puts $loSt "\nsecStructAv: Done."
  }   ;# *** secStructAv ***


# +++ trajInfo- and selInfo-independent procedures +++

#|-proc nI_rmsdFrame {molID molFrm molSelTxt refID refFrm refSelTxt} :
#|  -Returns the RMSD (mass weighted) for two selections of atoms .
#|  -the molecule will be alligned over the reference .
#|  -both selections must consider the same number of atoms .
#|  -used for a single frame .
#|  -it will create and delete the atom selection .
#|  -arguments :
#|    -molID, molFrm, molSelTxt :-ID, frame number and sel text of the molec ;
#|    -refID, refFrm, refSelTxt :-ID, frame number and sel text of the ref ;;
proc nI_rmsdFrame {molID molFrm molSelTxt refID refFrm refSelTxt} {
  set refSel [atomselect $refID $refSelTxt frame $refFrm]
  set molSel [atomselect $molID $molSelTxt frame $molFrm]
  set molTotSel [atomselect $molID all frame $molFrm]
  $molTotSel move [measure fit $molSel $refSel]
  set rVal [measure rmsd $molSel $refSel weight [$molSel get mass]]
  $refSel delete
  $molSel delete
  $molTotSel delete
  return $rVal
  }

#|-proc nI_gofr {molID selTxt1 selTxt2 dlt rmax pbc selupd fst lst stp} :
#|  -return the lists returned by the vmd command measure gofr .
#|  -all options of the command can be adjusted .
#|  -arguments :
#|    -molID :- ;
#|    -selTxt1,selTxt2 :- ;
#|    -dlt :- ;
#|    -rmax :- ;
#|    -pbc :- ;
#|    -selupd :- ;
#|    -fst :- ;
#|    -lst :- ;
#|    -stp :- ;;;
proc nI_gofr {molID selTxt1 selTxt2 dlt rmax pbc selupd fst lst stp} {
  set sel1 [atomselect $molID $selTxt1]
  set sel2 [atomselect $molID $selTxt2]
  set rVal [measure gofr $sel1 $sel2 delta $dlt rmax $rmax usepbc $pbc \
                       selupdate $selupd first $fst last $lst step $stp]
  $sel1 delete
  $sel2 delete
  return $rVal
  }

#|-proc nI_hbCount {molId selTxtD selTxtA cutoff angle} :
#|  -return the nomber of hydogen bond found between two atom selections .
#|  -each atom selection is combined with other default selections to consider
#|   _ only real hydrogen donors or acceptors .
#|  -only the current frame of the molecule is considered .
#|  -a list of the indices of the atoms involved in the HB are printed .
#|  -arguments :- ;;
proc nI_hbCount {molId selTxtD selTxtA cutoff angle} {
  global selInfo
  set sel1 [atomselect $molId "($selTxtD) and ($selInfo(HBd,selTxt))"]
  set sel2 [atomselect $molId "($selTxtA) and ($selInfo(HBa,selTxt))"]
  set res [measure hbonds $cutoff $angle $sel1 $sel2]
  $sel1 delete
  $sel2 delete
  set rest [lassign $res l_donor l_acceptor l_hydrogen]
  set nHB [llength $l_donor]
  puts "\n $nHB hydrogen bonds in molecule $molId: [molinfo $molId get name]:"
  puts " Donors: \t Acceptors: \t Hydrogens:"
  puts "-------------------------------------------------------------------"
  for {set i 0} {$i < $nHB} {incr i} {
    puts "$i:\t[lindex $l_donor $i] \t [lindex $l_acceptor $i] \t\
               [lindex $l_hydrogen $i]"
    }
  return $nHB
  }

#|-proc nI_sasa {molId selTxt protSelTxt sRad} :
#|  -return the solvent-accesible surface area for a protein or section
#|   _ of a protein using the restrict option of the sasa proc in VMD .
#|  -arguments :
#|    -molId :-VMD molecule id already loaded ;
#|    -selTxt :-restrict selection considered for the calculation ;
#|    -protSelTxt :-global selection (see sasaTrajDat and its refs) ;
#|    -sRad :-sphere radius added to each atomic radius to find the points on
#|     _ a sphere that are exposed to solvent ;;;
proc nI_sasa {molId selTxt protSelTxt sRad} {
  set selRest [atomselect $molId $selTxt]
  set selProt [atomselect $molId $protSelTxt]
  set res [measure sasa $sRad $selProt -restrict $selRest]
  $selRest delete
  $selProt delete
  return $res
  }

#|-proc frame2time {frm trajLength numFrames} :
#|  -As simTime returns the time in ns equivalent to a frame number .
#|  -the trajectory length and the number of frames are not calculated but
#|   _ taken as arguments instead .
#|  -in principle it is not necessary to have the trajectory loaded in vmd .
#|  -arguments :
#|    -frm :-frame number ;
#|    -trajLength -Length of the trajectory in ns
proc frame2time {frm trajLength numFrames} {
  return [expr $frm*$trajLength/($numFrames-1.0)]
  }

#|-proc time2frame {nsTime trajLength numFrames} :
#|  -return the nearest frame equivalent to a specific tim in ns .
#|  -the result will be rounded to an integer ;
proc time2frame {nsTime trajLength numFrames} {
  return [expr round($nsTime*($numFrames-1)/$trajLength)]
  }


# +++ Utility procedures +++

#|-proc agrXY {ll_xd-yd-lbl prefix args} :
#|  -writes an .agr file for xmgrace plotting for one or more sets of xy data .
#|  -writes a .dat file for each of the data sets included .
#|  -arguments :
#|    -ll_xd-yd-lbl :
#|      -list of lists of xd-yd-lbl data sets {{l_xDat l_yDat lbl} ...} .
#|      -allows multiple sets to be plotted on the same graph ;
#|    -prefix :-Name prefix (without extension) used for agr and dat files .
#|      -the .dat file name also includes the corresponding data set label ;
#|    -args (variable arguments) :
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -"title", "titl", "plotTitle", "graphTitle" :
#|        -title for the xmgrace plot .
#|        -default value :-"Title: " ;;
#|      -"subtitle", "stitl", "plotSubtitle", "graphSubtitle" :
#|        -subtitle for the xmgrace plot .
#|        -default value :-"" ;;
#|      -"xLabel", "labelX", "xlbl", "xLbl" :-label for the x axis (abscissa) .
#|        -default value :-"x Label" ;;
#|      -"yLabel", "labelY", "ylbl", "yLbl" :-label for the y axis (ordinate) .
#|        -default value :-"y Label" ;;;
#|    -loSt :-output stream for log messages ;;
proc agrXY {ll_xd-yd-lbl prefix args} {
  global anMD_graphicsOn
# default values for arguments
  set procName [lindex [info level 0] 0]; set loSt stdout
  set titl "Title: "; set stitl ""; set xlbl "x Label"; set ylbl "y Label"
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "title" - "titl" - "plotTitle" - "graphTitle" - "pTitle" {
            set titl $val}
          "subtitle" - "stitl" - "plotSubitle" - "graphSubtitle" - "pSubtitle" {
            set stitl $val}
          "xLabel" - "labelX" - "xlbl" - "xLbl" {set xlbl $val}
          "yLabel" - "labelY" - "ylbl" - "yLbl" {set ylbl $val}
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {puts $loSt "$procName: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "$procName: Odd number of variable arguments! args: $args"
      return ""
      }
  puts $loSt "\nCreating a XY-data agr file for xmgrace..."
# opening agr file and writing header
  set agrFile [open $prefix.agr w]
  puts $agrFile "@map color 0 to (255, 255, 255), \"white\""
  puts $agrFile "@map color 1 to (0, 0, 0), \"black\""
  puts $agrFile "@map color 2 to (0, 0, 255), \"blue\""
  puts $agrFile "@map color 3 to (255, 0, 0), \"red\""
  puts $agrFile "@map color 4 to (0, 139, 0), \"green4\""
  puts $agrFile "@map color 5 to (255, 165, 0), \"orange\""
  puts $agrFile "@map color 6 to (148, 0, 211), \"violet\""
  puts $agrFile "@map color 7 to (220, 220, 220), \"grey\""
  puts $agrFile "@map color 8 to (255, 0, 255), \"magenta\""
  puts $agrFile "@map color 9 to (188, 143, 143), \"brown\""
  puts $agrFile "@map color 10 to (64, 224, 208), \"turquoise\""
  puts $agrFile "@map color 11 to (0, 255, 255), \"cyan\""
  puts $agrFile "@map color 12 to (0, 255, 0), \"green\""
  puts $agrFile "@map color 13 to (255, 255, 0), \"yellow\""
  puts $agrFile "@map color 14 to (114, 33, 188), \"indigo\""
  puts $agrFile "@map color 15 to (103, 7, 72), \"maroon\""
  puts $agrFile "@type xy"
  puts $agrFile "@title \"$titl\""
  puts $agrFile "@title size 1.50"
  puts $agrFile "@subtitle \"$stitl\""
  puts $agrFile "@subtitle size 1.20"
  puts $agrFile "@xaxis bar linewidth 2.5"
  puts $agrFile "@xaxis label \"$xlbl\""
  puts $agrFile "@xaxis label char size 1.80"
  puts $agrFile "@xaxis tick major linewidth 2.5"
  puts $agrFile "@xaxis tick minor linewidth 2.5"
  puts $agrFile "@xaxis ticklabel char size 1.60"
  puts $agrFile "@yaxis bar linewidth 2.5"
  puts $agrFile "@yaxis label \"$ylbl\""
  puts $agrFile "@yaxis label char size 1.80"
  puts $agrFile "@yaxis tick major linewidth 2.5"
  puts $agrFile "@yaxis tick minor linewidth 2.5"
  puts $agrFile "@yaxis ticklabel char size 1.60"
  puts $agrFile "@legend box linewidth 2.5"
  puts $agrFile "@legend char size 1.0000"
  puts $loSt "Writing XY data on file $prefix.agr..."
  puts $loSt "title: $titl"
  puts $loSt "subtitle: $stitl"
  puts $loSt "x and y labels: $xlbl, $ylbl"
# writing data for each one of the sets
  set s 0
  foreach dataSet ${ll_xd-yd-lbl} {
    set l_xDat [lindex $dataSet 0]
    set l_yDat [lindex $dataSet 1] 
    set setName [lindex $dataSet 2]
    set nDat [llength $l_xDat]
    set datFile [open ${prefix}_$setName.dat w]
    puts $agrFile "@s$s line linewidth 3.0"
    puts $agrFile "@s$s legend \"$setName\""
    puts $loSt "Including data set $setName and writing ${prefix}_$setName.dat..."
    for {set d 0} {$d < $nDat} {incr d} {
      puts $agrFile "[lindex $l_xDat $d] [lindex $l_yDat $d]"
      puts $datFile [lindex $l_xDat $d]\t[lindex $l_yDat $d]
      }
    puts $agrFile "&"
    close $datFile
    incr s
    }
  close $agrFile
  if {$anMD_graphicsOn} {exec xmgrace $prefix.agr &}
  }

#|-proc strNames {l_name} :
#|  -Return a string with the comma-separated names ;
proc strNames {l_name} {
  set names ""
  foreach name $l_name {
    if {$names == ""} {set names $name} else {append names ,$name}}
  return $names
  }

#|-proc infoAtom {ind molId} :
#|  -return a string with the information associated to an atom .
#|  -output format: molName.segName.resName-resId.name.index ;
#|  -requires the script aaCode001.tcl (already sourced) .
#|  -arguments :
#|    -ind :
#|      -VMD index of the atom in the loaded molecule ;
#|    -molId :
#|      -Id of the loaded molecule in VMD ;;;
proc infoAtom {ind molId} {
  set atSel [atomselect $molId "index $ind"]
  set molName [molinfo $molId get name]
  set segName [$atSel get segname]
  if {$segName == "{}"} {set segName [$atSel get chain]}
  set resName [$atSel get resname]
  set resId [$atSel get resid]
  set atName [$atSel get name]
  $atSel delete
  return "$molName.$segName.[aaCode_capToCam $resName]$resId.$atName.$ind"
  }

#|-proc aaCode_1to3 {resname} :
#|  -return an aminoacid 3-letter code from a 1-letter code resname ;
proc aaCode_1to3 {resname} {
  switch $resname {
    "A" {return "ALA"}
    "R" {return "ARG"}
    "N" {return "ASN"}
    "D" {return "ASP"}
    "B" {return "ASX"}
    "C" {return "CYS"}
    "E" {return "GLU"}
    "Q" {return "GLN"}
    "Z" {return "GLX"}
    "G" {return "GLY"}
    "H" {return "HIS"}
    "I" {return "ILE"}
    "L" {return "LEU"}
    "K" {return "LYS"}
    "M" {return "MET"}
    "F" {return "PHE"}
    "P" {return "PRO"}
    "S" {return "SER"}
    "W" {return "TRP"}
    "T" {return "THR"}
    "Y" {return "TYR"}
    "V" {return "VAL"}
    default {return $resname}
    }
  }

#|-proc aaCode_3to1 {resname} :
#|  -return an aminoacid 1-letter code from a 3-letter code resname .
#|  -the input must be in capitals ;
proc aaCode_3to1 {resname} {
  switch $resname {
    "ALA" { return "A" }
    "ARG" { return "R" }
    "ASN" { return "N" }
    "ASP" { return "D" }
    "ASX" { return "B" }
    "CYS" { return "C" }
    "GLU" { return "E" }
    "GLN" { return "Q" }
    "GLX" { return "Z" }
    "GLY" { return "G" }
    "HIS" { return "H" }
    "HSE" { return "H" }
    "HSD" { return "H" }
    "ILE" { return "I" }
    "LEU" { return "L" }
    "LYS" { return "K" }
    "MET" { return "M" }
    "PHE" { return "F" }
    "PRO" { return "P" }
    "SER" { return "S" }
    "THR" { return "T" }
    "TRP" { return "W" }
    "TYR" { return "Y" }
    "VAL" { return "V" }
    default { return $resname }
    }
  }

#|-proc aaCode_capToCam {resname} :
#|  -return an aminoacid 3-letter code in cammel style from an aa 3-letter code
#|   _ resname in capitals ;
proc aaCode_capToCam {resname} {
  switch $resname {
    "ALA" { return "Ala" }
    "ARG" { return "Arg" }
    "ASN" { return "Asn" }
    "ASP" { return "Asp" }
    "ASX" { return "AsX" }
    "CYS" { return "Cys" }
    "GLU" { return "Glu" }
    "GLN" { return "Gln" }
    "GLX" { return "GlZ" }
    "GLY" { return "Gly" }
    "HIS" { return "His" }
    "HSE" { return "Hse" }
    "HSD" { return "Hsd" }
    "ILE" { return "Ile" }
    "LEU" { return "Leu" }
    "LYS" { return "Lys" }
    "MET" { return "Met" }
    "PHE" { return "Phe" }
    "PRO" { return "Pro" }
    "SER" { return "Ser" }
    "THR" { return "Thr" }
    "TRP" { return "Trp" }
    "TYR" { return "Tyr" }
    "VAL" { return "Val" }
    default { return $resname }
    }
  }


# +++ testing procedures +++

