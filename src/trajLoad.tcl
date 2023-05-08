#|-proc trajLoad {id {l_fragId "name"} args} :
#|  -load the coordinates of a trajectory according to the information
#|   _ stored in the trajInfo array .
#|  -arguments :
#|    -id :
#|      -vmd molecule ID of the molecule corresponding to the trajectory .
#|      -the psf file must be already loaded by means of the trajInfo script
#|       _ specific for the trajectory .
#|      -may be specified as "top" to use the top molecule ;
#|    -l_fragId :
#|      -list of traj fragment specifications (decoded by proc trajFragSpec) .
#|      -acceptable values :-see trajFragSpec procedure ;
#|      -default value :-"name" ;;
#|    -args (variable arguments) :
#|      -"loadStepUsr", "loadStep", "loadstep", "loadFreq", "stride" :
#|        -frequency for loading frames from the dcd file .
#|        -acceptable values :
#|          -"trajFrag" :
#|            -loadStep value specified in the trajFrag lists is used .
#|            -this value is defaulted in trajInfo_* script files .
#|            -this value is convenient to process the traj fragments
#|             _ sequentially, i.e. to analyze each fragment one at a time ;
#|          -"trajInfo", "global" :
#|            -uses the value stored in trajInfo(<id>,loadStep) .
#|            -this value is defaulted in trajInfo_* script files .
#|            -this value is convenient to fit the whole trajectory in
#|             _ available RAM ;
#|          -a positive integer (excluding zero), i.e. 1, 2, ... ;
#|        -default value :-"trajInfo" ;
#|        -trajectory size is divided by this number, so for big trajectories a
#|         _ larger value may avoid memory overflow ;
#|      -exclude, except :-fragId exclude list for the trajFragSpec proc .
#|        -default value :-{} ;;
#|      -loSt, channelId, log :-output stream for log messages .
#|        -default value :-stdout ;;;;
#|  -notes :
#|    -each fragment coordinates are to be appended to the molecule .
#|    -it is up to the user to load the fragments in a correct order .
#|    -it is recomended to use trajDelete before loading small or
#|     _ non-contiguous fragments ;;
proc trajLoad {id {l_fragId "name"} args} {
  global trajInfo ind maxTrajSizeGB currTrajSize
  global trajFragList
# default values for optional arguments
  set loSt stdout; set exclude {}; set loadStepUsr "trajInfo"
  if {$id == "top"} {set id [molinfo top]}
  if {$id == -1} {puts $loSt "No trajInfo loaded."; return}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "loadStepUsr" - "loadStep" - "loadstep" - "loadFreq" - "stride" {
            set loadStepUsr $val
            }
          "exclude" - "except" {set exclude $val}
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {puts $loSt "trajLoad: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "trajLoad: Odd number of optional arguments! args: $args"
      return ""
      }
  set l_fragId [trajFragSpec $l_fragId $id exclude $exclude loSt $loSt]
  puts $loSt "\nLoading traj coordinates: id: $id name: $trajInfo($id,name)"
#  puts $loSt "  $trajInfo($id,desc)"
  for {set tfi 0} {$tfi < [llength $trajInfo($id,trajFrag)]} {incr tfi} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
# extract and assign values for new variables from the trajFrag list
    foreach elem $trajFragList {set $elem [lindex $l_trajFrag $ind($elem)]}
    if {($loadStepUsr == "trajInfo") || ($loadStepUsr == "global")} {
      set loadStep $trajInfo($id,loadStep)
    } elseif {$loadStepUsr != "trajFrag"} {
      set loadStep $loadStepUsr
      }
    foreach fragId $l_fragId {
      if {($fragName == $fragId)||($simName == $fragId)||($fragId == "all")} {
        puts $loSt "  trajectory fragment: $fragName (simName: $simName)"
        if {$iniFrame != "unk"} {
          puts $loSt "Coordinates of trajectory fragment already loaded"
          break
          }
# checks whether the trajectory will fit in memory
        set trajSize [expr {$dcdSize/$loadStep}]
        set totalSize [expr {$currTrajSize + $trajSize}]
        if {$totalSize > $maxTrajSizeGB} {
          puts $loSt "Not enough memory to load traj with load step: $loadStep"
# NOTE: code for automatically reducing the load step may be included
          return
          }
        if {$timeStep == 0} {
          set dt $dcdFreq} else {set dt [expr {$dcdFreq*$timeStep/1000000.0}]}
#       puts $loSt "Simulation time in fragment: [expr {$finTime-$iniTime+$dt}]"
# load trajectory fragment
        set iniFrame [expr {[molinfo $id get numframes]}]
        puts -nonewline $loSt "  file: $dcdFile loadStep: $loadStep "
        puts -nonewline $loSt "([expr {$finTime-$iniTime+$dt}] ns) "
        mol addfile $dcdFile step $loadStep waitfor all $id
        set currTrajSize $totalSize
        set finFrame [expr {[molinfo $id get numframes] - 1}]
        set loadedFrames [expr {$finFrame - $iniFrame + 1}]
        set frameTime [expr {$dt*$loadStep}]
# updates initial and final frames
#        set finTime [expr {$iniTime + $dt*$loadStep*$loadedFrames}]
# NOTE: the finTime is not changed (in case the last frame was not loaded)
#   only when the trajFrag is deleted partially this would be changed
#   thus, dt has to be calculated based on the dcdFreq and not from the sim time
        set l_trajFrag [lreplace $l_trajFrag $ind(iniFrame) \
                       $ind(frameTime) $iniFrame $finFrame $frameTime]
        set trajInfo($id,trajFrag) [lreplace $trajInfo($id,trajFrag) $tfi \
                                             $tfi $l_trajFrag]
        puts $loSt "frames: $loadedFrames"
        puts $loSt "  Memory used: $currTrajSize of $maxTrajSizeGB GB"
        }
      }
    }
  }

