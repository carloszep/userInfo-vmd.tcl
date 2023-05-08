#|-proc trajDelete {id {l_fragId "all"} args} :
#|  -removes traj coordinates from memory for the specified traj fragments .
#|  -arguments :
#|    -id :-VMD molecule ID that is refered in the trajInfo array ;
#|    -l_fragId :-list of traj fragments to load .
#|      -trajectory fragment specification(s) (decoded by proc trajFragSpec) .
#|      -acceptable values :-see trajFragSpec procedure ;
#|      -default value :-"all" ;;
#|    -args (variable arguments) :
#|      -loSt, channelId, log :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -exclude, except :-fragId exclude list for the trajFragSpec proc .
#|        -default value :-{} ;;;;
#|  -notes :
#|    -the trajFrag list is readed in reverse order ;;
proc trajDelete {id {l_fragId "all"} args} {
  global trajInfo ind maxTrajSizeGB currTrajSize
  global trajFragList
# default values for arguments
  set loSt stdout; set exclude {}
  if {$id == "top"} {set id [molinfo top]}
  if {$id == -1} {puts $loSt "trajDelete: No trajInfo loaded."; return}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "loSt" -  "channelId" - "log" {set loSt $val}
          "exclude" - "except" {set exclude $val}
          default {puts $loSt "trajDelete: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "trajDelete: Odd number of optional arguments! args: $args"
      return ""
      }
  set l_fragId [trajFragSpec $l_fragId $id exclude $exclude loSt $loSt]
  if {$l_fragId == "name"} {set l_fragId $trajInfo($id,name)}
  puts $loSt "\nDeleting traj coordinates: Id $id name: $trajInfo($id,name)"
  for {set tfi [expr {[llength $trajInfo($id,trajFrag)] - 1}]} {$tfi >= 0} \
                                                               {incr tfi -1} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
# extract and assign values for new variables from the trajFrag list
    foreach elem $trajFragList {set $elem [lindex $l_trajFrag $ind($elem)]}
    foreach fragId $l_fragId {
      if {($fragName == $fragId)||($simName == $fragId)||($fragId == "all")} {
        puts $loSt "  trajectory fragment: $fragName (simName: $simName)"
        if {$iniFrame == "unk"} {
          puts $loSt "Coordinates of traj fragment $fragName not loaded"
          continue
        }
        animate delete beg $iniFrame end $finFrame $id
# calcualtes the amount of memory released
        if {$timeStep == 0} {
#          set loadStepUsed [expr {round($frameTime/$dcdFreq)}]
           set loadStepUsed 1
        } else {
          set loadStepUsed \
            [expr {round($frameTime*1000000.0/$dcdFreq/$timeStep)}]}
        set trajSize [expr {$dcdSize/$loadStepUsed}]
        set currTrajSize [expr {$currTrajSize - $trajSize}]
        set framesDeleted [expr {$finFrame - $iniFrame + 1}]
        puts $loSt "  Deleted $framesDeleted frames."
        if {[expr {abs($currTrajSize) < 0.00000001}]} {
          set currTrajSize 0.0
          }
        puts $loSt \
        "  Memory used: $currTrajSize of $maxTrajSizeGB GB"
# deletes the previous frames values and updates the trajFrag list
        set l_trajFrag [lreplace $l_trajFrag $ind(iniFrame) $ind(frameTime) \
                                    "unk" "unk" "unk"]
        set trajInfo($id,trajFrag) [lreplace $trajInfo($id,trajFrag) $tfi \
                                             $tfi $l_trajFrag]
# shift the frame numbers of the fragments after the current just deleted
        for {set tfr 0} {$tfr < [llength $trajInfo($id,trajFrag)]} {incr tfr} {
          set l_trajFragR [lindex $trajInfo($id,trajFrag) $tfr]
          set iniFrameR [lindex $l_trajFragR $ind(iniFrame)]
          set finFrameR [lindex $l_trajFragR $ind(finFrame)]
          if {($iniFrameR != "unk")} {
            if {($iniFrameR > $finFrame)} {
              set l_trajFragR [lreplace $l_trajFragR \
                                        $ind(iniFrame) $ind(finFrame) \
                                        [expr {$iniFrameR - $framesDeleted}] \
                                        [expr {$finFrameR - $framesDeleted}]]
              set trajInfo($id,trajFrag) [lreplace $trajInfo($id,trajFrag) \
                                                   $tfr $tfr $l_trajFragR]
              }
            }
          }
        }
      }
    }
  }   ;# trajDelete

