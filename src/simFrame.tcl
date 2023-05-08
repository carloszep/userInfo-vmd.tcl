#|-proc simFrame {nsTime {id "top"} {fragId "name"}} :
#|  -returns the frame number corresponding to a time in ns for a simulation .
#|  -the result will be transformed to an integer .
#|  -if the result exeed the final frame of the specified traj frag the last
#|   _ frame of that fragment is returned (error for rounding) :
#|    -the behavior must be tested more extensively ;
#|  -arguments :
#|    -nsTime :-time in ns corresponding the the frame requested ;
#|    -id :-Id of the vmd trajectory included in trajInfo array .
#|      -the value "top" may be specified to use the top molecule ;
#|    -fragId :-specifies the traj fragments to be considered .
#|      -acceptable values :
#|        -"name" :
#|          -this is the default value .
#|          -the name of the simulation is used (from trajInfo(id,name)) .
#|          -this makes the program consider the last simName which
#|           _ generally corresponds to the production simulation ;
#|        -a fragName or simName included in the tragFrag lists ;;;
#|  -notes :
#|    -the requested simTime may be included in several fragments, but only
#|     _ the last occurrence in the trajFragments is considered .
#|    -trajFragSpec may be used, but this may be affect performance ;;
proc simFrame {nsTime {id "top"} {fragId "name"}} {
  global trajInfo ind
  if {$id == "top"} {set id [molinfo top]}
  set numFrames [molinfo $id get numframes]
  if {$numFrames == 0} {return -1}
  if {$fragId == "name"} {set fragId $trajInfo($id,name)}
#  set fragId [trajFragSpec $fragId $id]
  set numFrag [llength $trajInfo($id,trajFrag)]
#  set finTimeLast [lindex $trajInfo($id,trajFrag) end $ind(finTime)]
#  if {$finTimeLast > $nsTime} {
#    return -1
#    set nsTime $finTimeLast
#    }
# searches for the trajFrag where the specified simulation time is included
  for {set tfi [expr {$numFrag - 1}]} {$tfi >= 0} {incr tfi -1} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
    if {($fragId == [lindex $l_trajFrag $ind(fragName)]) || \
        ($fragId == [lindex $l_trajFrag $ind(simName)])} {
      set iniTime [lindex $l_trajFrag $ind(iniTime)]
      set finTime [lindex $l_trajFrag $ind(finTime)]
      if {($nsTime >= $iniTime) && ($nsTime <= $finTime)} {break}
      }
    }
  if {$tfi == -1} {return -1}
  set iniFrame [lindex $l_trajFrag $ind(iniFrame)]
  if {$iniFrame == "unk"} {return -1}
  set frameTime [lindex $l_trajFrag $ind(frameTime)]
  set frame [expr {round($iniFrame + ($nsTime - $iniTime)/$frameTime)}]
  if {$frame > [getTrajProp "finFrame" $id $fragId]} {
    set frame [getTrajProp "finFrame" $id $fragId]}   ;# to avoid rounding error
# note : 
#  return [expr {int($iniFrame + ($nsTime - $iniTime)/$frameTime)}]
  return $frame
  }

