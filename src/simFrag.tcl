#|-proc simFrag {nsTime {id "top"} {fragId "name"}} :
#|  -returns the fragName of the fragment containing a specific simulation
#|   _ time (in ns) considering a list of fragIds specifications .
#|  -returns -1 id the time is not included in the fragId specification .
#|  -arguments :
#|    -nsTime :-simulation time in nano secounds to search for .
#|      -acceptable values :
#|      -"now" :-the current frame is used to determine the sim time ;;;
#|    -id :-mol Id of the trajectory considered in the trajInfo array ;
#|    -fragId :-list of fragment specification to be considered in the search :
#|      -value "name" (default) may be used to consider the sim name in
#|       _ trajInfo(id,name) value ;;
#|   -NOTE :
#|     -it seems neccesary to preprocess fragId with trajFragSpec
#|      _ before calling this proc.) .
#|     -has to be considered to include trajFragSpec within this proc ;;;
proc simFrag {nsTime {id "top"} {fragId "name"}} {
  global trajInfo ind
  if {$id == "top"} {set id [molinfo top]}
  if {$fragId == "name"} {set fragId $trajInfo($id,name)}
  if {$nsTime == "now"} {set nsTime [simTime [molinfo $id get frame] $id]}
  set fragName "-1"
  set numFrag [llength $trajInfo($id,trajFrag)]
  # searches for the trajFrag where the specified simulation time is included
  foreach frag $fragId {
    for {set tfi [expr {$numFrag - 1}]} {$tfi >= 0} {incr tfi -1} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
      if {($frag == [lindex $l_trajFrag $ind(fragName)]) || \
          ($frag == [lindex $l_trajFrag $ind(simName)])} {
        set iniTime [lindex $l_trajFrag $ind(iniTime)]
        set finTime [lindex $l_trajFrag $ind(finTime)]
        if {($nsTime >= $iniTime) && ($nsTime <= $finTime)} {
          set fragName [lindex $l_trajFrag $ind(fragName)]
          break
          }
        }
      }
    if {$fragName != "-1"} {break}
    }
  return $fragName
  }   ;# simFrag

