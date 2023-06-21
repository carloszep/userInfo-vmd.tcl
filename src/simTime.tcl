#|-proc simTime {frm {id "top"}} :
#|  -Returns the corresponding time (in ns) of a frame in a trajectory .
#|  -if the trajectory is was not previously loaded returns -1 .
#|  -note that it calculates the trajectory length and the number of frames
#|   _ each time the procedure is called which could not be the ideal .
#|  -arguments :
#|    -frm :-frame number .
#|      -acceptable values :
#|        - "now" :-the current frame is used to determine the sim time ;;;
#|    -id :-Id if the vmd trajectory included in trajInfo array .
#|      -the value "top" may be specified to use the top molecule ;;;
proc simTime {frm {id "top"}} {
  global trajInfo ind
  if {$id == "top"} {set id [molinfo top]}
  if {$frm == "now"} {set frm [molinfo $id get frame]}
  set numFrames [molinfo $id get numframes]
  if {$numFrames == 0} {return -1}
  if {$frm >= $numFrames} {
    return -1
#    set frm [expr {[molinfo $id get numframes] - 1}]
    }
# searches the l_trajFrag that includes the indicated frame
  foreach l_trajFrag $trajInfo($id,trajFrag) {
    set iniFrame [lindex $l_trajFrag $ind(iniFrame)]
    set finFrame [lindex $l_trajFrag $ind(finFrame)]
    if {$iniFrame == "unk"} {continue}
    if {($frm >= $iniFrame) && ($frm <= $finFrame)} {break}
    }
# calculates the simulation time
  set iniTime [lindex $l_trajFrag $ind(iniTime)]
  set frameTime [lindex $l_trajFrag $ind(frameTime)]
  return [expr {$iniTime + $frameTime*($frm - $iniFrame)}]
  }

