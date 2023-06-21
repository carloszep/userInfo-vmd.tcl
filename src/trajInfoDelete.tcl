#|-proc trajInfoDelete {id {loSt stdout}} :
#|  -unset all the information associated to a molId in the trajInfo array  .
#|  -the structure and trajectory will also be deleted from VMD .
#|  -arguments :
#|    -id :-VMD mol ID which is refered in the trajInfo array .
#|      -may specified as "top" to use the top molecule ;
#|    -loSt :-output stream for log messages ;;;
proc trajInfoDelete {id {loSt stdout}} {
  global trajInfo selInfo
  if {$id == "top"} {set id [molinfo top]}
  puts $loSt \
    "\nDeleting info, struct and traj of $trajInfo($id,name) (Id: $id)"
  trajDelete $id "all" loSt $loSt
  puts $loSt "Deleting default selections associated to this molecule..."
  puts $loSt "  Original selInfo array keys:"
  puts $loSt "    [array names selInfo]"
  puts $loSt "  Array keys to be deleted:"
  puts $loSt "    [array names selInfo *-m$id,*]"
  array unset selInfo *-m$id,*
  puts $loSt "Deleting trajInfo information associated to this molecule..."
  array unset trajInfo $id,*
  mol delete $id
  puts $loSt "Simulation information, structure and trajectory deleted."
  }

