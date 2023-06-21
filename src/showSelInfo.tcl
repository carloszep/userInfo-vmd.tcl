#|-proc showSelInfo {selId {loSt stdout}} :
#|  -reports atom selection information stored in the selInfo array .
#|  -arguments :
#|    -selId :-array keys considered to show info ;
#|    -loSt :-optional output stream for log messages .-default is stdout ;;
#|  -Notes :-No error checking is included yet ;;
proc showSelInfo {selId {loSt stdout}} {
  global selInfo
  puts $loSt "\nshowSelInfo: Selection information (selId: $selId)"
  if {[info exists selInfo($selId,selTxt)]} {
      puts $loSt "  Selection text: $selInfo($selId,selTxt)"
      set selTxt $selInfo($selId,selTxt)
    } else {
      puts $loSt "  No selection text specified. Default is 'none'."
      set selTxt "none"
      }
  if {[info exists selInfo($selId,molId)]} {
      puts $loSt "  VMD molecule ID: $selInfo($selId,molId)"
      set molId $selInfo($selId,molId)
    } else {
      puts $loSt "  No molId associated. Default is the 'top' molecule."
      set molId "top"
      }
  if {[info exists selInfo($selId,fragId)]} {
      puts $loSt "  Trajectory fragment (fragId) specification:\
                    $selInfo($selId,fragId)"
      set fragId $selInfo($selId,fragId)
    } else {
      puts $loSt "  No trajectory fragment (fragId) specified.\
                    Default is 'name'."
      set fragId "name"
      }
  if {[info exists selInfo($selId,frame)]} {
      puts $loSt "  Selection frame: $selInfo($selId,frame)"
      set frm $selInfo($selId,frame)
    } else {
      puts $loSt "  No frame specified. Default is the current frame or 'now'."
      set frm "now"
      }
  if {[info exists selInfo($selId,first)]} {
      puts $loSt "  Selection range first frame: $selInfo($selId,first)"
    } else {
      puts $loSt "  No first time specified. Default is obtained from\
                    the fragId specification, or 'first' otherwise."}
  if {[info exists selInfo($selId,last)]} {
      puts $loSt "  Selection range last frame: $selInfo($selId,last)"
    } else {
      puts $loSt "  No last frame specified. Default is obtained from\
                    the fragId specification, or 'last' otherwise."}
  if {[info exists selInfo($selId,step)]} {
      puts $loSt "  Selection range step: $selInfo($selId,step)"
    } else {puts $loSt "  No step specified. Default is 1."}
  if {[info exists selInfo($selId,updateSel)]} {
      if {$selInfo($selId,updateSel)} {
          puts $loSt "  The selection will be updated every frame."
        } else {puts $loSt "  The selection will not be updated every frame."}
    } else {puts $loSt "  The selection will not be updated every frame."}
  if {[info exists selInfo($selId,loadStep)]} {
      puts $loSt "  Load step for the fragId: $selInfo($selId,loadStep)"
    } else {puts $loSt "  No load step specified. Default is 'trajFrag'"}
  if {[info exists selInfo($selId,desc)]} {
      puts $loSt "  Description text for the selId: $selInfo($selId,desc)"
    } else {puts $loSt "  No description text specified. Default is ''"}
  if {[info exists selInfo($selId,label)]} {
      puts $loSt "  Label (legend) specified: $selInfo($selId,label)"
    } else {puts $loSt "  No label specified. Default is the selId, '$selId'"}
  set tempSel [atomselect $molId $selTxt frame $frm]
  puts $loSt "The selection comprises [$tempSel num] atoms."
  $tempSel delete
  }

