#|-proc getSelInfo {selId selN idN fragN frmN frstN lstN stpN updN ldstpN
#|                                                            _ lblN descN} :
#|  -extract the selection information associated to a selId and assign the
#|   _ values to the variables selN, idN, fragN frmN, frstN, lstN, stpN,
#|   _ updN, ldstpN, and lblN .
#|  -it will assign default values if the keys have not been specified .
#|  -if the molId specified is already loaded the default values will be set
#|   _ accordingly .
#|  -arguments :
#|    -selId :-selId used as key in the selInfo array ;
#|    -names of the variables used in the calling procedure :
#|      -selN :-selection text variable name .
#|        -default value to return :-'none' ;;
#|      -idN :-VMD mol ID .
#|        -if "top" is specified the top molecule Id will be returned if
#|         _ exists .
#|        -default value to return :
#|          -the current top molecule .
#|          -'top' if no molecule was loaded into VMD ;;
#|      -fragN :-name of the variable with the fragId specification .
#|        -this fragId may be used to determine first and last frames if they
#|         _ are not specified in the selInfo array .
#|        -default value to return :-"name" ;;
#|      -frmN :-frame variable name .
#|        -if "now" is specified the curent frame will be returned if the
#|         _ trajectory is already loaded  .
#|        -default value to return :
#|          -the current frame if the trajectory is already loaded .
#|          -'now' if the trajectory has not been loaded into VMD .
#|          -the first frame of the specified fragId if the current frame
#|           _ is outsde the fragment's range ;;
#|      -frstN :-first frame variable name .
#|        -default value to return :-the value obtained from the trajFrag
#|         _specification, or "first" otherwise ;;
#|      -lstN :-last frame variable name .
#|        -if "last" is specified the last frame number is returned if the
#|         _ trajectory was already loaded .
#|        -default value to return :
#|          -the value obtained from the trajFrag lists according to the
#|           _ fragId specification .
#|          -the last frame of the specified molecule .
#|          -'last' if the trajectory has not been loaded into VMD ;;
#|      -stpN :-step (stride) variable name .
#|        -default value to return :-1 ;
#|      -updN :-update selection variable name .
#|        -default value to return :-0 ;;
#|      -ldstpN :-loadStep variable name .
#|        -default value to return :-'trajInfo' ;;
#|      -lblN :-label variable name .
#|        -default value to return :-the specified 'selId' ;;
#|      -descN :-varible name for the description text of the selId .
#|        -default value to return :-"" ;;;;
#|  -notes :
#|    -the default values will be more accurated if the trajectories are loaded
#|     _ before calling this proc .
#|    -new selInfo array key 'desc' is not managed by this proc yet ;;
proc getSelInfo {selId selN idN fragN frmN frstN lstN stpN updN ldstpN \
                                                                lblN descN} {
  global selInfo trajInfo
  upvar $selN sel
  upvar $idN id
  upvar $fragN frag
  upvar $frmN frm
  upvar $frstN frst
  upvar $lstN lst
  upvar $stpN stp
  upvar $updN upd
  upvar $ldstpN ldstp
  upvar $lblN lbl
  upvar $descN desc
  if {[info exists selInfo($selId,selTxt)]} {
    set sel $selInfo($selId,selTxt)} else {set sel "none"}
  if {[info exists selInfo($selId,molId)]} {
      if {$selInfo($selId,molId) == "top"} {
          if {[molinfo num] > 0} {set id [molinfo top]} else {set id "top"}
        } else {set id $selInfo($selId,molId)}
    } else {if {[molinfo num] > 0} {set id [molinfo top]} else {set id "top"}}
  if {[info exists selInfo($selId,fragId)]} {
      if {$selInfo($selId,fragId) == "name"} {
          if {$id != "top"} {
            set frag $trajInfo($id,name)} else {set frag "name"}
        } else {set frag $selInfo($selId,fragId)}
    } else {
      if {$id == "top"} {set frag "name"} else {set frag $trajInfo($id,name)}}
  if {[info exists selInfo($selId,frame)]} {
      if {$selInfo($selId,frame) == "now"} {
          if {[info exists trajInfo($id,name)]} {
              if {[molinfo $id get numframes] != 0} {
                  set frm [simTime [molinfo $id get frame] $id]
                } else {set frm "now"}
            } else {set frm "now"}
        } else {set frm $selInfo($selId,frame)}
    } else {
      if {[info exists trajInfo($id,name)]} {
          if {[getTrajProp "loaded" $id $frag]} {
            set frm [simTime [molinfo $id get frame] $id]
            if {[simFrag $frm $id [trajFragSpec $frag $id]] == -1} {
              set frm [getTrajProp "iniTime" $id $frag]}
          } else {set frm "now"}
        } else {set frm "now"}}
  if {[info exists selInfo($selId,first)]} {
      set frst $selInfo($selId,first)
      if {($frst == "first") || ($frst == "iniTime")} {
        set frst [getTrajProp "iniTime" $id $frag]
        if {($frst == "unk")||($frst == "")} {set frst "first"}
        }
    } else {
      if {$id == "top"} {set frst "first"} else {
        set frst [getTrajProp "iniTime" $id $frag]
        if {($frst == "unk")||($frst == "")} {set frst "first"}
        }
      }
  if {[info exists selInfo($selId,last)]} {
      if {$selInfo($selId,last) == "last"} {
          if {[info exists trajInfo($id,name)]} {
              set lst [getTrajProp "finTime" $id $frag]
                if {($lst == "unk") || ($lst == "")} {set lst "last"}
            } else {set lst "last"}
        } else {set lst $selInfo($selId,last)}
    } else {
      if {[info exists trajInfo($id,name)]} {
          set lst [getTrajProp "finTime" $id $frag]
          if {($lst == "unk") || ($lst == "")} {set lst "last"}
        } else {set lst "last"}
      }
  if {[info exists selInfo($selId,step)]} {
    set stp $selInfo($selId,step)} else {set stp 1}
  if {[info exists selInfo($selId,updateSel)]} {
    set upd $selInfo($selId,updateSel)} else {set upd 0}
  if {[info exists selInfo($selId,loadStep)]} {
    set ldstp $selInfo($selId,loadStep)} else {set ldstp "trajInfo"}
  if {[info exists selInfo($selId,label)]} {
    set lbl $selInfo($selId,label)} else {set lbl $selId}
  if {[info exists selInfo($selId,desc]} {
    set desc $selInfo($selId,desc)} else {set desc ""}
  }

