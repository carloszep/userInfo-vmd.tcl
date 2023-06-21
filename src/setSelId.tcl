#|-proc setSelId {selIds args} :
#|  -declares an entry for the selInfo array to specify an atom selection
#|   _ information from a trajectory .
#|  -arguments :
#|    -selIds :-list of identificator names used to access or point to
#|     _ paricular elements of the selInfo array (the first index) .
#|      -if the specified selId is already used in the selInfo array the
#|       _ specified elements of the array will be updated (overwritten) .
#|      -acceptable values and keywords :
#|        -list of selId names either new or already used in the selInfo array .
#|        -"unset", "delete" :
#|          -unsets all the array keys including the specified selId .
#|          -a list od selIds is acceptable .
#|          -the second argument is then expected to be the selId(s) to delete .
#|          -usage for this case :-'setSelId unset <l_selId>' ;;;;
#|    -args (variable arguments) :
#|      -"selTxt", "sel", "selection", "atomselect", "atomSel" :
#|        -specify the "selTxt" key in the second index of the selInfo array :
#|          -i.e. selInfo(<selId>,selTxt) .
#|          -specify a selection text for the "atomselect" proc of VMD ;
#|        -acceptable values :
#|          -VMD's atom selection texts (see VMD's user guide) .
#|          -a list of values is acceptable corresponding to each selId
#|           _ specified :
#|            -for a single selection text it may be necessary to enclose in
#|             _ curly braces, i.e. {"protein and name CA"} ;
#|          -a single value common to several selIds specified ;
#|        -default :-if no selTxt is specified the array key is not declared ;;
#|      -"molId", "molid", "id", "ID", "mol", "idMol", "vmdMol", "VMDmol" :
#|        -VMD mol ID, preferently already loaded from a trajInfo script .
#|        -acceptable values :
#|          -an integer number corresponding to a VMD mol ID .
#|          -a list of Ids corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"top" :-refers to VMD's top molecule ;
#|          -"unset", "delete" :-the array element ("<selId>,molId") is unset ;;
#|        -default :-the array entry is not declared ;;
#|      -"fragId", "frag", "trajFrag", "fragSpec", "simFrag", "fragment",
#|        _ "fragments" :
#|        -list of trajectory fragment specifications .
#|        -acceptable values :
#|          -a valid fragId specification (see "trajFragSpec" procedure) .
#|          -"unset", "delete" :-the selInfo(<selId>,fragId) key is unset ;
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified ;
#|        -default :-the variable (selInfo array entry) will not be set ;;
#|      -"exclude", "except", "excl", "exFragId" :
#|        -fragId specification "exclude" list .
#|        -this exclude specification will be added to fragId spec if it does
#|         _ not includes already any "exclude" keywords, otherwise is ignored .
#|        -acceptable values :
#|          -see "trajFragSpec" procedure .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to all the specified selIds ;
#|        -default value :-{} ;
#|        -note :-this feature is not functional yet ;;
#|      -"first", "ini", "start", "begin", "iniRange" :
#|        -left-bound (starting) range of the simulation time interval to be
#|         _ considered in the selection (to be used by some analysis procs) .
#|        -specifies the selInfo(<selId>,first) array element .
#|        -acceptable values :
#|          -a real/integer number corresponding to a sim time or frame .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"now" :-the current frame is taken as "first" value ;
#|          -"unset", "delete" :
#|            -the selInfo(<selId>,first) array key is unset ;;
#|        -default :-the selInfo(<selId>,first) key will not be set ;;;
#|      -"last", "end", "final", "stop", "until", "finRange" :
#|        -right-bound (ending) range of the simulation time interval to be
#|         _ considered in the selection (to be used by some analysis procs) .
#|        -specifies the selInfo(<selId>,last) array element .
#|        -acceptable values :
#|          -a real/integer number corresponding to a sim time or frame .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"now" :-the current frame is taken as "last" value ;
#|          -"unset", "delete" :
#|            -the selInfo(<selId>,last) array key is unset ;;
#|        -default :-the selInfo(<selId>,last) key will not be set ;;
#|      -"frame", "time", "simTime", "current", "currFrame", "currTime" :
#|        -frame/time specification (used for some analysis procedures) .
#|        -it is expected to refer to a simulation time in ns .
#|        -specify the 'selInfo(<selId>,frame)' array element .
#|        -acceptable values :
#|          -real (float) value in ns or an integer for a sim frame .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"now" :-refers to the time equivalent to the current frame ;
#|          -"first", "iniTime" :-points to the "iniTime" property of the
#|           _ fragId specification (see "getTrajProp" procedure) ;
#|          -"last", "finTime" :-point to the "finTime" property of the fragId
#|           _ specification (see "getTrajProp" procedure) ;
#|          -"unset", "delete" :
#|            -the selInfo array key ("<selId>,frame") should be unset ;;
#|        -default :-the variable (array entry) will not be set ;;
#|      -"step", "stride", "frameStep", "iterateStep", "sampling" :
#|        -step to be used in an interation loop over the loaded traj frames .
#|        -will affect the duration time and results of some analysis procs .
#|        -specify the 'selInfo(<selId>,step)' array element .
#|        -acceptable values :
#|          -an integer number .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"unset", "delete" :-the selInfo(<selId>,step) key will be unset ;;
#|        -default :-the 'selInfo(<selId>,step)' key will not be set ;;
#|      -"updateSel", "update", "updSel", "selUpdate" :
#|        -flag specifying whether the selection has to be updated every frame .
#|        -specify the 'selInfo(<selId>,updateSel)' array element .
#|        -acceptable values :
#|          -"0" :-the selection will not be updated ;
#|          -"1" :-the selection will be updated every frame ;
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified ;
#|        -default :-the 'selInfo(<selId>,updateSel)' key will not be set ;;
#|      -"loadStep", "dcdLoadStep", "trajLoadStep" :
#|        -trajectory (dcd) load step to be used for analysis procs .
#|        -specifies the 'selInfo(<selId>,loadStep)' array element .
#|        -acceptable values :
#|          -a positive integer value .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"unset", "delete" :
#|            -the selInfo(<selId>,loadStep) key will be unset ;;
#|        -default :-the 'selInfo(<selId>,loadStep)' key will not be set ;;
#|      -"desc", "description" :
#|        -description text for the selection to be used for analysis procs .
#|        -specified the 'selInfo(<selId>,desc)' array element .
#|        -acceptable values :
#|          -text string .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -for a single selId the description text may be eclosed in curly
#|           _ braces besides quotation marks, i.e. {"description text"} .
#|          -"unset", "delete" :-the selInfo(<selId>,desc) key will be unset ;;
#|        -default :-the 'selInfo(<selId>,desc)' key will not be set ;;
#|      -"label", "selIdLabel", "selInfoLabel", "lbl", "legend" :
#|        -short string to be used as data set label for graphs .
#|        -specifies the 'selInfo(<selId>,label)' array element .
#|        -acceptable values :
#|          -short string (name) (no blank spaces preferently) .
#|          -a list of values for each selId specified is acceptable .
#|          -a single value common to several selIds specified .
#|          -"unset", "delete" :-the selInfo(<selId>,label) key will be unset ;;
#|        -default :-the 'selInfo(<selId>,label)' key will not be set ;;
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -any other arbitrary token to be used as new user-defined array key .
#|        -specifies array element 'selInfo(<selId>,<newDefKey>)' .
#|        -must be followed by the respective array key value .
#|        -acceptable values :
#|          -general tokens .
#|          -a list of values corresponding to each selId specified .
#|          -a single value common to several selIds specified .
#|          -"unset", "delete" :
#|            -the selInfo(<selId>,<newDefKey>) key will be unset ;;;;;
#|  -notes :
#|    -this proc is intended to ease selId setting up to the user .
#|    -"getSelInfo" proc will still be used to manage default values of selIds
#|     _ by other procedures .
#|    -the exclude lists as arguments independent from the fragIds lists were
#|     _ not implemented, use exclude lists within the fragId specification ;;
proc setSelId {selIds args} {
  global selInfo
# default values for variables
  set loSt stdout; set excl {}; set procName [lindex [info level 0] 0]
  set sel {}; set id {}; set frag {}; set excl {}; set frst {}; set lst {}
  set frm {}; set stp {}; set upd {}; set ldStp {}; set lbl {}; set desc {}
  set delSelId {}; set def {}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "selTxt" - "sel" - "selection" - "atomselect" - "atomSel" {
            set sel $val}
          "molId" - "molid" - "Id" - "id" - "ID" - "mol" - "idMol" - \
            "vmdMol" - "VMDmol" {set id $val}
          "fragId" - "frag" - "trajFrag" - "fragSpec" - "simFrag" - \
            "fragment" - "fragments" {set frag $val}
          "first" - "ini" - "start" - "begin" - "iniRange" {set frst $val}
          "last" - "end" - "final" - "stop" - "until" - "finRange" {
            set lst $val}
          "frame" - "time" - "simTime" - "current" - "currFrame" - "currTime" {
            set frm $val}
          "step" - "stride" - "frameStep" - "iterateStep" - "sampling" {
            set stp $val}
          "updateSel" - "update" - "updSel" - "selUpdate" {set upd $val}
          "loadStep" - "dcdLoadStep" - "trajLoadStep" {set ldStp $val}
          "desc" - "description" {set desc $val}
          "label" - "selIdLabel" - "selInfoLabel" - "lbl" - "legend" {
            set lbl $val}
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {
            lappend def $arg
            set $arg $val
            }
          }
        }
      }
    } else {   ;# odd number of arguments
# check for 'unset <selId>' case. selId to delete in the first pos in args
      if {($selIds == "unset") || ($selIds == "delete")} {
        set delSelId [lindex $args 0]
        set args [lrange $args 1 end]
      } else {
        puts $loSt "$procName: Odd number of variable arguments! args: $args"
        return ""}
# allows selecting the loSt when using keyword unset/delete
      foreach {arg val} $args {
        switch $arg {
          "loSt" -  "channelId" - "log" {set loSt $val}
          default {puts $loSt "$procName: argument out of context: $arg"}
          }
        }
      }
# check for 'unset <selId>' case, and return if requested
  if {$delSelId != {}} {
    foreach del $delSelId {
      array unset selInfo $del,*
      puts $loSt "$procName: deleted selId: $del,*"
      }
    return ""
    }
# show info
  puts $loSt "\n$procName: setting up selection Ids for the selInfo array..."
  puts $loSt "  List of selIds: $selIds"
  puts $loSt "  Variable arguments: $args"
  if {$def != {}} {puts $loSt "  Unknown array keys to be defined: $def"}
# analyze input arguments to create lists with the same length as selIds
  set nSel [llength $selIds]
  set lvar [concat [list sel id frag frst lst frm stp upd ldStp desc lbl] \
                   [set def]]
  foreach var $lvar {
    set val [set $var]
    set nArg [llength $val]
    if {$nArg == 1} {
      for {set i 1} {$i < $nSel} {incr i} {lappend $var [lindex $val 0]}
    } elseif {(($nArg > 1) && ($nArg < $nSel)) || ($nArg > $nSel)} {
      puts $loSt "$procName: Incorrect number of values for argument $var"
      return ""
      }
    }
# correct fragIds if exclude lists were specified (NOT IMPLEMENTED YET)
  if {[llength $excl] > 0} {
    puts $loSt \
      "WARNING: exclude lists independent from fragIds not implemented yet!"
    }
# set the specified selInfo array elements
  for {set s 0} {$s < [llength $selIds]} {incr s} {
    set selId [lindex $selIds $s]
    foreach var $lvar {
      set val [set $var]
      if {[llength $val] > 0} {
        set selVal [lindex $val $s]
        switch $var {
          sel {
            puts $loSt "  setting 'selInfo($selId,selTxt)' to '$selVal'"
            set selInfo($selId,selTxt) $selVal
            }
          id {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,molId
              puts $loSt "  'selInfo($selId,molId)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,molId)' to '$selVal'"
              set selInfo($selId,molId) $selVal
              }
            }
          frag {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,fragId
              puts $loSt "  'selInfo($selId,fragId)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,fragId)' to '$selVal'"
              set selInfo($selId,fragId) $selVal
              }
            }
          frst {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,first
              puts $loSt "  'selInfo($selId,first)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,first)' to '$selVal'"
              set selInfo($selId,first) $selVal
              }
            }
          lst {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,last
              puts $loSt "  'selInfo($selId,last)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,last)' to '$selVal'"
              set selInfo($selId,last) $selVal
              }
            }
          frm {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,frame
              puts $loSt "  'selInfo($selId,frame)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,frame)' to '$selVal'"
              set selInfo($selId,frame) $selVal
              }
            }
          stp {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,step
              puts $loSt "  'selInfo($selId,step)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,step)' to '$selVal'"
              set selInfo($selId,step) $selVal
              }
            }
          upd {
            puts $loSt "  setting 'selInfo($selId,updateSel)' to '$selVal'"
            set selInfo($selId,updateSel) $selVal
            }
          ldStp {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,loadStep
              puts $loSt "  'selInfo($selId,loadStep)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,loadStep)' to '$selVal'"
              set selInfo($selId,loadStep) $selVal
              }
            }
          desc {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,desc
              puts $loSt "  'selInfo($selId,desc)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,desc)' to '$selVal'"
              set selInfo($selId,desc) $selVal
              }
            }
          lbl {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,label
              puts $loSt "  'selInfo($selId,label)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,label)' to '$selVal'"
              set selInfo($selId,label) $selVal
              }
            }
          default {
            if {($selVal == "unset") || ($selVal == "delete")} {
              array unset selInfo $selId,$var
              puts $loSt "  'selInfo($selId,$var)' deleted"
            } else {
              puts $loSt "  setting 'selInfo($selId,$var)' to '$selVal'"
              set selInfo($selId,$var) $selVal
              }
            }
          }
        }
      }
    }
  }

