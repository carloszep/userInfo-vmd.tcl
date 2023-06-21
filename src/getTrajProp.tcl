#|-proc getTrajProp {prop {id "top"} {l_fragId "name"} args} :
#|  -return the value of a specified property stored in the trajFrag lists .
#|  -the returned values may vary depending on the current status of the
#|   _ trajectory, e.g. whether or not the coordinates were loaded .
#|  -arguments :
#|    -prop :-property whose value is returned for the specified trajFrag
#|     _ list(s) .
#|      -acceptable values :
#|        -any element of the trajFragList .
#|        -"fragName", "dcdFile" :
#|          -return a list of the respective data for the number of fragments
#|           _ refered through the l_fragId argument ;
#|        -"simName", "timeStep", "dcdFreq", "loadStep", "frameTime" :
#|          -return a list of values of the property if these are heterogeneous
#|           _ among the matching trajFrag lists refered by the l_fragId arg .
#|          -return a single value of the property if the value is the same for
#|           _ the trajFrag lists refered by the l_fragId argument :
#|            -for these properties the single-structure trajFrag lists are
#|             _ excluded (e.g. same iniTime and finTime) .
#|            -thus, single-structure coordinates may be considered as a part
#|             _ of the subsequent traj fragment (i.e. initial structure) .
#|            -useful for a proc to know whether a set of fragments are part of
#|             _ a single sim set-up or they have heterogeneous sim/min runs ;;
#|        -"dcdSize" :
#|          -return the sum of the dcd sizes considered in the trajFrag lists
#|           _ matching the l_fragId specification ;
#|        -"iniTime" :
#|          -return the iniTime property value from the first trajFrag list
#|           _ matching the specification by the l_fragId argument ;
#|        -"finTime" :
#|          -returns the finTime property value from the last trajFrag list
#|           _ matching the specification by the l_fragId argument ;
#|        -"iniFrame" :
#|          -returns the iniFrame property value from the first loaded trajFrag
#|           _ list matching the l_fragId specification .
#|          -returns "unk" if all the fragFrag lists specifyied by l_fragId
#|           _ were not loaded ;
#|        -"finFrame" :
#|          -returns the finFrame property value from the last loaded trajFrag
#|           _ list matching the l_fragId specification .
#|          -returns "unk" if none of the traj fragments specifyed by l_fragId
#|           _ were no loaded ;
#|        -"numTrajFrag", "numFrag" :
#|          -returns the number of trajFrag lists considered in
#|           _ the l_fragId specification ;
#|        -"listTrajFrag" :
#|          -returns a list of full trajFrag lists containing fragments
#|           _ matching the specification by the l_fragId argument .
#|          -may be used to extract sublists of the trajFrag stored in the
#|           _ trajInfo array ;
#|        -"seqTraj", "seqTime" :
#|          -meaning "are the traj fragments sequential in simulation time?" .
#|          -returns "1" if the list of trajFrags specified by l_fragId have
#|           _ a sequential (growing) trend in time; returns "0" otherwise :
#|            -i.e.: an equilibration sequence would return value of 0 ;;
#|        -"seqFrame", "seqFrames" :
#|          -meaning "are the fragment's frames sequential or consecutive?" .
#|          -returns "1" if the coordinates of the trajFrag lists specified
#|           _ in l_fragId were loaded consecutively and in order .
#|          -returns "0" when the fin and iniFrames between contiguous
#|           _ trajFrags specified by l_trajFrag are not consecutive ;
#|        -"loaded" :
#|          -returns "1" if all trajFrags specified by l_fragId were loaded,
#|           _ returns "0" otherwise ;
#|        -"loadedFrag", "fragLoaded" :
#|          -returns a list of fragNames that are loaded in the fragID spec ;;;
#|    -id :
#|      -acceptable values :
#|        -a VMD molid of a loaded molecule considered in the trajInfo array .
#|        -"top" :-the current top molecule in VMD is used ;;
#|      -default value :-"top" ;;
#|    -l_fragId :
#|      -list of traj fragment specifications (decoded by proc trajFragSpec) .
#|      -acceptable values :-see trajFragSpec procedure ;
#|      -default value :-"name" ;;
#|    -args (variable arguments) :
#|      -loSt, channelId, log :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -exlude, except, excl, exFragId :
#|        -list of fragIds to be excluded from the fragment specification .
#|        -acceptable values :
#|          -the same as l_fragId except for "all", "exclude", and
#|           _ "except" keywords ;;;;
#|  -notes :
#|    -the proc retuns strictly the values stored in the trajFrag lists, so
#|     _ the overrided default properties (i.e. loadStep) may have different
#|     _ values from those actually used .
#|    -usage examples... :- ;;;
proc getTrajProp {prop {id "top"} {l_fragId "name"} args} {
  global ind trajInfo trajFragList
# default values for arguments
  set loSt stdout; set exclude {}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "loSt" - "channelId" - "log" {set loSt $val}
          "exclude" - "except" - "excl" - "exFragId" {set exclude $val}
          default {puts $loSt "trajFragSpec: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "trajFragSpec: Odd number of optional arguments! args: $args"
      return ""
      }
  if {$id == "top"} {set id [molinfo top]}
  if {$id == -1} {puts $loSt "trajDelete: No trajInfo loaded."; return}
  set l_fragId [trajFragSpec $l_fragId $id exclude $exclude loSt $loSt]
# temporal variables to store the value(s) of the property requested
  set l_valProp {}
  set sumValProp 0.0
  set homoValProp ""
  set retValProp ""
  for {set tfi 0} {$tfi < [llength $trajInfo($id,trajFrag)]} {incr tfi} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
    foreach fragId $l_fragId {
      if {$fragId == [lindex $l_trajFrag $ind(fragName)]} {
# the current trajFrag is included in the l_fragId specification
        switch $prop {
          "fragName" - "dcdFile" {
            lappend l_valProp [lindex $l_trajFrag $ind($prop)]
            set retValProp $l_valProp
            }
          "simName" - "timeStep" - "dcdFreq" - "loadStep" - "frameTime" {
            set valProp [lindex $l_trajFrag $ind($prop)]
            lappend l_valProp $valProp
# exclude the single-structure trajFrags as homoValProp
            if {[lindex $l_trajFrag $ind(iniTime)] != \
                [lindex $l_trajFrag $ind(finTime)]} {
              if {$homoValProp == ""} {
                set homoValProp $valProp
              } else {
                if {$valProp == $homoValProp} {
                  set homoValProp $valProp
                } else {
                  set homoValProp "none"
                  }
                }
              }
# evaluates whether the prop values are homogeneous up to this trajFrag
            if {$homoValProp == "none"} {
              set retValProp $l_valProp
            } else {
              set retValProp $valProp
              }
            }
          "dcdSize" {
            set valProp [lindex $l_trajFrag $ind($prop)]
            set sumValProp [expr {$sumValProp + $valProp}]
            set retValProp $sumValProp
            }
          "iniTime" {return [lindex $l_trajFrag $ind($prop)]}
          "finTime" {set retValProp [lindex $l_trajFrag $ind($prop)]}
          "iniFrame" {
            set valProp [lindex $l_trajFrag $ind($prop)]
            set retValProp $valProp
            if {$valProp != "unk"} {
              return $valProp
              }
            }
          "finFrame" {
            set valProp [lindex $l_trajFrag $ind($prop)]
            if {$retValProp == ""} {
              set retValProp $valProp
            } else {
              if {$valProp != "unk"} {
                set retValProp $valProp
                }
              }
            }
          "numTrajFrag" - "numFrag" {
            if {$retValProp == ""} {set retValProp 1} else {incr retValProp}
            }
          "listTrajFrag" {lappend retValProp $l_trajFrag}
          "seqTraj" - "seqTime" {
            if {$retValProp == ""} {
              set prevFinTime [lindex $l_trajFrag $ind(finTime)]
              if {$prevFinTime == "unk"} {return 0}
              set retValProp 1
            } else {
              if {$prevFinTime >= [lindex $l_trajFrag $ind(iniTime)]} {
                return 0
              } else {
                set prevFinTime [lindex $l_trajFrag $ind(finTime)]
                }
              }
            }
          "seqFrame" - "seqFrames" {
            if {$retValProp == ""} {
              set prevFinFrame [lindex $l_trajFrag $ind(finFrame)]
              if {$prevFinFrame == "unk"} {return 0}
              set retValProp 1
            } else {
              set currIniFrame [lindex $l_trajFrag $ind(iniFrame)]
              if {$currIniFrame == "unk"} {return 0} 
              if {$prevFinFrame != [expr {$currIniFrame - 1}]} {
                return 0
              } else {
                set prevFinFrame [lindex $l_trajFrag $ind(finFrame)]
                }
              }
            }
          "loaded" {
            if {[lindex $l_trajFrag $ind(iniFrame)] == "unk"} {
              return 0} else {set retValProp 1}
            }
          "loadedFrag" - "fragLoaded" {
            if {[lindex $l_trajFrag $ind(iniFrame)] != "unk"} {
              lappend retValProp $fragId
              }
            }
          default {puts "Unknown property! ($prop)"; return ""}
          }
        break
        }
      }
    }
  return $retValProp
  }   ;# getTrajProp

