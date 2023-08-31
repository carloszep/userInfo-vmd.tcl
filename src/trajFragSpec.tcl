#|-trajFragSpec.tcl :
#|  -procedure trajFragSpec implemented for the userInfoLib namespace .

#|  -namespace userInfoLib :
namespace eval userInfoLib {
#|    -export :
#|      -trajFragSpec ;
namespace export trajFragSpec
#|    -namespace commands :

#|-proc trajFragSpec {l_fragId {id "top"} args} :
#|  -interpreter for user-defined trajectory fragment specifications (fragId) .
#|  -returns a valid specification of trajectory fragments for a molecule Id :
#|    -the specification is a sublist of fragNames refering to particular
#|     _ elements of the full trajFrag list stored in the trajInfo array .
#|    -the returned list will be correctly sorted with no repeated elements ;
#|  -returns an empty string if the fragId specification is not valid .
#|  -the fragIds included in the "exclude" list as well as the unrecognized
#|   _ fragIds will be excluded from the returned list .
#|  -arguments :
#|    -l_fragId :-list of user-provided traj fragment specifications .
#|      -acceptable values in the list :
#|        -list of fragNames or simNames :
#|          -each refer to names at positions 0 and 1 of the trajFrag list
#|           _ in the trajInfo array ;
#|        -"all" :-a list with all the fragNames is returned .
#|          -all fragments in the trajFrag list are specified ;
#|        -"name" :-the trajInfo(id,name) value will be considered .
#|          -generally will consider the fragments for the last production sim ;
#|        -"first" :-the frag in position 0 of the trajFrag list is specified ;
#|        -"end", "last" :-the last traj fragment is specified ;
#|        -"end-<n>" :-the fragment at position last minus n is specified ;
#|        -"<fragId> to <fragId>", "<fragId> - <fragId>" :
#|          -range of consecutive elements ;
#|        -"exclude", "except", "excl", "exFragId" :
#|          -keyword to include the exclude list in the same fragId spec .
#|          -fragIds after the keyword will be considered as an exclude list .
#|          -overrides the "exclude" argument .
#|          -in the position 0 of the l_fragId specification (as first fragId)
#|           _ will exclude subsequent fragIds from the "all" specification ;
#|        -an integer refering to the fragment position in the trajFrag list ;;
#|    -id :-VMD mol Id pointing to a particular trajInfo array element .
#|      -acceptable values :-"top" :-VMD's current top molecule is used ;;
#|      -default value :-"top" ;;
#|    -args (variable arguments) :
#|      -exclude, except :
#|        -list of fragIds to be excluded from the fragment specification .
#|        -acceptable values :
#|          -the same as l_fragId except for "all" and "exclude" keywords ;;
#|      -loSt, channelIdi, log :-output stream for log messages .
#|        -default value :-stdout ;;;;;
proc trajFragSpec {l_fragId {id "top"} args} {
# global variables
  global trajInfo
# namespace variables
  variable indTFL
# default values for variables and arguments
  set loSt stdout; set exclude {}
  if {$id == "top"} {set id [molinfo top]}
  if {$id == -1} {puts $loSt "trajFragSpec: No trajInfo loaded."; return ""}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch [string tolower $arg] {
          "exclude" - "except" - "excl" - "exfragid" {set exclude $val}
          "lost" -  "channelid" - "log" {set loSt $val}
          default {puts $loSt "trajFragSpec: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "trajFragSpec: Odd number of optional arguments! args: $args"
      return ""
      }
# checks for "exclude" keywords
  set pos -1
  if {[lsearch $l_fragId "exclude"] >= 0} {
    set pos [lsearch $l_fragId "exclude"]
  } elseif {[lsearch $l_fragId "except"] >= 0} {
    set pos [lsearch $l_fragId "except"]
  } elseif {[lsearch $l_fragId "excl"] >= 0} {
    set pos [lsearch $l_fragId "excl"]
  } elseif {[lsearch $l_fragId "exFragId"] >= 0} {
    set pos [lsearch $l_fragId "exFragId"]
    }
  if {$pos == 0} {
    set exclude [lrange $l_fragId 1 end]
    set l_fragId "all"
  } elseif {$pos > 0} {
    set exclude [lrange $l_fragId [expr {$pos + 1}] end]
    set l_fragId [lrange $l_fragId 0 [expr {$pos - 1}]]
    }
#  puts $loSt "\ntrajFragSpec: user fragId: $l_fragId"
#  if {$exclude != {}} {puts $loSt "trajFragSpec: user exclude list: $exclude"}
# build up fragNamesTmp, simNamesTmp, and namesTmp lists
  set fragNamesTmp {}
  set simNamesTmp {}
  set namesTmp {}
  for {set tfi 0} {$tfi < [llength $trajInfo($id,trajFrag)]} {incr tfi} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
    lappend fragNamesTmp [lindex $l_trajFrag $indTFL(fragName)]
    lappend simNamesTmp [lindex $l_trajFrag $indTFL(simName)]
    if {($trajInfo($id,name) == [lindex $l_trajFrag $indTFL(fragName)]) || \
        ($trajInfo($id,name) == [lindex $l_trajFrag $indTFL(simName)])} {
        lappend namesTmp [lindex $l_trajFrag $indTFL(fragName)]
      }
    }
# decode exclude list into fragNames
  set l_fragNameEx {}
  foreach excl $exclude {
    if {$excl == "name"} {
      set l_fragNameEx [concat $l_fragNameEx $namesTmp]
    } elseif {$excl == "first"} {
      lappend l_fragNameEx [lindex [lindex $trajInfo($id,trajFrag) 0] \
                                   $indTFL(fragName)]
    } elseif {($excl == "last") || ($excl == "end")} {
      lappend l_fragNameEx [lindex [lindex $trajInfo($id,trajFrag) end] \
                                   $indTFL(fragName)]
    } elseif {($excl == "to") || ($excl == "-")} {
      lappend l_fragNameEx "to"
    } elseif {[string is integer $excl]} {
      if {($excl >= 0) && ($excl < [llength $trajInfo($id,trajFrag)])} {
        lappend l_fragNameEx [lindex [lindex $trajInfo($id,trajFrag) $excl] \
                                     $indTFL(fragName)]}
    } elseif {[string first "end-" $excl] == 0} {
      lappend l_fragNameEx [lindex [lindex $trajInfo($id,trajFrag) $excl] \
                                   $indTFL(fragName)]
    } elseif {[lsearch $simNamesTmp $excl] >= 0} {
# no keyword used, then checks for a simName specified in excl
      for {set isn 0} {$isn < [llength $simNamesTmp]} {incr isn} {
        if {[lindex $simNamesTmp $isn] == $excl} {
          lappend l_fragNameEx [lindex $fragNamesTmp $isn]
          }
        }
    } else {
# otherwise a fragName is assumed
      lappend l_fragNameEx $excl
      }
    }
# check for ranges with the "to" keyword for excluded fragIds
  set pos [lsearch $l_fragNameEx "to"]
  while {$pos >= 0} {
# positions for the starting and final range elements in the full list of names
    set prevPos [lsearch $namesTmp [lindex $l_fragNameEx [expr {$pos - 1}]]]
    set postPos [lsearch $namesTmp [lindex $l_fragNameEx [expr {$pos + 1}]]]
    set l_fragNameEx [concat [lrange $l_fragNameEx 0 [expr {$pos - 1}]] \
                             [lrange $fragNamesTmp $prevPos $postPos] \
                             [lrange $l_fragNameEx [expr {$pos + 1}] end]]
    set pos [lsearch $l_fragNameEx "to"]
    }
# decode l_fragId list into fragNames
  set l_fragName {}
  foreach fragId $l_fragId {
    if {$fragId == "all"} {
      set l_fragName $fragNamesTmp
      break
    } elseif {$fragId == "name"} {
      set l_fragName [concat $l_fragName $namesTmp]
    } elseif {$fragId == "first"} {
      lappend l_fragName [lindex [lindex $trajInfo($id,trajFrag) 0] \
                                 $indTFL(fragName)]
    } elseif {($fragId == "last") || ($fragId == "end")} {
      lappend l_fragName [lindex [lindex $trajInfo($id,trajFrag) end] \
                                 $indTFL(fragName)]
    } elseif {($fragId == "to") || ($fragId == "-")} {
      lappend l_fragName "to"
    } elseif {[string is integer $fragId]} {
      if {($fragId >= 0) && ($fragId < [llength $trajInfo($id,trajFrag)])} {
        lappend l_fragName [lindex [lindex $trajInfo($id,trajFrag) $fragId] \
                                   $indTFL(fragName)]}
    } elseif {[string first "end-" $fragId] == 0} {
      lappend l_fragName [lindex [lindex $trajInfo($id,trajFrag) $fragId] \
                                 $indTFL(fragName)]
    } elseif {[lsearch $simNamesTmp $fragId] >= 0} {
# no keyword used, then checks for a simName specified in fragId
      for {set isn 0} {$isn < [llength $simNamesTmp]} {incr isn} {
        if {[lindex $simNamesTmp $isn] == $fragId} {
          lappend l_fragName [lindex $fragNamesTmp $isn]
          }
        }
    } else {
# otherwise a fragName is assumed
      lappend l_fragName $fragId
      }
    }
# check for ranges with the "to" keyword
  set pos [lsearch $l_fragName "to"]
  while {$pos >= 0} {
# positions for the starting and final range elements in the full list of names
    set prevPos [lsearch $fragNamesTmp [lindex $l_fragName [expr {$pos - 1}]]]
    set postPos [lsearch $fragNamesTmp [lindex $l_fragName [expr {$pos + 1}]]]
    set l_fragName [concat [lrange $l_fragName 0 [expr {$pos - 1}]] \
                           [lrange $fragNamesTmp $prevPos $postPos] \
                           [lrange $l_fragName [expr {$pos + 1}] end]]
    set pos [lsearch $l_fragName "to"]
    }
# build final list running over the trajFrag list
  set retList {}
  foreach fragName $fragNamesTmp {
    if {[lsearch $l_fragNameEx $fragName] < 0} {
      if {[lsearch $l_fragName $fragName] >= 0} {
        lappend retList $fragName
        }
      }
    }
  return $retList
  }   ;# trajFragSpec

}   ;# namespace eval userInfoLib

#|    - ;;
::userInfoLib::add_commands trajFragSpec
::userInfoLib::logMsg "added commands in trajFragSpec.tcl to userInfoLib namespace" 1



