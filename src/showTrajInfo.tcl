#|-proc showTrajInfo {id {l_fragId "all"} args} :
#|  -reports the information of a trajectory stored in trajInfo array .
#|  -a list of fragName or simName names can be specified to delimit info .
#|  -an optional return value can be choose .
#|  -arguments :
#|    -id :-vmd mol ID associated to the trajectory .
#|      -may be specified as "top" to use the top molecule ;
#|    -l_fragId :
#|      -user trajectory fragment(s) specification .
#|      -acceptable values :-see trajFragSpec procedure ;
#|      -default value :-"all" ;;
#|    -args (variable arguments) :
#|      -loSt, channelId, log :-output stream for log messages .
#|        -default value :-stdout ;;
#|      -"exlude", "except", "excl", "exFragId" :
#|        -fragId exclude list for the trajFragSpec proc .
#|        -default value :-{} ;;;
#|  -notes :
#|    -it could report also the selIds associated to this molId ;;
proc showTrajInfo {id {l_fragId "all"} args} {
  global ind trajInfo trajFragList currTrajSize
  global maxTrajSizeGB
# default values for arguments
  set loSt stdout; set exclude {}; set procName [lindex [info level 0] 0]
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "loSt" -  "channelId" - "log" {set loSt $val}
          "exclude" - "except" {set exclude $val}
          default {puts $loSt "$procName: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "$procName: Odd number of variable arguments! args: $args"
      return ""
      }
  if {$id == "top"} {set id [molinfo top]}
  if {$id == -1} {puts $loSt "procName: No trajInfo loaded."; return}
  set l_fragId [trajFragSpec $l_fragId $id exclude $exclude loSt $loSt]
# print information
  puts $loSt "\n+++ Global simulation parameters: +++"
  puts $loSt "  VMD Id: $id"
  puts $loSt "  sim name: $trajInfo($id,name)"
  puts $loSt "  description: $trajInfo($id,desc)" 
  puts $loSt "  psf file: $trajInfo($id,psf)"
  puts $loSt "  total traj size (GB): $trajInfo($id,dcdSize)"
  puts $loSt "  global load step: $trajInfo($id,loadStep)"
  puts $loSt "  frames currently loaded: [molinfo $id get numframes]"
  puts $loSt "  memory used: $currTrajSize of $maxTrajSizeGB GB"
  set nFrag [llength $trajInfo($id,trajFrag)]
  puts $loSt "  total number of trajectory fragments: $nFrag"
  puts -nonewline $loSt "\nfrag:"
  foreach tfli $trajFragList {puts -nonewline $loSt "  $tfli"}
  puts $loSt ""
  set fragIdIniFrame ""
  set fragIdDcdSize 0.0
  set loaded 0
  for {set tfi 0} {$tfi < $nFrag} {incr tfi} {
    set l_trajFrag [lindex $trajInfo($id,trajFrag) $tfi]
# extract and assign values for new variables from the trajFrag list
    foreach elem $trajFragList {set $elem [lindex $l_trajFrag $ind($elem)]}
    foreach fragId $l_fragId {
      if {($fragId == $fragName)||($fragId == $simName)||($fragId == "all")} {
        puts -nonewline $loSt "$tfi"
        foreach fragElem $l_trajFrag {puts -nonewline $loSt "  $fragElem"}
        puts $loSt ""
        if {$iniFrame == "unk"} {break}
        if {$fragIdIniFrame == ""} {
# this is only set the first time a requested fragment was found to be loaded
          set fragIdIniFrame $iniFrame
          set fragIdIniTime $iniTime
          set loaded 1
          }
        set fragIdFinFrame $finFrame
        set fragIdFinTime [simTime $finFrame $id]
# calcualtes the amount of memory released
        set timeStep [lindex $l_trajFrag $ind(timeStep)]
        if {$timeStep == 0} {
           set loadStepUsed 1
        } else {
          set loadStepUsed \
            [expr {round($frameTime*1000000.0/$dcdFreq/$timeStep)}]}
        set fragIdDcdSize \
              [expr {$fragIdDcdSize + $dcdSize/$loadStepUsed}]
        break
        }
      }
    }
# report information only for the specified trajFrags if loaded
  if {$loaded} {
    puts $loSt "\nSimulation parameters for loaded fragment(s)"
    puts $loSt "  iniFrame: $fragIdIniFrame"
    puts $loSt "  finFrame: $fragIdFinFrame"
    puts $loSt "  iniTime (ns): $fragIdIniTime"
    puts $loSt "  finTime (ns): $fragIdFinTime"
    puts $loSt "  accumulated DCD size (GB): $fragIdDcdSize"
    }
  }   ;# showTrajInfo

