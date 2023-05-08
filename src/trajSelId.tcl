#|-namespace eval trajSelId :
#|  -contains variables and procedures used to iterate along trajectory frames
#|   _ specified by a selId as an selInfo array key .
#|  -usage :
#|    -while {[trajSelId::iterate] <selId> <varArgs>]} { ... } :
#|      -the while loop will iterate the number once for each frame in the
#|       _ trajectory according to the selInfo specifications .
#|      -to access the current frame number or the equivalent simulation time
#|       _ value inside the while loop use commands [trajSelId::getFrame] and
#|       _ [trajSelId::getDatX] ;;
#|  -list of namespace variables :
#|    -selId :-selection Id array key to access the selInfo array .
#|      -identifies the user-defined selection (of frames) to go over ;
#|    -id :-VMD mol Id specified in the selId ;
#|    -loSt :-output stream for log messages .
#|      -default value :-stdout ;;
#|    -errorFlag :-if "1" no iteration will be done and iterate returns "0" ;
#|    -finishedFlag :-if "1" no iteration will be done and iterate returns "0" ;
#|    -loadTraj :-if "1" each trajFrag have to be loaded ;
#|    -frameCount :-frame iteration counter ;
#|    -xDatComm :-command used to calculate the X data value ;
#|    -fragIds :-list of fragIds considered ;
#|    -currFragInd :-index of the current fragment ;
#|    -currFrag :-fragName of the current fragment ;
#|    -currFrame :-current frame along the trajectory .
#|      -this will be consecutive within a traj fragment, but will change
#|       _ every time the current fragment changes ;
#|    -finFrame :-points to the last frame of the current trajFrag or to the
#|     _ last frame indicated in the selId ;
#|    -usrFinFrame :-points to the final frame corresponding to the "last"
#|     _ selInfo array key specified in the selId .
#|      -this has a value different from -1 only when the current fragment
#|       _ includes the indicated frame ;
#|    -iniTime :-simulation time corresponding to the first frame to consider ;
#|    -finTime :-simulation time corresponding to the last frame to consider ;
#|    -loadStep :-loadStep as indicated in the selId ;
#|    -keepFragIds :-list of traj fragments to be kept (avoid deleting them) ;
#|    -upd :-update otion in the selInfo array ;
#|    -step :-step option in the selInfo array .
#|      -affect the frequency and number of frames gone over ;
#|    -seqTime :-indicates that the traj is sequential in time .
#|      -if "0" (as in an equilibration) :
#|        -frame number is reported instead of time as X data .
#|        -No user-specified first and last times are considered ;
#|      -note that this value is only set as default in the namespace variable
#|       _ declaration, but is not default-reset in the init proc :
#|        -it is reset to 1 at the end of the iteration cycle ;;;
#|  -list of namespace procedures :
#|    -init :-initializes namespace varibles ;
#|    -changeFrag :-updates frame pointers according to a new traj fragment ;
#|    -iterate :-perform a iteration pointing to the next traj frame ;
#|    -setSeqTime :-sets sequential-time frags evaluation to use the sim time
#|     _ (if possible) or the frame number in the getDatX command ;
#|    -setKeepFragIds :-sets user-specified keepFragIds ;
#|    -getFrame :-returns the current frame ;
#|    -getDatX :-returns the value of the simTime or simFrame corresponding
#|     _ to the current frame ;
#|    -getErrorFlag :-returns the value of the "errorFlag" variable ;;;
namespace eval trajSelId {
  global trajInfo selInfo
  namespace export getFrame getDatX iterate getErrorFlag
  variable selId ""
  variable id -1
  variable loSt stdout
  variable errorFlag 0
  variable finishedFlag 0
  variable loadTraj 0
  variable frameCount 0
  variable xDatComm {simTime $currFrame $id}
  variable fragIds {}
  variable currFragInd 0
  variable currFrag ""
  variable currFrame -1
  variable finFrame -1
  variable usrFinFrame -1
  variable iniTime ""
  variable finTime ""
  variable loadStep "trajFrag"
  variable keepFragIds {}
  variable upd 0
  variable step 1
  variable seqTime 1

#|-proc trajSelId::init {elIdUsr {keepFragIdsUsr {}} {exFragIds {}}
#|                                                 _ {loStUsr stdout}} :
#|  -initializes most of the variables according to the user-specified selId .
#|  -notes :-the 'lbl' and 'desc' varaibles from getSelInfo is not used ;;
  proc init {selIdUsr {keepFragIdsUsr {}} {exFragIds {}} {loStUsr stdout}} {
    variable selId $selIdUsr
    variable id -1
    variable loSt $loStUsr
    variable errorFlag 0
    variable finishedFlag 0
    variable loadTraj 0
    variable frameCount 0
    variable xDatComm {simTime $currFrame $id}
    variable fragIds {}
    variable currFragInd 0
    variable currFrag ""
    variable currFrame 0
    variable finFrame -1
    variable usrFinFrame -1
    variable iniTime ""
    variable finTime ""
    variable loadStep "trajFrag"
    variable keepFragIds {}
    variable upd 0
    variable step 1
    variable seqTime
    if {$selId == ""} {return}
    getSelInfo $selId selTxt id frag frm iniTime finTime step upd loadStep \
                      lbl desc
    set keepFragIds [trajFragSpec $keepFragIdsUsr $id loSt $loSt]
    if {$id == "top"} {
      puts $loSt "Error: no loaded trajInfo"
      set errorFlag 1
      return
      }
    if {[getTrajProp "loaded" $id $frag exclude $exFragIds loSt $loSt]} {
      set loadTraj 0} else {set loadTraj 1}
    if {$seqTime} {
      set seqTime [getTrajProp "seqTime" $id $frag excl $exFragIds loSt $loSt]
      }
    if {!$seqTime} {
      set xDatComm {set frameCount}
      puts $loSt "Non-sequential time: using frame number instead."
      }
# setting up current fragment, current frame, and final frame
    set fragIds [trajFragSpec $frag $id exclude $exFragIds loSt $loSt]
    changeFrag $currFragInd
    }   ;# trajSelId::init

#|-proc changeFrag {fragInd} :
#|  -updates the pointers to initial and final frames and to fragments .
#|  -loads the updated traj fragment if necessary ;
  proc changeFrag {fragInd} {
    variable loadTraj
    variable id
    variable fragIds
    variable loadStep
    variable loSt
    variable currFragInd $fragInd
    variable currFrag
    variable iniTime
    variable finTime
    variable currFrame
    variable finFrame
    variable usrFinFrame
    variable seqTime
# seek the fragment containing the user-specified iniTime
    if {$seqTime} {
      for {set fragi $currFragInd} {$fragi < [llength $fragIds]} {incr fragi} {
        if {($iniTime>=[getTrajProp "iniTime" $id [lindex $fragIds $fragi]])&&\
            ($iniTime<=[getTrajProp "finTime" $id [lindex $fragIds $fragi]])} {
          set currFragInd $fragi; break}}}
    set currFrag [lindex $fragIds $currFragInd]
    if {$loadTraj} {trajLoad $id $currFrag loadStep $loadStep loSt $loSt}
    set currFrame [getTrajProp "iniFrame" $id $currFrag loSt $loSt]
    set finFrame [getTrajProp "finFrame" $id $currFrag loSt $loSt]
# correct the initial and final frames according to the user-specified selId
    if {$seqTime} {
      set usrIniFrame [simFrame $iniTime $id $currFrag]
      set usrFinFrame [simFrame $finTime $id $currFrag]
      if {$usrIniFrame != -1} {
        if {$currFrame < $usrIniFrame} {set currFrame $usrIniFrame}
        }
      if {$usrFinFrame != -1} {
        if {$finFrame > $usrFinFrame} {set finFrame $usrFinFrame}
        }
      }
    }   ;# changeFrag

#|-proc trajSel::iterate {selIdUsr args} :
#|  -increases the currFrame variable except if it is the first frame .
#|  -arguments :
#|    -selIdUsr :-key name to access the selInfo array .
#|      -array keys considered :
#|        -molId, fragId, first, last, step, and loadStep ;;
#|    -args (variable arguments) :
#|      -"keepFragId", "keep" :
#|        -fragId specification to avoid deleting specific traj fragments .
#|        -default value :-{} ;;
#|      -"exFragId", "exclude", "except", "excl" :
#|        -"exclude" fragId specification for the trajFragSpec proc .
#|        -default value :-{} ;;
#|      -"loSt", "channelId", "log" :-output stream for log messages .
#|        -default value :-stdout ;;;;
#|  -note :-not documented in detail yet ;;
  proc iterate {selIdUsr args} {
    variable selId
    variable id
    variable loSt
    variable step
    variable errorFlag
    variable finishedFlag
    variable loadTraj
    variable frameCount
    variable currFragInd
    variable currFrag
    variable currFrame
    variable finFrame
    variable fragIds
    variable usrFinFrame
    variable keepFragIds
# default values for arguments
  set loStUsr stdout; set keepFragIdUsr {}; set exFragId {}
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    if {[llength $args] > 0} {
      foreach {arg val} $args {
        switch $arg {
          "keepFragId" - "keep" {set keepFragIdUsr $val}
          "exFragId" - "exclude" - "except" - "excl" {set exFragId $val}
          "loSt" -  "channelId" - "log" {set loStUsr $val}
          default {puts $loSt "trajSelId::iterate: argument unkown: $arg"}
          }
        }
      }
    } else {   ;# odd number of arguments
      puts $loSt "trajSelId::iterate: Odd number of arguments! args: $args"
      return ""
      }
    if {$errorFlag || $finishedFlag} {puts $loSt "Error Flag"; return 0}
    if {$selId != $selIdUsr}  {
      init $selIdUsr $keepFragIdUsr $exFragId $loStUsr
      return 1
      }
    incr currFrame $step
    incr frameCount
    if {$currFrame <= $finFrame} {   ;# running over the current frag
      return 1
    } else {   ;# end of the current frag exceeded
# delete the current (just completed) fragment if it was not previouly loaded
#  and it was not included in the keepFragIds list
      if {($loadTraj) && \
        ([lsearch $keepFragIds $currFrag] < 0)} {
          trajDelete $id $currFrag loSt $loSt
          }
      if {$usrFinFrame == -1} {
# the frag with the user finTime was not reached yet
        if {$currFragInd == [expr {[llength $fragIds] - 1}]} {
# no more fragments to load
          set seqTime 1
          set selId ""
          if {$loadTraj} {trajDelete $id $keepFragIds loSt $loSt}
          return 0
        } else {
# load the next fragment
          changeFrag [expr {$currFragInd + 1}]
          return 1
          }
      } else {
# the current frag contains the user-specified finTime
        set seqTime 1
        set selId ""
        if {$loadTraj} {trajDelete $id $keepFragIds loSt $loSt}
        return 0
        }
      }
    }   ;# trajSelId::iterate

#|-proc setSeqTime {val} :
#|  -sets seqTime usage for getDatX to report frame numbers instead of
#|   _ simulation time if val is "0" ;
  proc setSeqTime {val} {
    variable seqTime $val
    }

#|-proc trajSelId::getFrame {} :
#|  -returns the current frame ;
  proc getFrame {} {
    variable currFrame
    return $currFrame
    }

#|-proc trajSelId::getDatX {} :
#|  -returns the current simTime or simFrame in the trajectory ;
  proc getDatX {} {
    variable xDatComm
    variable id
    variable currFrame
    variable frameCount
    return [eval $xDatComm]
    }

#|-proc trajSelId::setKeepFragIds {l_frag} :
#|  -sets the keepFragIds variable .
#|  -used when a particular fragment or frame must be loaded during the whole
#|   _ iteration process (i.e. a reference frame for analysis
#|   _ of a long trajectory comprised of several fragments) ;
  proc setKeepFragIds {l_frag} {
    variable keepFragIds $l_frag
#    puts $loSt "traj fragments to be kept: $keepFragIds"
    }

#|-proc trajSelId::getErrorFlag {} :
#|  -returns the value of the "errorFlag" variable ;
  proc getErrorFlag {} {
    variable errorFlag
    return $errorFlag
    }   ;# trajSelId::getErrorFlag

  }   ;# namespace trajSelId



proc testNS {selId} {
  namespace import trajSelId::getFrame trajSelId::getDatX
  puts "testing ::trajSelId:: selId: $selId"
  puts "current namespace: [namespace current]"
  puts "  going over all frames..."
  while {[trajSelId::iterate $selId keep iniProdSim]} {
    puts "frame: [getFrame]   xDat: [getDatX]"
    }
  puts "Error flag: [trajSelId::getErrorFlag]"
  }

