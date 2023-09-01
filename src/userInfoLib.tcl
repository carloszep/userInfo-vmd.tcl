#|-userInfoLib.tcl :| 
#|  -tcl scripts library to manage VMD trajectories from MD simulations .
#|  -assists in incorporating user-specified information from simulations .
#|  -declares the namespace containing the varaibles and commands to be
#|   _ used by other analysis tools .
#|  -started from the userInfoLib_v.0.0.5 library .
#|  -dates :
#|    -created :-2023-05-03.Wed ;
#|    -modified :-2023-09-01.Fri ;;
#|  -authors and contributors :
#|    -Carlos Z. Gómez Castro ;
#|  -public software repositories :
#|    -https://github.com/carloszep/anMD-vmd.tcl ;
#|  -version information :
#|    -current version :-0.0.2 ;
#|    -changes in progress :
#|      -definition of the userInfoLib namespace .
#|      -added trajFragSpec.tcl to the userInfoLib namespace ;;
#|  -notes :
#|    -for the moment this library uses the global array selInfo ;
#|  -source :
#|    -logLib.tcl ;
source logLib.tcl

#|  -namespace eval userInfoLib :
namespace eval userInfoLib {

#|    -import :
#|      -::logLib::* ;
  namespace import ::logLib::*

#|    -export commands :
#|      -get_maxTrajSizeGB .
#|      -set_maxTrajSizeGB  .
#|      -loadedTrajSize ;
namespace export get_maxTrajSizeGB set_maxTrajSizeGB loadedTrajSize

#|    -variables :
#|      -maxTrajSizeGB .
#|      -currTrajSizeGB .
#|      -trajFragList .
#|      - ;
  variable maxTrajSizeGB 6.0
  variable currTrajSizeGB 0.0
  variable trajFragList [list fragName simName dcdFile timeStep dcdFreq \
               loadStep dcdSize iniTime finTime iniFrame finFrame frameTime]
  variable indTFL

#|    -commands :
#|      -proc init {} :
#|        -initialization of the userInfoLib namespace .
#|        -it also initializes the imported logLib namespace ;
  proc init {} {
    variable indTFL
    variable trajFragList
# initializes of the name and version of the logLib namespace :
    set_logName "userInfoLib"
    set_logVersion "0.0.2"
    set_logLevel 3
    set_logFileName ""
# configure string prefix for all log messages
    set_logPrefixStr "userInfo: "
# incorportates the userInfo list of commands to the logLib list
    add_commands [list get_maxTrajSizeGB set_maxTrajSizeGB \
                       loadedTrajSize get_trajFragList]
    set i 0
    foreach elem $trajFragList {
      set indTFL($elem) $i
      incr i
      }
    logMsg "initialized indices array indTFL from trajFragList..." 3
    logMsg "trajFragList: [get_trajFragList]" 2
    }
#|      -proc get_maxTrajSizeGB {} :
#|        - ;
  proc get_maxTrajSizeGB {} {
    variable maxTrajSizeGB
    return $maxTrajSizeGB
    }

#|      -proc set_maxTrajSizeGB {trajSize} :
#|        - ;
  proc set_maxTrajSizeGB {trajSize} {
    variable maxTrajSizeGB
    set maxTrajSizeGB $trajSize
    }

#|      -proc loadedTrajSize {} :
#|        - ;
  proc loadedTrajSize {} {
    variable currTrajSizeGB
    return $currTrajSizeGB
    }

#|      -proc get_trajFragList {} :
#|        - ;
  proc get_trajFragList {} {
    variable trajFragList
    return $trajFragList 
    }

#|      - ;
#|    - ;
  }   ;# namespace eval userInfoLib

#|    -reports to log the initial configuration .
::userInfoLib::logMsg "initialized [::userInfoLib::get_logName_version]" 1
::userInfoLib::logMsg "Manage VMD trajectories from MD simulations." 1
::userInfoLib::logMsg "log path: [::userInfoLib::get_logPath]" 2
::userInfoLib::logMsg "log output to: [::userInfoLib::get_logOutputStream]" 2
::userInfoLib::logMsg "output level: [::userInfoLib::get_logLevel]" 1
::userInfoLib::logMsg "output file for log: [::userInfoLib::get_logFileName]" 2
::userInfoLib::logMsg "print to screen: [::userInfoLib::get_logScreen]" 2
::userInfoLib::logMsg "list of commands: [::userInfoLib::list_commands]" 2
#|    -flush the output buffer .
::userInfoLib::logFlush
#|    -run the initialization proc .
::userInfoLib::init
#|    -external procedures to be added to userInfoLib :
#|      -trajFragSpec.tcl ;
source trajFragSpec.tcl

#|  - ;| eof



