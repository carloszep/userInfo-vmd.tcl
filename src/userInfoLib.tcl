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

#|    -namespace variables :
#|      -maxTrajSizeGB :
#|        -maximum size of memory used by loaded dcd trajectores .
#|        -depends upon the SO and the physical RAM memory available .
#|        -it is recommended to not using more than ~75% of the total RAM :
#|          -for instance, with 4 GB of RAM, set maxTrajSizeGB to 3.0 ;;
  variable maxTrajSizeGB 6.0
#|      -currTrajSizeGB :
#|        -registers the amout of memory occupated by dcd trajectories
#|         _ loaded into VMD using the userInfoLib utilities .
#|        -before a traj fragment is loaded (i.e. by trajLoad), it is
#|         _verified that currTrajSizeGB does not exceed maxTrajSizeGB ;
  variable currTrajSizeGB 0.0
#|      -trajFragList :
#|        -list of trajectory-related property names used by userInfoLib .
#|        -default value :
#|          -{fragName simName dcdFile  timeStep dcdFreq   loadstep dcdSize
#|           _ iniTime finTime iniFrame finFrame frameTime} ;;
  variable trajFragList [list fragName  simName   dcdFile  timeStep  dcdFreq  \
                              loadStep  dcdSize   iniTime  finTime   iniFrame \
                              finFrame  frameTime ]
#|      -indTFL :
#|        -maps the index of the property names in the trajFragList .
#|        -it is initialized in the init command ;;
  variable indTFL

#|    -commands :
#|      -proc init {} :
#|        -initialization of the userInfoLib namespace .
#|        -it also initializes the imported logLib namespace ;
  proc init {} {
    variable indTFL
    variable trajFragList
# initializes of the name and version of the logLib namespace :
    set_logName "userInfo"
    set_logVersion "0.0.2"
    set_logLevel 3
    set_logFileName "stdout"
    logMsg "" 2
    set_logPrefixStr "[get_logName_version]: "
# printing introduction
    logMsg "userInfo library to manage Molecular Dynamics trajectories " 2
    logMsg "  and set up structural analysis in VMD." 2
    logMsg "Usage:" 2
    logMsg " 1. source the trajInfo script for a specific simulation." 2
    logMsg " 2. set the array selInfo for each atoms selections to be used." 2
    logMsg "" 2
# incorportates the userInfo list of commands to the logLib list
    logMsg "adding userInfo variables and commands to logLib:" 3
    logMsg "  commands: get_maxTrajSizeGB set_maxTrajSizeGB" 3
    logMsg "            loadedTrajSize get_trajFragList" 3
    logMsg "  variables: maxTrajSizeGB" 3
    add_commands [list get_maxTrajSizeGB set_maxTrajSizeGB \
                       loadedTrajSize get_trajFragList]
    add_variables [list maxTrajSizeGB]
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

#|    -run the initialization proc .
::userInfoLib::init
#|    -external procedures to be added to userInfoLib :
#|      -trajFragSpec.tcl ;
source trajFragSpec.tcl

#|  - ;| eof



