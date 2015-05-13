# TEST THE CIRCULAR FLOW FIELD

set box_length 800
setmd box_l $box_length $box_length $box_length
set particle1 100
set particle 1000000

set x 1
set y 1
set z 1

        set partid 0
for {set k 0} {$k<$particle1} {incr k} {
	for {set j 0} {$j<$particle1} {incr j} {
	    	for {set i 0} {$i<$particle1} {incr i} {
			part $partid pos [expr $i*$x] [expr $j*$y] [expr $k*$z] type 0
	        	set partid [expr $partid+1]
 		}
	}
}


setmd time_step 0.01
setmd skin 1
#------------------LANGEVIN THERMOSTAT--------------------#
set temp 1
set gamma 0.1027
thermostat langevin $temp $gamma
#---------------BOUNDARY_CONDITIONS-----------------------#
#...Periodic boundary cond. Later, this may be changed....#
setmd periodic 1 1 1
#---------------PARTICLE MOTION---------------------------#

integrate 0


exit

