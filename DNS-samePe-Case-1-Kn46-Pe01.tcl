#--------------------------------------------------------------------------#
# CASE-1: dp=5 nm, T=600K
# gamma=0.103
#--------------------------------------------------------------------------#
# PART-1---> DEFINE THE MODELS AND CONSTANTS OF EQUATIONS
#--------------------------------------------------------------------------#
set timediff [clock microseconds]

# Environment variables that are used:
# --NICHT-- TCLNUMPARTS : number of particles per core (def: 1000)
# TCLOUTDIR : output directory (prefix)
# TCLSKIN : value for ESPResSo variable skin (def: 0.4)
# TCLMAXNUMCELLS: value for ESPResSo variable man_num_cells (def: 1000)
# TCLSTEPS: value for n_steps
# TCLNCYCLE: value for n_cycle


#-------------------BOX-----------------------------------#
# Box length should be same to keep l_vortex=100 sigma

set box_length 800
puts "box length: $box_length"
setmd box_l $box_length $box_length $box_length
# Particle volume densitiy should be same for all simulations.
set density 0.00625
set num_nodes [setmd n_nodes]
set particle 3200000

#-------------------OUTPUT FILE-----------------------------------#
set outdir DNS-samePe_c-1_Kn46-Pe01
puts "output directory: $outdir"


#-------------------CHECK POINTS-----------------------------------#
set checkpoint "checkpoint-"
set checkpoint_each_n_cycle 1

#-------------------PERFORMANCE---------------------------#
if { [info exists ::env(TCLMAXNUMCELLS) ] } {
    setmd max_num_cells $::env(TCLMAXNUMCELLS)
} else {
    setmd max_num_cells [expr $particle/$num_nodes]
}
puts -nonewline "max_num_cells: "
puts [setmd max_num_cells]

#-------------------VVERLET INTEGRATOR--------------------#
if { [info exists ::env(TCLNCYCLE) ] } {
    set n_cycle $::env(TCLNCYCLE)
} else {
    set n_cycle 100
}
if { [info exists ::env(TCLNSTEPS) ] } {
    set n_steps $::env(TCLNSTEPS)
} else {
    set n_steps 10
}
puts "n_cycle: ${n_cycle}, n_steps: ${n_steps}"
setmd time_step 0.001
puts "time_step [setmd time_step]"
# to adjust parameter skin for performance reasons
if { [info exists ::env(TCLSKIN) ] } {
    setmd skin $::env(TCLSKIN)
} else {
    setmd skin 0.8
}
puts -nonewline "skin: "
puts [setmd skin]

#---------------BOUNDARY_CONDITIONS-----------------------#
setmd periodic 1 1 1
#---------------LANGEVIN/LB THERMOSTAT--------------------#
set temp 1
set gamma 0.1027
thermostat off

#----------LENNARD JONES SHORT RANGED POTENTIAL-----------#
set eps 1.0
set sigma 1.0
set rcut 2.5
inter 0 0 lennard-jones $eps $sigma $rcut auto

#---------------------BONDED POTENTIALS-------------------#

inter 0 harmonic 800 1.0
# set of angle potentials
set angleres 181
for {set a 0} {$a < $angleres} {incr a} {
  inter [expr 1 + $a] angle_harmonic 800 [expr $a * [PI] / ($angleres-1)]
}

#------OUTPUT TRAJECTORY VOR VIEWING IN VMD---------------#
#set fl [open "|gzip -c - > ${outdir}case1.vtf.gz" a]

#---------------------OUTPIT-------------------#
set iter 0
if {[info exists checkpoint]} {
    set files [glob -nocomplain "$outdir$checkpoint\[0-9\]*.dat.gz"]
    foreach f $files {
	set replacement [subst {$outdir$checkpoint "" .dat.gz ""}]
	set number [string map $replacement $f]
	if {$iter<$number} {set iter $number}
    }
    if {[llength $files]>0} {
	puts "restarting with checkpoint $outdir$checkpoint$iter.dat.gz"
    }
}

#--------------------------------------------------------------------------#
# PART-2---> CEATE RANDOMLY PLACED PARTICLES
#--------------------------------------------------------------------------#
set timediff [clock microseconds]

if {(! [info exists checkpoint]) || (! [file exists "$outdir$checkpoint$iter.dat.gz"])} {
    #set particle [expr int(($box_length**3)*$density)]
    set total_time_random [time { for {set i 0} {$i<$particle} {incr i} {
	set posx [expr $box_length*[t_random]]
	set posy [expr $box_length*[t_random]]
	set posz [expr $box_length*[t_random]]
	part $i pos $posx $posy $posz type 0
    } } ]
}    

# write structure information
if {[info exists fl]} {
    if {(! [info exists checkpoint]) || (! [file exists "$outdir$checkpoint$iter.dat.gz"])} {
	writevsf $fl
    }
    # fix set of particles to output
    for {set i 0} {$i <$particle} {incr i} {
	lappend pidlist $i
	vtfpid $i
    }
}

if {(! [info exists checkpoint]) || (! [file exists "$outdir$checkpoint$iter.dat.gz"])} {
    #-------------------WARMUP--------------------------------#
    
    set total_time_warmup [time { 
	#set min 0
	set cap 20
	set warmupiter 0
	#while { $min < 1.0 && $warmupiter <= 10 } { }
	while { $warmupiter <= 10 } {
           
	    # limit forces to allow smooth adaption
	    inter forcecap $cap
	    integrate 40
	    #set min [analyze mindist 0 0]
	    incr cap 20
	    incr warmupiter 1
	    #puts "Warmup, mindist $min, cap $cap"
	    puts "Warmup, cap $cap"
	}
	# turn off the ljforcecap
	inter forcecap 0
	#puts "minimal particle distance is [analyze mindist 0 0]" # check Axel
    } ]

    set checkpoint_just_read 0
} else {
    # read checkpoint
    puts "reading checkpoint $outdir$checkpoint$iter.dat.gz"
    set f [open "|gzip -cd $outdir$checkpoint$iter.dat.gz" "r"]
    while {[blockfile $f read auto] != "eof"} {}
    catch {close $f}

    # update particles with binary stored positions and velocities
    #puts "reading checkpoint $outdir${checkpoint}binary-$iter.dat.gz"
    #set f [open "|gzip -cd $outdir${checkpoint}binary-$iter.dat.gz" "r"]
    #readmd $f
    #catch {close $f}

    # avoid overwriting initial checkpoint
    set checkpoint_just_read 1
}

proc write_checkpoint {iter} {
    global tcl_precision
    global outdir checkpoint
    set timediffchkpt [clock microseconds]
    set tcl_precision 10
    set f [open "|gzip -c - > $outdir$checkpoint$iter.dat.gz" "w"]
    blockfile $f write tclvariable iter
    blockfile $f write variable box_l
    blockfile $f write variable skin
    blockfile $f write variable max_num_cells
    blockfile $f write variable time_step
    blockfile $f write random
    blockfile $f write particles {id pos v}
    blockfile $f write bonds all
    close $f
    #set f [open "|gzip -c - > $outdir${checkpoint}binary-$iter.dat.gz" "w"]
    #writemd $f posx posy posz vx vy vz
    #close $f

    if {[info exists fl]} {
	set f [open "$outdir${checkpoint}$iter.vtf" "w"]
	writevsf $f
	writevcf $f folded
	close $f
    }
    puts "time writing checkpoint: [expr [clock microseconds] - $timediffchkpt]"
}

#--------------------------------------------------------------------------#
# PART-3---> MAIN SIMULATION
#--------------------------------------------------------------------------#
set timediff [clock microseconds]

thermostat langevin $temp $gamma

on_collision bind_three_particles 1 0 1 180

# looop
if {[catch {
    while { $iter<$n_cycle } {
	# writes current checkpoint
	if {[info exists checkpoint] && ($iter % $checkpoint_each_n_cycle == 0) && ($checkpoint_just_read == 0)} {
	
	    write_checkpoint $iter

	    puts "sort_particles $iter"
	    #sort_particles
	} else {
	    set checkpoint_just_read 0
	}

	puts [time {integrate $n_steps}]

	# write out coordinates for VMD
	if {[info exists fl]} {
	    writevcf $fl pids $pidlist
	}

	incr iter
    }
} err]} {
    puts "aborted with error $err"
    write_checkpoint "crash"
    exit -1
}

#--------------------------------------------------------------------------#
# PART-4---> POSTPROCESSING 
#--------------------------------------------------------------------------#
set timediff [clock microseconds]

# write final checkpoint:
if {[info exists checkpoint]} {
    write_checkpoint $iter
}
if {[info exists fl]} {
    close $fl
}

exit

