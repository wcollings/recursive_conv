# SARA_v7
# How this is going to work is, we'll have a "step" function
# that computes the next value, but doesn't save it anywhere.
# instead, it just returns the result and moves on. The solver
# is going to try a bunch of different values to see which one
# converges mathematically the best, or is the most stable, etc.
# So our framework needs to not save those intermediary results
# instead, there's a function at the bottom called "accept"
# that actually saves that result.

package require math::complexnumbers
namespace import ::math::complexnumbers::*

proc SARA:Init {instance impedanceFlag k0 k0j k1 kj1 k2 kj2 k3 kj3 k4 kj4 sig1 sigj1 sig2 sigj2 sig3 sigj3 sig4 sigj4 {v_IC 0}} {
    global SARA
    global solver
	global outFile
    set sig [list $sig1 $sig2 $sig3 $sig4]
    set sigj [list $sigj1 $sigj2 $sigj3 $sigj4]
    set i 1
    if {$k1 == 0 && $kj1 == 0} {
		puts "k1 is 0, not allowed"
	}
    set order 3
    while {$i < 4} {
		if {[lindex $sig $i] == 0 && [lindex $sigj $i] == 0} {
			set order [- $i 1]
			break
		}
		incr i
	}
    puts "Found order: $order"
	if {impedanceFlag == 0}
	{
		puts "Using inductance calculations..."
		set SARA($instance,"impedanceFlag") 0
	}
	else {
		puts "Using impedance calculations..."
		set SARA($instance,"impedanceFlag") 1
	}
	set SARA($instance,"order") $order
	set SARA($instance,"k0") [complex $k0 $k0j]
	set SARA($instance,"k1") [complex $k1 $kj1]
	set SARA($instance,"k2") [complex $k2 $kj2]
	set SARA($instance,"k3") [complex $k3 $kj3]
	set SARA($instance,"k4") [complex $k4 $kj4]
	set SARA($instance,"sig1") [complex $sig1 $sigj1]
	set SARA($instance,"sig2") [complex $sig2 $sigj2]
	set SARA($instance,"sig3") [complex $sig3 $sigj3]
	set SARA($instance,"sig4") [complex $sig4 $sigj4]

    set solver($instance,"num_iter") 0
	set solver($instance,"prev_time") [complex 0 0]
	set solver($instance,"tt") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Integ") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Prev") [complex 0 0]
	set solver($instance,"xx_0") [complex 0 0]
	set solver($instance,"yy") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"yy_Final") [complex 0 0]
	
	if {$v_IC != 0} {
		puts "Found voltage I.C.: $v_IC"
		set temp [complex $v_IC 0]
		set solver($instance,"xx_0") $temp
		for {set i 0} {$i < 4} {incr i} {
			set [lindex solver($instance,"xx") $i] $temp
		}
		set temp [complex 1 0]
		for {set i 0} {$i < 4} {incr i} {
			set [lindex solver($instance,"tt") 0] $temp
		}
		
		SARA:Accept $instance -2 $v_IC 
		puts "Accept 1"
		SARA:Accept $instance -1 $v_IC 
		puts "Accept 2"
		SARA:Accept $instance 0 $v_IC 
		puts "Accept 3"
	}
	puts "done with Init"
}

proc SARA:Phi {sigma_i delta_n} {
	return [exp [* $sigma_i $delta_n]]
	# expr {exp($i*$n)}
}

proc SARA:zeta {sigma_i delta_n} {
	return [* [* $sigma_i $delta_n] [complex -1 0]]
	# expr {-$i*$n}
}

proc SARA:q0 {sigma_i delta_n n} {
	if {[real $delta_n] == 0 } {
		return $delta_n
	}
	return [* [/ $delta_n [SARA:zeta $sigma_i $delta_n]] [- [complex 1 0] [SARA:Phi $sigma_i $delta_n]]]
}

proc SARA:q1 {sigma_i delta_n n} {
	if {[real $delta_n] == 0} {
		return $delta_n;
	}
	set phi [SARA:Phi $sigma_i $delta_n]
	set zi [SARA:zeta $sigma_i $delta_n]
	set c1 [complex 1 0]
	set c1n [complex -1 0] 
	set c2 [complex 2 0]
	set c3 [complex 3 0]
	switch $n {
	0 	{
			#return [/ [* $delta_n [* $c2 [- $zi $c1]]] [* $zi $zi]]
			return [/ [* $delta_n [+ $c1n [+ $zi $phi]]] [* $zi $zi]]
			#return [/ [* [complex 1 0] $delta_n] $c2]
			#return [complex 0 0]
		}
	1	{
			#return $delta_n
			return [/ [* $delta_n [- $c1 [* [+ $c1 $zi] $phi]]] [* $zi $zi]]
			#return [* [/ [* [complex 1 0] $delta_n] $c2] [- $c1 $zi]]
			#return [complex 0 0]
		}
	}
}

proc SARA:q2 {sigma_i delta_n n} {
	if {$delta_n==[complex 0 0]} {
		return $delta_n
	}
	set phi [expr {SARA:Phi $sigma_i $delta_n}]
	set zi [expr {SARA:zeta $sigma_i $delta_n}]
	switch $n {
		0 	{
				return [expr {$delta_n/[2 * pow $zi 3]*[2-[3*$zi]+[2*pow $zi 2] - [2-$zi]*$phi]}]
			}
		1	{
				return [expr {$delta_n/[pow $zi 3]*[-2*[1-$zi]+[2-pow $zi 2]*$phi]}]
			}
		2	{
				return [expr {$delta_n/[2*pow $zi 3]*[2-$zi-[2+$zi]*$phi]}]
			}
	}
}

proc SARA:Step {instance time input} {  
	#puts "step"
    global SARA
    global solver
	set impedanceFlag $SARA($instance,"impedanceFlag")
    set time [complex $time 0]
    set xx_Curr [complex $input 0]
    set delta_n [- $time $solver($instance,"prev_time")];      # delta_n = current_time - solver(last_recorded_time)
    set xx_Integ_Section [/ [* $delta_n [+ $xx_Curr $solver($instance,"xx_Prev")]] [complex 2 0]];    # integral from n-1 to n = delta_n * (x[n] + x[n-1]) / 2
    set xx_Integ_Full [+ $xx_Integ_Section [lindex $solver($instance,"xx_Integ") 0]];  	             # integral from 0 to n = xx_Integ (0 to n-1) + xx_Integ_Section (n-1 to n)
	set order $SARA($instance,"order")
	if {impedanceFlag == 0} {
    	set final [* $xx_Integ_Full $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)
    	for {set i 1} {$i <= [+ $order 1]} {incr i} {
			set sig_i $SARA($instance,"sig$i")
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			set temp [+ $temp [* $k_i [* $q $xx_Integ_Full]]];        # calculate using most recent x integral, which is not stored
			for {set j 1} {$j <= $order} {incr j} {
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") [- $j 1]]]]];     # use $j-1 since $j starts from 1
        }
        set final [+ $final $temp]
    }
	else {
    	set final [* $xx_Curr $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)
    	for {set i 1} {$i <= [+ $order 1]} {incr i} {
			set sig_i $SARA($instance,"sig$i")
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			set temp [+ $temp [* $k_i [* $q $xx_Integ_Full]]];        # calculate using most recent x integral, which is not stored
			for {set j 1} {$j <= $order} {incr j} {
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") [- $j 1]]]]];     # use $j-1 since $j starts from 1
        }
        set final [+ $final $temp]
	}
	}
	set order $SARA($instance,"order")
    set final [* $xx_Integ_Full $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)
    for {set i 1} {$i <= [+ $order 1]} {incr i} {
        set sig_i $SARA($instance,"sig$i")
        set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
        set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
        set q_function "SARA:q$order"
        set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
        set temp [+ $temp [* $k_i [* $q $xx_Integ_Full]]];        # calculate using most recent x integral, which is not stored
        for {set j 1} {$j <= $order} {incr j} {
            set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
            set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") [- $j 1]]]]];     # use $j-1 since $j starts from 1
        }
        set final [+ $final $temp]
    }
    return [real $final]
}

proc SARA:Accept {instance time input} {
	#puts "accept"
    global SARA
    global solver
	global outFile
    set time [complex $time 0]
    set xx_Curr [complex $input 0]
    
	incr solver($instance,"num_iter")
    set delta_n [- $time $solver($instance,"prev_time")];     # delta_n = current_time - solver(last_recorded_time)
    set xx_Integ_Section [/ [* $delta_n [+ $xx_Curr $solver($instance,"xx_Prev")]] [complex 2 0]];    # integral from n-1 to n = delta_n * (x[n] + x[n-1]) / 2
    set xx_Integ_Full [+ $xx_Integ_Section [lindex $solver($instance,"xx_Integ") 0]];               # integral from 0 to n = xx_Integ (0 to n-1) + xx_Integ_Section (n-1 to n)
    set order $SARA($instance,"order")
    set tempList [lrange [linsert $solver($instance,"xx_Integ") 0 $xx_Integ_Full] 0 [+ $order 1]];   # add most recent integral to front of list and shift all elements back
	set solver($instance,"xx_Integ") $tempList
	set tempList [lrange [linsert $solver($instance,"tt") 0 $delta_n] 0 [+ $order 1]];   # add most recent time to front of list and shift all elements back
	set solver($instance,"tt") $tempList
    set solver($instance,"prev_time") $time
    
    set final [* $xx_Integ_Full $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)
	# puts $outFile "[real $time],[real $xx_Integ_Full]"
    for {set i 1} {$i <= [+ $order 1]} {incr i} {
        set sig_i $SARA($instance,"sig$i")
        set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
        set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
        set q_function "SARA:q$order"
        for {set j 0} {$j <= $order} {incr j} {
            set q [$q_function $sig_i $delta_n $j]
			#puts $q
            set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") $j]]]];     # use $j since most recent value is already stored
        }
        set solver($instance,"yy") [lreplace $solver($instance,"yy") [- $i 1] [- $i 1] $temp];              # lreplace using $i-1 because $i starts from 1, not 0
        set final [+ $final $temp]
    }
    set solver($instance,"xx_Prev") $xx_Curr;       # Update xx_Prev for next iteration
    set solver($instance,"yy_Final") $final;        # Store final calculated value
}