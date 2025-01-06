# SARA_v8
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
interp alias {} ? {} set errorInfo

proc SARA:Init {instance impedanceFlag k0 k0j k1 kj1 k2 kj2 k3 kj3 k4 kj4 kn kjn sig1 sigj1 sig2 sigj2 sig3 sigj3 sig4 sigj4 {v_IC 0}} {
    global SARA
    global solver
	global outFile
	puts "[info nameofexecutable]"
	set outFile [open "SARA_v8.txt" a+]
	puts $outFile "_____________________new iteration - start of Init_____________________"
	puts $outFile "SARA:Init $instance $impedanceFlag $k0 $k0j $k1 $kj1 $k2 $kj2 $k3 $kj3 $k4 $kj4 $kn $kjn $sig1 $sigj1 $sig2 $sigj2 $sig3 $sigj3 $sig4 $sigj4 $v_IC"
	
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
	if {$impedanceFlag == 0} {
		puts "Using inductance calculations..."
		set SARA($instance,"impedanceFlag") 0
	} else {
		puts "Using impedance calculations..."
		set SARA($instance,"impedanceFlag") 1
	}
	set SARA($instance,"num_terms") $order
	set SARA($instance,"order") 1
	set SARA($instance,"k0") [complex $k0 $k0j]
	set SARA($instance,"k1") [complex $k1 $kj1]
	set SARA($instance,"k2") [complex $k2 $kj2]
	set SARA($instance,"k3") [complex $k3 $kj3]
	set SARA($instance,"k4") [complex $k4 $kj4]
    set SARA($instance,"kn") [complex $kn $kjn]
	set SARA($instance,"sig1") [complex $sig1 $sigj1]
	set SARA($instance,"sig2") [complex $sig2 $sigj2]
	set SARA($instance,"sig3") [complex $sig3 $sigj3]
	set SARA($instance,"sig4") [complex $sig4 $sigj4]

    set solver($instance,"num_iter") 0
	set solver($instance,"tt_Prev") [complex 0 0]
	set solver($instance,"tt_Delta") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Integ") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Prev") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_0") [complex 0 0]
    set solver($instance,"poly_coeffs") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]];     #for current timestep
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
			set [lindex solver($instance,"tt_Delta") 0] $temp
		}
		
		SARA:Accept $instance -2 $v_IC 
		puts "Accept 1"
		SARA:Accept $instance -1 $v_IC 
		puts "Accept 2"
		SARA:Accept $instance 0 $v_IC 
		puts "Accept 3"
	}
	puts "done with Init"
	close $outFile
}

proc SARA:Phi {sigma_i delta_n} {
	# puts "sigma_i: $sigma_i; delta_n: $delta_n"
	# puts "[* $sigma_i $delta_n]"
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
			# return [/ [* [complex 1 0] $delta_n] $c2]
			#return [complex 0 0]
		}
	1	{
			#return $delta_n
			return [/ [* $delta_n [- $c1 [* [+ $c1 $zi] $phi]]] [* $zi $zi]]
			# return [* [/ [* [complex 1 0] $delta_n] $c2] [- $c1 $zi]]
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

# Needs to be called BEFORE the solver object update in "SARA:Accept"
proc SARA:poly_coeffs {instance time input} {
    global solver
    set t0 [complex $time 0];                                       #  most recent time, or t_n
    set t1 [lindex $solver($instance,"tt_Delta") 0];          #  t_(n-1)
	if {[real $t1] == 0 && [real [lindex $solver($instance,"tt_Delta") 1]] == 0 } {
		set t2 [complex -[real $t0] 0]
	} else {
    	set t2 [lindex $solver($instance,"tt_Delta") 1];          #  t_(n-2)
	}
    set x0 [complex $input 0];                                      #  most recent input, or x_n
    set x1 [lindex $solver($instance,"xx_Prev") 0];     #  x_(n-1)
    set x2 [lindex $solver($instance,"xx_Prev") 1];     #  x_(n-2)
	#puts "poly_coeffs"
	#puts "t0: $t0; t1: $t1; t2: $t2"
	#puts "[* [* [complex 2 0] [- $t1 $t2]] [* [- $t0 $t1] [- $t0 $t2]]]"
	#puts "[* [complex 3 0] [- [* [- $t1 $t2] [- $x0 $x1]] [* [- $t0 $t1] [- $x1 $x2]]]]"
    set a $x1
    set b [+ [* [* [- $t1 $t2] [- $t1 $t2]] [- $x0 $x1]] [* [* [- $t0 $t1] [- $t0 $t1]] [- $x1 $x2]]]
    set c [/ [* [complex 3 0] [- [* [- $t1 $t2] [- $x0 $x1]] [* [- $t0 $t1] [- $x1 $x2]]]] [* [* [complex 2 0] [- $t1 $t2]] [* [- $t0 $t1] [- $t0 $t2]]]]
    set d [/ [* [complex [- 0 1] 0] [- [* [- $t1 $t2] [- $x0 $x1]] [* [- $t0 $t1] [- $x1 $x2]]]] [* [* [complex 2 0] [- $t1 $t2]] [* [* [- $t0 $t1] [- $t0 $t1]] [- $t0 $t2]]]]
    set poly_coeffs [list $a $b $c $d]
    #set solver($instance,"poly_coeffs") $poly_coeffs
    return $poly_coeffs
}

# Needs to be called BEFORE the solver object update in "SARA:Accept"
proc SARA:poly_deriv {instance time input} {
    global solver
	if {$time == 0} {
		return [complex 0 0]
	} 
	set poly_coeffs [SARA:poly_coeffs $instance $time $input]
    set delta_t [- [complex $time 0] $solver($instance,"tt_Prev")]
    set powTerms [complex 1.0 0.0]
    set numTerms [llength $solver($instance,"poly_coeffs")]
    set output [complex 0 0]

	for {set i 1} {$i < $numTerms} {incr i} {
		set output [+ $output [* [* [lindex $poly_coeffs $i] $powTerms] [complex $i 0]]]
		set powTerms [* $powTerms $delta_t]
	}

    return $output
}

proc SARA:Step {instance time input} {  
	#puts "step"
    global SARA
    global solver
	global outFile
	set outFile [open "SARA_v8.txt" a+]
	puts $outFile "SARA:Step $instance $time $input"
	close $outFile
	set impedanceFlag $SARA($instance,"impedanceFlag")
    set c_time [complex $time 0]
    set xx_Curr [complex $input 0]
    set delta_n [- $c_time $solver($instance,"tt_Prev")];      # delta_n = current_time - solver(last_recorded_time)
   	set order $SARA($instance,"order")
	set num_terms $SARA($instance,"num_terms")

	if {$impedanceFlag == 0} {
        set xx_Integ_Section [/ [* $delta_n [+ $xx_Curr [lindex $solver($instance,"xx_Prev") 0]]] [complex 2 0]];    # integral from n-1 to n = delta_n * (x[n] + x[n-1]) / 2
        set xx_Integ_Full [+ $xx_Integ_Section [lindex $solver($instance,"xx_Integ") 0]];  	             # integral from 0 to n = xx_Integ (0 to n-1) + xx_Integ_Section (n-1 to n)
    	set final [* $xx_Integ_Full $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)
    	for {set i 1} {$i <= [+ $num_terms 1]} {incr i} {
			set sig_i $SARA($instance,"sig$i")
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			set temp [+ $temp [* $k_i [* $q $xx_Integ_Full]]];        # calculate using most recent x integral, which is not stored
			for {set j 1} {$j <= $order} {incr j} {
				set delta_n [lindex $solver($instance,"tt_Delta") [- $j 1]];
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") [- $j 1]]]]];     # use $j-1 since $j starts from 1
            }
        set final [+ $final $temp]
        }
    } else {
		# puts "step"
		if {[real $SARA($instance,"kn")]==0 &&[imag $SARA($instance,"kn")]==0} {
			set final [* $xx_Curr $SARA($instance,"k0")];
		} else {
    		set final [+ [* $xx_Curr $SARA($instance,"k0")] [* [SARA:poly_deriv $instance $time $input] $SARA($instance,"kn")]];       # Needs to be called BEFORE the solver object update in "SARA:Accept"
		}
		for {set i 1} {$i <= [+ $num_terms 1]} {incr i} {

			set sig_i $SARA($instance,"sig$i")
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			puts "[lindex $solver($instance,"yy") [- $i 1]]"
			puts "$sig_i $delta_n"
			puts "[SARA:Phi $sig_i $delta_n]"
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			set temp [+ $temp [* $k_i [* $q $xx_Curr]]];        # calculate using most recent x value, which is not stored
			for {set j 1} {$j <= $order} {incr j} {
				set delta_n [lindex $solver($instance,"tt_Delta") [- $j 1]];
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Prev") [- $j 1]]]]];     # use $j-1 since $j starts from 1
            }
        set final [+ $final $temp]
	    }
	}
	#puts "Step return: [real $final]"
    return [real $final]

}

proc SARA:Accept {instance time input} {
	#puts "SARA:Accept $instance $time $input"
	#puts "accept"
    global SARA
    global solver
	global outFile
	set outFile [open "SARA_v8.txt" a+]
	puts $outFile "SARA:Accept $instance $time $input"
	close $outFile
    set impedanceFlag $SARA($instance,"impedanceFlag")
	#global outFile
    set c_time [complex $time 0]
    set xx_Curr [complex $input 0]
    
	incr solver($instance,"num_iter")
    set delta_n [- $c_time $solver($instance,"tt_Prev")];     # delta_n = current_time - solver(last_recorded_time)
    set order $SARA($instance,"order")
	set num_terms $SARA($instance,"num_terms")

    if {$impedanceFlag == 0} {
        set xx_Integ_Section [/ [* $delta_n [+ $xx_Curr [lindex $solver($instance,"xx_Prev") 0]]] [complex 2 0]];    # integral from n-1 to n = delta_n * (x[n] + x[n-1]) / 2
        set xx_Integ_Full [+ $xx_Integ_Section [lindex $solver($instance,"xx_Integ") 0]];  	             # integral from 0 to n = xx_Integ (0 to n-1) + xx_Integ_Section (n-1 to n)
    	set final [* $xx_Integ_Full $SARA($instance,"k0")];       # first part of calculation: k0 * integral(x|0 to n)

        set tempList [lrange [linsert $solver($instance,"xx_Integ") 0 $xx_Integ_Full] 0 [+ $order 1]];   # add most recent integral to front of list and shift all elements back
	    set solver($instance,"xx_Integ") $tempList
	    set tempList [lrange [linsert $solver($instance,"tt_Delta") 0 $delta_n] 0 [+ $order 1]];   # add most recent time to front of list and shift all elements back
	    set solver($instance,"tt_Delta") $tempList
        set solver($instance,"tt_Prev") $c_time
        

    	for {set i 1} {$i <= [+ $num_terms 1]} {incr i} {
			set sig_i $SARA($instance,"sig$i")
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			set temp [+ $temp [* $k_i [* $q $xx_Integ_Full]]];        # calculate using most recent x integral, which is not stored
			for {set j 1} {$j <= $order} {incr j} {
				set delta_n [lindex $solver($instance,"tt_Delta") [- $j 1]];  
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Integ") [- $j 1]]]]];     # use $j-1 since $j starts from 1
            }
        set final [+ $final $temp]
        }
    } else {
		# puts "accept"
		if {[real $SARA($instance,"kn")]==0 &&[imag $SARA($instance,"kn")]==0} {
			set final [* $xx_Curr $SARA($instance,"k0")];
		} else {
    		set final [+ [* $xx_Curr $SARA($instance,"k0")] [* [SARA:poly_deriv $instance $time $input] $SARA($instance,"kn")]];       # Needs to be called BEFORE the solver object update in "SARA:Accept"
		}   	
        set tempList [lrange [linsert $solver($instance,"xx_Prev") 0 $xx_Curr] 0 [+ $order 1]];   # add most recent x value to front of list and shift all elements back
	    set solver($instance,"xx_Prev") $tempList
	    set tempList [lrange [linsert $solver($instance,"tt_Delta") 0 $delta_n] 0 [+ $order 1]];   # add most recent time to front of list and shift all elements back
	    set solver($instance,"tt_Delta") $tempList
        set solver($instance,"tt_Prev") $c_time

        for {set i 1} {$i <= [+ $num_terms 1]} {incr i} {
			set sig_i $SARA($instance,"sig$i");
			set k_i $SARA($instance,"k$i");                             # select appropriate k and sigma
			set temp [* [lindex $solver($instance,"yy") [- $i 1]] [SARA:Phi $sig_i $delta_n]];        # lindex using $i-1 because $i starts from 1, not 0
			set q_function "SARA:q$order"
			# set q [$q_function $sig_i $delta_n 0];                            # find the q value for time = n
			# set temp [+ $temp [* $k_i [* $q $xx_Curr]]];        # calculate using most recent x value, which is not stored
			for {set j 0} {$j <= $order} {incr j} {
				set delta_n [lindex $solver($instance,"tt_Delta") $j];
				set q [$q_function $sig_i $delta_n $j];                       # $j starts from 1 since we already used q0 and need to do q1,q2,etc.
				set temp [+ $temp [* $k_i [* $q [lindex $solver($instance,"xx_Prev") $j]]]];     # use $j-1 since $j starts from 1
            }
        	# set final [+ $final $temp]
			set solver($instance,"yy") [lreplace $solver($instance,"yy") [- $i 1] [- $i 1] $temp]
	    }
	}
	#puts "Accept return: [real $final]"
    # set solver($instance,"yy_Final") $final;        # Store final calculated value
}