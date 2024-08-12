# How this is going to work is, we'll have a "step" function
# that computes the next value, but doesn't save it anywhere.
# instead, it just returns the result and moves on. The solver
# is going to try a bunch of different values to see which one
# converges mathematically the best, or is the most stable, etc.
# So our framework needs to not save those intermediary results
# instead, there's a function at the bottom called "accept"
# that actually saves that result.

# TODO: figure out how to initialize arrays, including multidimentional ones
# TODO: what args does Init need? It may be a fair number, for instance
#      the entirety of the pade approximant coefficients and the order
# TODO: translate Step, Accept, and the q functions (+zeta and phi)
# TODO: do we need to do something different if the coefficients are
#      complex numbers? How do we deal with that?

package require qcomplex

proc SARA:Init {instance args} {
    global SARA
	 global solver
    
    set SARA($instance,order) 3
    set SARA($instance,p0) 1.2
    set SARA($instance,p1) -2
    set SARA($instance,p2) 1.255
    set SARA($instance,p3) 1

    set SARA($instance,q0) 1.2
    set SARA($instance,q1) -2
    set SARA($instance,q2) 1.255
    set SARA($instance,q3) 1
set solver("num_iter") 0
	 set solver("prev_time") 0
	 set solver("tt") {{0 0 0 0}}   
	 set solver("xx") {{0 0 0 0}}
	 set solver("yy","out") [list 0 0 0 0 0 0 0 0 0]
	 set solver("yy","out") [list 0 0 0 0 0 0 0 0 0]
	 for {set i 0} {$i < $order} {incr i} {
		 lappend solver("yy")
	 }
	 set solver(yy) {{0 0 0 0}}
}

proc SARA:Phi {i n} {
	expr {exp($i*$n)}
}

proc SARA:zeta {i n} {
	expr {-$i*$n}
}

proc SARA:q0 {sigma_i delta_n n} {
	expr {[$delta_n/[zeta $sigma_i $delta_n]]*[1-[Phi $sigma_i $delta_n]]}
}

proc SARA:q1 {sigma_i delta_n n} {
	set phi [expr {Phi $sigma_i $delta_n}]
	set zi [expr {zeta $sigma_i $delta_n}]
	if {$delta_n==0} {
		return 0;
	}
	else switch $n {
		0 	{
			return [expr {$delta_n/[qcomplex::pow $zi 2]*[-1+$zi+$phi]}]			;# supposed to be cpow
			}
		1	{
			return [expr {$delta_n/[qcomplex::pow $zi 2]*[1-[1+$zi]*$phi]}] 		;# supposed to be cpow
			}
	}
}

proc SARA:q2 {sigma_i delta_n n} {
	set phi [expr {Phi $sigma_i $delta_n}]
	set zi [expr {zeta $sigma_i $delta_n}]
	if {$zi==0} {
		return 0;
	}
	else switch $n {
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


proc SARA:Step {instance t x} {
    global SARA
	 global solver

	 set order $SARA($instance,order)
	 set final 0
	 set delta_n [expr {t-$solver(prev_time)}]
	 for {set i 0} {$i < $order} {incr i} {
		 set sigma $SARA($instance,"q$i")
		 set a $SARA($instance, "p$i")
		 set temp [expr {[lindex $solver(yy) i 0]*[SARA:Phi $sigma delta_n]}]
		 for {set j 0} {$j < $order} {incr j} {
			 set q [qq sigma delta_n]
			 set temp [expr {$temp + $a*$q*[lindex $solver("xx") $j}]
		 }
		 incr final $temp
	 }
    return $final
}

# TODO: make this save the results
proc SARA:Accept {instance time input output} {
    global SARA
    global SARA
	 global solver

	 set order $SARA($instance,order)
	 set final 0
	 set delta_n [expr {t-$solver(prev_time)}]
	 for {set i 0} {$i < $order} {incr i} {
		 set sigma $SARA($instance,"q$i")
		 set a $SARA($instance, "p$i")
		 set temp [expr {[lindex $solver(yy) i 0]*[SARA:Phi $sigma delta_n]}]
		 for {set j 0} {$j < $order} {incr j} {
			 set q [qq sigma delta_n]
			 set temp [expr {$temp + $a*$q*[lindex $solver("xx") $j}]
		 }
		 incr final $temp
	 }

}
