/* How this is going to work is, we'll have a "step" function
 * that computes the next value, but doesn't save it anywhere.
 * instead, it just returns the result and moves on. The solver
 * is going to try a bunch of different values to see which one
 * converges mathematically the best, or is the most stable, etc.
 * So our framework needs to not save those intermediary results
 * instead, there's a function at the bottom called "accept"
 * that actually saves that result.
*/

//TODO: figure out how to initialize arrays, including multidimentional ones
//TODO: what args does Init need? It may be a fair number, for instance
//      the entirety of the pade approximant coefficients and the order
//TODO: translate Step, Accept, and the q functions (+zeta and phi)
//TODO: do we need to do something different if the coefficients are
//      complex numbers? How do we deal with that?
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

	 set solver(num_iter) 0
	 set solver(prev_time) 0
	 set solver(tt) {}   //this should be a list
	 set solver(xx) {}   //this should be a list
	 set solver(yy) {}   //this should be a multidimentional list
}

proc SARA:Phi {i n} {
	expr {exp($i*$n)}
}

proc SARA:zeta {i n} {
	expr {-$i*$n}
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
		 set temp [expr {$solver(yy,i,0)*[SARA:Phi $sigma delta_n]}]
		 for {set j 0} {$j < $order} {incr j} {
			 set q [qq sigma delta_n]
			 set temp [expr {$temp + $a*$q*$solver(xx,$j)}]
		 }
		 set final [expr {$final + $temp}]
	 }
    return $final
}


proc SARA:Accept {instance time input output} {
    global SARA
    
    # ...
    
    set SARA($instance,internal_data) 1.2
}
