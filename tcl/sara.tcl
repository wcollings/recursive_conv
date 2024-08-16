# How this is going to work is, we'll have a "step" function
# that computes the next value, but doesn't save it anywhere.
# instead, it just returns the result and moves on. The solver
# is going to try a bunch of different values to see which one
# converges mathematically the best, or is the most stable, etc.
# So our framework needs to not save those intermediary results
# instead, there's a function at the bottom called "accept"
# that actually saves that result.

# TODO: translate Step, Accept, and the q functions (+zeta and phi)

interp alias {} ? {} set errorInfo
package require math::complexnumbers
namespace import ::math::complexnumbers::*

proc SARA:Init {instance k0 k0j p0 p0j p1 p1j p2 p2j p3 p3j q0 q0j q1 q1j q2 q2j q3 q3j} {
	global SARA
	global solver
	set i 1
	if {$p0 == 0 && $p0j == 0} {
		puts "p0 is 0, not allowed"
	}
	set order 3
	while {$i < 4} {
		set realpart_name q$i
		set imagpart_name [join [list "q" $i "j"] ""]
		upvar 0 $realpart_name realpart
		upvar 0 $imagpart_name imagpart
		if {$realpart == 0 && $imagpart == 0} {
			set order [- $i 1]
			break
		}
		incr i
	}
	puts "Found order: $order"

	set SARA($instance,"order") $order
	set SARA($instance,"k0") [complex $k0 $k0j]
	set SARA($instance,"p0") [complex $p0 $p0j]
	set SARA($instance,"p1") [complex $p1 $p1j]
	set SARA($instance,"p2") [complex $p2 $p2j]
	set SARA($instance,"p3") [complex $p3 $p3j]
	set SARA($instance,"q0") [complex $q0 $q0j]
	set SARA($instance,"q1") [complex $q1 $q1j]
	set SARA($instance,"q2") [complex $q2 $q2j]
	set SARA($instance,"q3") [complex $q3 $q3j]
	
	set solver($instance,"num_iter") 0
	set solver($instance,"prev_time") [complex 0 0]
	set solver($instance,"tt") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Prev") [complex 0 0]
	set solver($instance,"yy") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"yy_Final") [complex 0 0]
}

proc SARA:Phi {i n} {
	exp [* $i $n]
	# expr {exp($i*$n)}
}

proc SARA:zeta {i n} {
	- [* $i $n]
	# expr {-$i*$n}
}

proc SARA:q0 {sigma_i delta_n n} {
	expr {[$delta_n/[zeta $sigma_i $delta_n]]*[1-[Phi $sigma_i $delta_n]]}
}

proc SARA:q1 {sigma_i delta_n n} {
	set phi [SARA:Phi $sigma_i $delta_n]
	set zi [SARA:zeta $sigma_i $delta_n]
	if {$delta_n==0} {
		return 0;
	} else {
		set c1 [complex 1 0]
		set c1n [complex -1 0] 
		set c2 [complex 2 0]
		switch $n {
		0 	{
				return [* [/ $delta_n [pow $zi $c2]] [+ $c1n [+ $zi $phi]]]
			# return [expr {$delta_n/[pow $zi 2]*[-1+$zi+$phi]}]			;# supposed to be cpow
			}
		1	{
				return [* [/ $delta_n [pow $zi $c2]] [* [- $c1 [+ $c1 $zi]] $phi]]
			# return [expr {$delta_n/[pow $zi 2]*[1-[1+$zi]*$phi]}] 		;# supposed to be cpow
			}
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

proc SARA:Step {instance time input}   {   
	global SARA
	global solver
	set time [complex $time 0]
	set input [complex $input 0]

	set delta_n [- $time  [lindex $solver($instance,"tt") 0]]
	set currx [+ [lindex $solver($instance,"xx") 0] [* $delta_n [- $input $solver($instance,"xx_Prev")]]]

	set order $SARA($instance,"order")
	set final [* $SARA($instance,"k0") $currx]
	
	for {set i 0} {$i < $order} {incr i} {
		set sigma $SARA($instance,"q$i")
		set a $SARA($instance,"p$i")
		set temp [* [lindex $solver($instance,"yy") $i] [SARA:Phi $sigma $delta_n]]
		set qq "SARA:q$order"
		
		set q [$qq $sigma $delta_n 0]
		set temp [+ $temp [* $a [* $q $currx]]]
		for {set j 0} {$j < [- $order 1]} {incr j} {
			set q [qq $sigma $delta_n [expr{$j+1}]]
			set temp [+ $temp [* $a [* $q [lindex $solver($instance,"xx") $j]]]]
		}
		set $solver($instance,"yy") [lreplace $solver($instance,"yy") $i $i $temp]
		set final [+ $final $temp]
	}
	return [real $final]
}

proc SARA:Accept {instance time input} {
	global SARA
	global solver
	set time [complex $time 0]
	set input [complex $input 0]

	set delta_n [- $time  [lindex $solver($instance,"tt") 0]]
	set currx [+ [lindex $solver($instance,"xx") 0] [* $delta_n [- $input $solver($instance,"xx_Prev")]]]

	set tempList [list 0 0 0]
	set order $SARA($instance,"order")
	set final [* $SARA($instance,"k0") $currx]
	
	for {set i 0} {$i <= $order} {incr i} {
		set sigma $SARA($instance,"q$i")
		set a $SARA($instance,"p$i")
		set temp [* [lindex $solver($instance,"yy") $i] [SARA:Phi $sigma $delta_n]]
		set qq "SARA:q$order"
		
		set q [$qq $sigma $delta_n 0]
		set temp [+ $temp [* $a [* $q $currx]]]
		for {set j 0} {$j < [- $order 1]} {incr j} {
			set q [qq $sigma $delta_n [expr{$j+1}]]
			set temp [+ $temp [* $a [* $q [lindex $solver($instance,"xx") $j]]]]
		}
		set solver($instance,"yy") [lreplace $solver($instance,"yy") $i $i $temp]
		set final [+ $final $temp]
	}
	set tempList [lrange $solver($instance,"xx") 0 2]
	set solver($instance,xx) [concat xx_Curr tempList]
	set tempList [lrange $solver($instance,"tt") 0 2]
	set solver($instance,tt) [concat time tempList]
	set solver($instance,xx_Prev) $input
	set solver($instance,yy_Final) $final
}

SARA:Init "inst" 279695688.23288864 0 -740340728373193.2 152205548356904.16 -740340728373193.2 -152205548356904.16 0 0 0 0 -2923683.340113021 3792673.1864589 -2923683.340113021 -3792673.1864589 0 0 0 0 
SARA:Accept "inst" 1e-8 2
puts "[SARA:Step "inst" 2e-8 2]"
parray solver
