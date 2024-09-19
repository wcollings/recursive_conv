# How this is going to work is, we'll have a "step" function
# that computes the next value, but doesn't save it anywhere.
# instead, it just returns the result and moves on. The solver
# is going to try a bunch of different values to see which one
# converges mathematically the best, or is the most stable, etc.
# So our framework needs to not save those intermediary results
# instead, there's a function at the bottom called "accept"
# that actually saves that result.

package require math::complexnumbers
interp alias {} complex {} math::complexnumbers::complex
interp alias {} ccos {} math::complexnumbers::cos
interp alias {} cexp {} math::complexnumbers::exp
interp alias {} real {} math::complexnumbers::real
interp alias {} imag {} math::complexnumbers::imag
interp alias {} c+ {} math::complexnumbers::+
interp alias {} c- {} math::complexnumbers::-
interp alias {} c* {} math::complexnumbers::*
interp alias {} c/ {} math::complexnumbers::/
interp alias {} cpow {} math::complexnumbers::pow


proc fac {n} {
	if { $n <= 1} {
		return 1
	}
	return [tcl::mathop::* $n [fac [tcl::mathop::- $n 1]]]
}
proc sign {n} {
	if {$n >= 0} {
		return 1
	} else {
		return -1
	}
}
proc r_pow {a b} {
	if {$b == 0 } {
		return 1
	} elseif {$b==1 } {
		return $a
	} else {
		return [tcl::mathop::* $a [r_pow $a [tcl::mathop::- $b 1]]]
	}
}
proc expm1 {zeta} {
	set num_terms 6
	set result 0
	set sgn [sign $zeta]
	for {set i 1} {$i <= $num_terms } {incr i} {
		set p [expr {$zeta**$i}]
		set f [fac $i]
		set result [expr {$p/$f}]
	}
	return $result
}

proc c_exp { zeta } {
	set r [expm1 [real $zeta]]
	set i [expm1 [imag $zeta]]
	return [complex $r $i]
}

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
			set order [tcl::mathop::- $i 1]
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
	cexp [c* $i $n]
	# expr {exp($i*$n)}
}

proc SARA:zeta {i n} {
	c- [c* $i $n]
	# expr {-$i*$n}
}

proc SARA:q0 {sigma_i delta_n n} {
	if {[real $delta_n] == 0 } {
		return $delta_n
	}
	return [* [/ $delta_n [SARA:zeta $sigma_i $delta_n]] [- [complex 1 0] [SARA:Phi $sigma_i $delta_n]]]
	# expr {[$delta_n/[SARA:zeta $sigma_i $delta_n]]*[1-[SARA:Phi $sigma_i $delta_n]]}
}

proc SARA:q1 {sigma_i delta_n n} {
	if {[real $delta_n] == 0} {
		return $delta_n;
	}
	# set phi [SARA:Phi $sigma_i $delta_n]
	set zi [SARA:zeta $sigma_i $delta_n]
	set c1 [complex 1 0]
	set c1n [complex -1 0] 
	set c2 [complex 2 0]
	switch $n {
	0 	{
			return [c* [c/ $delta_n [c* $zi $zi]] [c- [c_exp [c* $sigma_i $delta_n]] $zi]]
			# return [* [/ $delta_n [pow $zi $c2]] [+ $c1n [+ $zi $phi]]]
		# return [expr {$delta_n/[pow $zi 2]*[-1+$zi+$phi]}]			;# supposed to be cpow
		}
	1	{
			return $delta_n
			# return [* [/ $delta_n [pow $zi $c2]] [* [- $c1 [+ $c1 $zi]] $phi]]
		# return [expr {$delta_n/[pow $zi 2]*[1-[1+$zi]*$phi]}] 		;# supposed to be cpow
		}
	}
}

proc SARA:q2 {sigma_i delta_n n} {
	if {$delta_n==[complex 0 0]} {
		return $delta_n
	}
	set phi [expr {Phi $sigma_i $delta_n}]
	set zi [expr {zeta $sigma_i $delta_n}]
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

proc SARA:Step {instance time input}   {   
	global SARA
	global solver
	set time [complex $time 0]
	set input [complex $input 0]

	set delta_n [- $time  $solver($instance,"prev_time")]
	set integ [* $delta_n $input]
	set currx [+ [lindex $solver($instance,"xx") 0] $integ ]

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
proc abs {x} {expr {abs($x)}}


proc SARA:Accept {instance time input} {
	global fp
	global SARA
	global solver
	set time [complex $time 0]
	set input [complex $input 0]
	set order $SARA($instance,"order")
	incr solver($instance,"num_iter")

	set delta_n [c- $time $solver($instance,"prev_time")]
	set integ [c* $delta_n $input]
	set currx [c+ [lindex $solver($instance,"xx") 0] $integ ]
	# Enqueue
	set tempList [lrange [linsert $solver($instance,"xx") 0 $currx] 0 [expr {$order+1}]]
	set solver($instance,"xx") $tempList
	set tempList [lrange [linsert $solver($instance,"tt") 0 $delta_n] 0 [expr {$order+1}]]
	set solver($instance,"tt") $tempList
	set solver($instance,"prev_time") $time
	puts -nonewline $fp "[real $currx],"

	set final [c* $SARA($instance,"k0") $currx]
	
	for {set i 0} {$i <= $order} {incr i} {
		set sigma $SARA($instance,"q$i")
		set a $SARA($instance,"p$i")
		set t1 [lindex $solver($instance,"yy") $i]
		set t2 [SARA:Phi $sigma $delta_n]
		set temp [c* [lindex $solver($instance,"yy") $i] [SARA:Phi $sigma $delta_n]]
		set qq "SARA:q$order"
		set q [$qq $sigma $delta_n 0]
		for {set j 0} {$j <= $order} {incr j} {
			set q [$qq $sigma $delta_n $j]
			set temp [c+ $temp [c* $a [c* $q [lindex $solver($instance,"xx") $j]]]]
		}
		set solver($instance,"yy") [lreplace $solver($instance,"yy") $i $i $temp]
		puts -nonewline $fp "[real $temp],"
		set final [c+ $final $temp]
	}
	set solver($instance,"xx_Prev") $input
	set solver($instance,"yy_Final") $final
	puts $fp [real $final]
}

SARA:Init "inst" 279695688.23288864 0 -740340728373193.2 152205548356904.16 -740340728373193.2 -152205548356904.16 0 0 0 0 -2923683.340113021 3792673.1864589 -2923683.340113021 -3792673.1864589 0 0 0 0 
set fp [open "output.csv" w]
puts $fp "t,v,dv,y0,y1,out"
set pi 3.14
set freq 1e6
set tstep [expr {1/$freq/20}]
set tend [expr {$tstep * 999}]
for {set t 0} {$t < $tend } {set t [expr {$t+$tstep}]} {
	# set t [expr {$i+$step}]
	# set t [expr {$i*$step}]
	set tau [expr {2*$freq*$t*$pi}]
	set v [real [ccos [complex $tau 0]]]
	# set i [real [SARA:Accept "inst" $tau $v]]
	puts -nonewline $fp "$t,$v,"
	SARA:Accept "inst" $t $v
}
