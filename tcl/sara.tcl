#SARA_v1
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
	# global outFile
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
	set scalingK [complex 1 0]
	set scalingSigma [complex 1 0]
	set SARA($instance,"order") $order
	set SARA($instance,"k0") [complex $k0 $k0j]
	set SARA($instance,"p0") [* $scalingK [complex $p0 $p0j]]
	set SARA($instance,"p1") [* $scalingK [complex $p1 $p1j]]
	set SARA($instance,"p2") [* $scalingK [complex $p2 $p2j]]
	set SARA($instance,"p3") [* $scalingK [complex $p3 $p3j]]
	set SARA($instance,"q0") [* $scalingSigma [complex $q0 $q0j]]
	set SARA($instance,"q1") [* $scalingSigma [complex $q1 $q1j]]
	set SARA($instance,"q2") [* $scalingSigma [complex $q2 $q2j]]
	set SARA($instance,"q3") [* $scalingSigma [complex $q3 $q3j]]
	
	set solver($instance,"num_iter") 0
	set solver($instance,"prev_time") [complex 0 0]
	set solver($instance,"tt") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"xx_Prev") [complex 0 0]
	set solver($instance,"yy") [list [complex 0 0] [complex 0 0] [complex 0 0] [complex 0 0]]
	set solver($instance,"yy_Final") [complex 0 0]

	# set outFile [open "integral_values.csv" w+]
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
	if {[real $delta_n] == 0 } {
		return $delta_n
	}
	return [* [/ $delta_n [SARA:zeta $sigma_i $delta_n]] [- [complex 1 0] [SARA:Phi $sigma_i $delta_n]]]
	# return $delta_n
	# expr {[$delta_n/[SARA:zeta $sigma_i $delta_n]]*[1-[SARA:Phi $sigma_i $delta_n]]}
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
			# return [* $c1n [/ $delta_n $c2]]
			# return [* [/ $delta_n [* $zi $zi]] [- [c_exp [* $sigma_i $delta_n]] $zi]]
			# return [* $delta_n [complex -0.07685 5e-11]]
			# return [complex -2.68e-11 5.15e-11]
			return [complex 0 0]
			# return [/ $c1 [* [pow $sigma_i $c2] $delta_n]]
		 	# return [* [/ $delta_n [pow $zi $c2]] [+ $c1n [+ $zi $phi]]]			;# supposed to be cpow
		 	
		}
	1	{
			# return [* [/ $delta_n $c2] [+ $c3 $zi]]
			return $delta_n
			# return [* [/ $delta_n [pow $zi $c2]] [- $c1 [* [+ $c1 $zi] $phi]]]
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

	set delta_n [- $time $solver($instance,"prev_time")]
	# set integ [* $delta_n $input]
	# set integ [* $delta_n $solver($instance,"xx_Prev")]
	set integ [* [/ $delta_n [complex 2 0]] [+ $input $solver($instance,"xx_Prev")]]
	set currx [+ [lindex $solver($instance,"xx") 0] $integ ]

	set order $SARA($instance,"order")
	set final [* $SARA($instance,"k0") $currx]
	
	for {set i 0} {$i <= $order} {incr i} {
		set sigma $SARA($instance,"q$i")
		set a $SARA($instance,"p$i")
		set temp [* [lindex $solver($instance,"yy") $i] [SARA:Phi $sigma $delta_n]]
		set qq "SARA:q$order"
		set q [$qq $sigma $delta_n 0]
		set temp [+ $temp [* $a [* $q $currx]]]
		for {set j 0} {$j < $order} {incr j} {
			set q [$qq $sigma $delta_n [+ $j 1]]
			set temp [+ $temp [* $a [* $q [lindex $solver($instance,"xx") $j]]]]
		}
		# set $solver($instance,"yy") [lreplace $solver($instance,"yy") $i $i $temp]
		set final [+ $final $temp]
	}
	return [real $final]
}
proc abs {x} {expr {abs($x)}}


proc SARA:Accept {instance time input} {
	global fp
	global SARA
	global solver
	# global outFile
	set time [complex $time 0]
	set input [complex $input 0]
	set order $SARA($instance,"order")
	incr solver($instance,"num_iter")
	set prev_time $solver($instance,"prev_time")

	set delta_n [- $time $prev_time]
	# set integ [* $delta_n $input]
	# set integ [* $delta_n $solver($instance,"xx_Prev")]
	set integ [* [/ $delta_n [complex 2 0]] [+ $input $solver($instance,"xx_Prev")]]
	set currx [+ [lindex $solver($instance,"xx") 0] $integ ]
	# Enqueue

	set tempList [lrange [linsert $solver($instance,"xx") 0 $currx] 0 [expr {$order + 1}]]
	set solver($instance,"xx") $tempList
	set tempList [lrange [linsert $solver($instance,"tt") 0 $delta_n] 0 [expr {$order + 1}]]
	set solver($instance,"tt") $tempList
	set solver($instance,"prev_time") $time

	set final [* $SARA($instance,"k0") $currx]
	
	for {set i 0} {$i <= $order} {incr i} {
		set sigma $SARA($instance,"q$i")
		set a $SARA($instance,"p$i")
		set t1 [lindex $solver($instance,"yy") $i]
		set t2 [SARA:Phi $sigma $delta_n]
		set temp [* [lindex $solver($instance,"yy") $i] [SARA:Phi $sigma $delta_n]]
		set qq "SARA:q$order"
		# set q [$qq $sigma $delta_n 0]
		# puts $outFile "q$i, $q"
		for {set j 0} {$j <= $order} {incr j} {
			set q [$qq $sigma $delta_n $j]
			# puts $outFile "q$i$j, $q"
			set temp [+ $temp [* $a [* $q [lindex $solver($instance,"xx") $j]]]]
		}
		set solver($instance,"yy") [lreplace $solver($instance,"yy") $i $i $temp]
		set final [+ $final $temp]
	}
	set solver($instance,"xx_Prev") $input
	set solver($instance,"yy_Final") $final
}
proc callAccept {instance time input} {
	if {[catch {SARA:Accept $instance $time $input } errmsg]} {
		puts "ErrorMsg: $errmsg"
		# puts "ErrorCode: $errorCode"
		# puts "ErrorInfo:\n$errorInfo\n"
	}
}

proc callStep {instance time input} {
	if {[catch {SARA:Step $instance $time $input } errmsg]} {
		puts "ErrorMsg: $errmsg"
		# puts "ErrorCode: $errorCode"
		# puts "ErrorInfo:\n$errorInfo\n"
	}
}

# SARA:Init "inst" 0 0 947766.7355944953 0 52233.26440550487 0 0 0 0 0 -2104987.5621120883 0 -95012.43788791099 0 0 0 0 0 
# SARA:Init "inst" 279695688.23288864 0 -740340728373193.2 152205548356904.16 -740340728373193.2 -152205548356904.16 0 0 0 0 -2923683.340113021 3792673.1864589 -2923683.340113021 -3792673.1864589 0 0 0 0 
# set pi 3.14
# set step [expr {$pi/1e3}]
# for {set t 0} {$t < $pi} {set t [expr {$t+$step}]} {
	# set t [expr {$i+$step}]
	# set t [expr {$i*$step}]
	# set v [real [cos [complex $t 0]]]
	# real [SARA:Accept "inst" $t $v]
	# puts $fp "$t,$v,[real [SARA:Accept "inst" $t $v]]"
# }
