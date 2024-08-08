proc SARA:Init {instance args} {
    global SARA
    
    # ...
    
    set SARA($instance,order) 3
    set SARA($instance,p0) 1.2
    set SARA($instance,p1) -2
    set SARA($instance,p2) 1.255
    set SARA($instance,p3) 1

    set SARA($instance,q0) 1.2
    set SARA($instance,q1) -2
    set SARA($instance,q2) 1.255
    set SARA($instance,q3) 1
}

proc SARA:Step {instance time input} {
    global SARA
    
    # ...
    
    set output 5.2
    return $output
}


proc SARA:Accept {instance time input output} {
    global SARA
    
    # ...
    
    set SARA($instance,internal_data) 1.2
}
