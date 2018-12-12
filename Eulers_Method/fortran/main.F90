!	Nina Haas
!	Calls Euler's method on function with varying number of steps

program main
	implicit none
	integer i
	real(kind=8) :: t,y,y_0
	call eulers_method(0.0,-2.0,10.0,8.0,'output_8.txt')
	call eulers_method(0.0,-2.0,10.0,16.0,'output_16.txt')
	call eulers_method(0.0,-2.0,10.0,32.0,'output_32.txt')
	call eulers_method(0.0,-2.0,10.0,64.0,'output_64.txt')
end program main	
