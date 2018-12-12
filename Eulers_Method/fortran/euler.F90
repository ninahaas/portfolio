!	Nina Haas
!	Includes specific ODE and subroutine that uses Euler's
!	method to approximate coordinate points and prints to output file

!	dydt()
!	Input: time, value
!	Output: value of differential equation
real(kind=8) function dydt(t,y)
	implicit none
	real(kind=8), intent(in) :: t,y
	dydt = 2.0*t/(y*(1.0+t**2.0))
end function dydt

!	eulers_method()
!	Input: beginning time, initial value, end time,
!	number of data points, output file name
!	Output: None
subroutine eulers_method(t_0,y_0,t_f,N,file_name)
	implicit none
	real(kind=8), external :: dydt
	real(kind=8), intent(in) :: y_0,t_0,t_f,N
	character(len=*), intent(in) :: file_name
	real(kind=8) :: y,h,t
	
	h = (t_f - t_0)/(N - 1)
	y = y_0	
	t = t_0
	open(unit=20, file=file_name)

	do while(t<=t_f)
		write(20,*) t, y
		y = y + h*dydt(t,y)
		t = t + h
	end do

	close(20)
end subroutine eulers_method
				
