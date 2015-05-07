!!!!!!!!!!!!! 3D GRADIENT MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Author : Bryan Gorman
!!!!!!!
!!!!!!!          u(y,z)   d C(x,y,z)    d**2 C   d**2 C   d**2 C
!!!!!!! PDE  Pe ------   ---------- = ------ + ----- + -----
!!!!!!!            U        d x         d x**2   d y**2   d z**2
!!!!!!!
!!!!!!!
!!!!!!! BCs:   C(0,y,z) = f(x)
!!!!!!!        dC/dx (inf,y,z) = 0
!!!!!!!        (expressed in this model in an 'outflow' BC)
!!!!!!!        dC/dy (x,0,z) = 0
!!!!!!!        dC/dy (x,1,z) = 0
!!!!!!!        dC/dz (x,y,0) = 0
!!!!!!!        dC/dz (x,y,e) = 0
!!!!!!!
!!!!!!!               U W          H
!!!!!!! where Pe = ---------, e = ---
!!!!!!!                D           W

!!! Conventions: 'East' (E) is next index value in x direction
!!!              'West' (W) is previous index value in x
!!!              'North' (N) is next index value in y
!!!              'South' (S) is previous index value in y
!!!              'Top' (T) is next index value in z
!!!              'Bottom' (B) is previous index value in z
!!!
!!!               i is used as x index
!!!               j is used as y index
!!!               k is used as z index


program model
implicit none

character(LEN=*), parameter :: &
path = 'vel_test400x40_test2_Pe5000\'

!!!!!!!!! INITIALIZE VARIABLES AND PARAMETERS !!!!!!!!!!!!!

integer, parameter :: &

!!!!!!!! Mesh refinement !!!!!!!!!!!!!!!!!!!!!!!!!!!
Nx = 1000, &       ! x mesh divisions
Ny = 40, &       ! y mesh divisions
Nz = 40          ! z mesh divisions
! Note that due to mesh generating scheme, Ny and Nz (but not Nx)
! must be even

real*8, parameter :: &

residual_drop = 1e-6, & ! Magnitude of decrease in overall residual before termination of iterations
residual_drop_2D = 1e-2, & ! Magnitude of decrease in plane residual before next plane is called

mifx = 1, &  ! MIFx == Mesh Increase Factor in x
mify = 1+3.0D0/Ny, &!  + 1.0D0/Ny, &  ! MIFy == Mesh Increase Factor in y
mifz = 1, &! + 1.0D0/Nz, &  ! MIFz == Mesh Increase Factor in z

pi = 3.14159265358979323846264338327950288419716939937510582097494, &

alpha = 1, &     ! Relaxation factor where alpha > 1 = overrelaxation,
                ! and alpha < 1 = underrelaxation.
                ! alpha > 1 will probably make convergence unstable

! !!!!!!!!! User input of physical parameters !!!!!!!!!!!!
! Q_in = 1, &      ! microliter/min (volumetric flow rate)
! W_in = 500, &    ! micron (width)
! H_in = 100, &    ! micron (height)
! D_in = 8e-7, &   ! cm^2/s  (diffusion coeff.)
! L_in = 9000, &   ! micron (length of mixing channel)

! !!!!! Conversion of physical parameters to SI units !!!!!
! Q_conv = 1.6667e-11*Q_in, &   !m^3/s
! W_conv = 1e-6*W_in, &         !m
! H_conv = 1e-6*H_in, &         !m
! D_conv = 1e-4*D_in, &         !m^2/s
! L_conv = 1e-6*L_in, &         !m
! U_conv = Q_conv/(W_conv*H_conv), &           !m/s  ! Calculation of average velocity

! !!!!!!!!!! Important parameters !!!!!!!!!!!!!!!!!!!!!!!
! Pe = U_conv*W_conv/D_conv, &         ! Peclet number
! e = H_conv/W_conv, &            ! epsilon = H/W (aka aspect ratio)

Pe = 5000, &
e = 0.1, &

!!!!!!!!!! Taylor Dispersion Constants !!!!!!!!!!!!!!!!!
k_tay = 1 + Pe**2*e**2/210, & ! Taylor Dispersion Constant
k_enh = 1 + 7.951*Pe**2*e**2/210, & ! 'Enhanced' Taylor Dispersion Constant

!!!!!!!!!! Dimensionless lengths !!!!!!!!!!!!!!!!!!!!!!!
! l = L_conv/W_conv, &  ! Dimensionless length scaled by W
l = .01*Pe, &
w = 1, &    ! Dimensionless width scaled by W
h = e, &    ! Dimensionless height scaled by W

! Input Boundary Condition for Analytic model
! C1 = 0, &
! C2 = 1

C1 = 0, &
C2 = 0.25, &
C3 = 0.5, &
C4 = 0.75, &
C5 = 1

integer :: &
i, & ! x index
j, & ! y index
k, & ! z index
n, & ! velocity counter
m, & ! analytical model counter
ni, & ! sweeping iteration counter
nj ! plane iteration counter

real*8 :: &
u_ave = 1, &!u_ave = 0,&
u_ave_old = 0, &
u_old = 0, &
!Baselength variables
baselength_x = 0, &
baselength_y = 0, &
baselength_z = 0, &
!Summation terms
sum_x = 0, &
sum_y = 0, &
sum_z = 0, &
sum_u_ave = 0, &
sum_u = 0, &
sum_ave_y = 0, &
sum_ave_z = 0, &
sum_term = 0, &
sum_model = 0, &
sum_model_old = 0, &
sum_model_nodiff = 0, &
sum_model_taylor = 0, &
! sum_model_taylor_enh = 0, &
mean_tot = 0, &
Am = 0, &
Ao = 1/2, &
!Error terms
error = 0, &
error_2D = 0, &
error_orig = 0, &
errorsum = 0, &
errorsum_2D = 0,&
r_2D = 0, &
r = 0, &
!Local diffusion and convection terms
D_e = 0, &
D_w = 0, &
D_n = 0, &
D_s = 0, &
D_t = 0, &
D_b = 0, &
P_e = 0, & 
P_w = 0, &
F_e = 0, &
F_w = 0

real*8, dimension(1:Nx,1:Ny,1:Nz) :: &
a_P = 0, &
a_E = 0, &
a_W = 0, &
a_N = 0, &
a_S = 0, &
a_T = 0, &
a_B = 0, &
C = 0

real*8, dimension(1:Nx,1:Ny) :: &
C_model = 0, &
C_model_nodiff = 0, &
C_model_taylor = 0, &
! C_model_taylor_enh = 0, &
C_ave_z = 0

real*8, dimension(1:Nx,1:Nz) :: &
C_ave_y = 0

real*8, dimension(1:Ny,1:Nz) :: &
u = 0

real*8, dimension(1:Nx) :: &
x = 0, &
deltax = 0

real*8, dimension(1:Ny) :: &
y = 0, &
deltay = 0

real*8, dimension(1:Nz) :: &
z = 0, &
deltaz = 0, &
b_P = 0, &
a = 0, &
b = 0, &
c_it = 0, &
d = 0, &
P = 0, &
Q = 0


print *, mifx
print *, mify
print *, mifz

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'k_tay')
WRITE (11,*) k_tay
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'k_enh')
WRITE (11,*) k_enh
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'epsilon')
WRITE (11,*) e
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'Pe')
WRITE (11,*) Pe
CLOSE (11)

!!!!!!!!!!!!!! Generate x mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!

! sum_x = 0
! do i = 2,Nx-1
! 	sum_x = sum_x + mifx**(i-2)
! end do
! baselength_x = l/(sum_x)

! deltax(1) = 0
! do i = 2,Nx-1
! 	deltax(i) = mifx**(i-2) * baselength_x
! end do
! deltax(Nx) = 0

! do i = 2,Nx
! 	x(i) = x(i-1)+deltax(i-1)/2 + deltax(i)/2
! end do

! write(6,*) "x mesh defined"

deltax = l/Nx
do i = 2,Nx
x(i) = x(i-1) + deltax(i)
end do


!!!!!!!!!!!!!!! Generate y mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!


sum_y = 0
do j = 2,Ny/2
	sum_y = sum_y+mify**(j-2)
end do
baselength_y = w/(2*sum_y)

deltay(1) = 0
do j = 2,Ny/2
	deltay(j) = mify**(j-2) * baselength_y
	deltay(Ny-j+1) = mify**(j-2) * baselength_y
end do
deltay(Ny) = 0

do j = 2,Ny
	y(j) = y(j-1)+deltay(j-1)/2 + deltay(j)/2
end do

write(6,*) "y mesh defined"

!!!!!!!!!!!!!!! Generate z mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!

sum_z = 0
do k = 2,Nz/2
    sum_z = sum_z+mifz**(k-2)
end do
baselength_z = h/(2*sum_z)

deltaz(1) = 0
do k = 2,Nz/2
	deltaz(k) = mifz**(k-2) * baselength_z
	deltaz(Nz-k+1) = mifz**(k-2) * baselength_z
end do
deltaz(Nz) = 0

do k = 2,Nz
	z(k) = z(k-1)+deltaz(k-1)/2 + deltaz(k)/2
end do

write(6,*) "z mesh defined"

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'x')
WRITE (11,*) x
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'y')
WRITE (11,*) y
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'z')
WRITE (11,*) z
CLOSE (11)

!!! Compute analytically-predicted concentration profile,
!!! both with and without Taylor Dispersion (parallel-
!!! plate and 'enhanced' rectangular geometry) effects.

!  Ao = (C1+C2+C3+C4+C5)/5
Ao = .5

do i = 1,Nx
	do j = 1,Ny
		sum_model = 0
! 		sum_model_nodiff = 0
! 		sum_model_old = 0
		sum_model_taylor = 0
! 		sum_model_taylor_enh = 0
		do m = 1,1000
!  			Am = (SqRt(2.0D0)/(m*pi))*((C1-C2)*sin(m*pi/5)+(C2-C3)*sin(2*m*pi/5) &
!  				+(C3-C4)*sin(3*m*pi/5)+(C4-C5)*sin(4*m*pi/5)+C5*sin(m*pi))
! 			Am = (SqRt(2.0D0)/((2*m-1)*pi))*((C1-C2)*sin((2*m-1)*pi/2))
			Am = SqRt(2.0D0)*(-1+cos(m*pi)+m*pi*sin(m*pi))/(m**2*pi**2)
			sum_model = sum_model + Am*exp(-2*m**2*pi**2*x(i) &
				/(Pe + SqRt(Pe**2 + 4*m**2*pi**2)))*cos(m*pi*y(j))
			sum_model_taylor = sum_model_taylor + Am*exp(-2*m**2*pi**2*x(i) &
				/(Pe + SqRt(Pe**2 + 4*k_tay*m**2*pi**2)))*cos(m*pi*y(j))
! 			sum_model_nodiff = sum_model_nodiff + Am*cos(m*pi*y(j))* &
! 				exp(-m**2*pi**2*x(i)/Pe)
! 			sum_model = sum_model + Am*exp(-2*(2*m-1)**2*pi**2*x(i) &
! 				/(Pe + SqRt(Pe**2 + 4*(2*m-1)**2*pi**2)))*cos((2*m-1)*pi*y(j))
! 			sum_model_nodiff = sum_model_nodiff + Am*cos((2*m-1)*pi*y(j))* &
! 				exp(-(2*m-1)**2*pi**2*x(i)/Pe)
! 			if (sum_model - sum_model_old <= 1e-10) then 
! 				exit
! 			end if
! 			sum_model_old = sum_model
		end do
		C_model(i,j) = Ao + SqRt(2.0D0)*sum_model
! 		C_model_nodiff(i,j) = Ao + Sqrt(2.0D0)*sum_model_nodiff
 		C_model_taylor(i,j) = Ao + SqRt(2.0D0)* sum_model_taylor
! 		C_model_taylor_enh(i,j) = Ao + SqRt(2.0D0) * sum_model_taylor_enh
    end do
end do

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_model')
WRITE (11, *) C_model
CLOSE (11)

! OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_model_nodiff')
! WRITE (11, *) C_model_nodiff
! CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_model_taylor')
WRITE (11, *) C_model_taylor
CLOSE (11)

! OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_model_taylor_enh')
! WRITE (11, *) C_model_taylor_enh
! CLOSE (11)

!!!! Calculate normalized velocities from classic analytical expression
!!!! for rectangular channel and store in matrix

! sum_u_ave = 0
! u_ave_old = 0
! do n = 1,1000
! 	sum_u_ave = sum_u_ave + tanh((2*n-1)*pi*w/(2*h))/(2*n-1)**5
! 	u_ave = h**2/12 * (1 - (192*h/(pi**5*w))*sum_u_ave)
! 	if (abs((u_ave - u_ave_old)/u_ave) < 1e-6) then
! 		exit
! 	end if
! 	u_ave_old = u_ave
! end do

! do j = 2,Ny - 1
! 	do k = 2,Nz - 1
! 		sum_u = 0
! 		u_old = 0
! 		do n = 1,1000
! 		! Necessary condition to account for because cosh rises to infinity extremely quickly, both terms rise at same rate so cosh(inf)/cosh(inf) ~= 1
! 			if (((2*n-1)*pi*w/(2*h)) < 700) then
! 				sum_term = (cosh((2*n-1)*pi*(y(j)-w/2)/h)/((2*n-1)**3 &
! 					* cosh((2*n-1)*pi*w/(2*h))))*sin((2*n-1)*pi*z(k)/h)				
! 			else
! 				sum_term = sin((2*n-1)*pi*z(k)/h)/(2*n-1)**3
! 			end if
! 			sum_u = sum_term + sum_u
! 			u(j,k) = (h*z(k)-(z(k))**2)/2 - (4*h**2/pi**3)*sum_u
! 			if (abs((u(j,k)-u_old)/u(j,k)) < 1e-6) then
! 				exit
! 			end if
! 			u_old = u(j,k)
! 		end do
! 	end do
! end do

do j = 2,Ny - 1
	do k = 2,Nz - 1
		u(j,k) = (6/h**2)*(h*z(k)-(z(k))**2)
	end do
end do

!u = u/u_ave

write(6,*) "velocity arrays calculated"

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'u')
WRITE (11, "(F8.6)") u
CLOSE (11)

!! Define local diffusion and convection at Control Vol. boundaries !!

do i = 2,Nx-1
	do j = 1,Ny
		do k = 1,Nz

		! Local diffusion condunctance
			if (i /= Nx-1) then
				D_e = 2*deltay(j)*deltaz(k)/(deltax(i)+deltax(i+1))
			else
				D_e = 0
			end if
			D_w = 2*deltay(j)*deltaz(k)/(deltax(i)+deltax(i-1))
			if (j /= Ny) then
				D_n = 2*deltaz(k)*deltax(i)/(deltay(j)+deltay(j+1))
			else
				D_e = 0
			end if
			if (j /= 1) then
				D_s = 2*deltaz(k)*deltax(i)/(deltay(j)+deltay(j-1)) 
			else
				D_s = 0
			end if
			if (k /= Nz) then
				D_t = 2*deltax(i)*deltay(j)/(deltaz(k)+deltaz(k+1))
			else
				D_t = 0
			end if
			if (k /= 1) then
				D_b = 2*deltax(i)*deltay(j)/(deltaz(k)+deltaz(k-1))
			else
				D_b = 0
			end if

			! Local convection conductance (0 when not in x direction)
			F_e = Pe*u(j,k)*deltay(j)*deltaz(k)
			F_w = Pe*u(j,k)*deltay(j)*deltaz(k)
			
			! Local Peclet
			if (i /= Nx-1 .and. j /= 1 .and. j /= Ny .and. k /= 1 & 
				.and. k /= Nz) then
				P_e = F_e/D_e
			else
				P_e = 0
			end if
			if (j /= 1 .and. j /= Ny .and. k /= 1 .and. k /= Nz) then
				P_w = F_w/D_w
			else
				P_w = 0
			end if
            
            ! Calculation of coefficients
			a_E(i,j,k) = D_e * max(0.0D0,(1-0.1*abs(P_e))**5) + max(-F_e,0.0D0)
			a_W(i,j,k) = D_w * max(0.0D0,(1-0.1*abs(P_w))**5) + max(F_w,0.0D0)
			a_N(i,j,k) = D_n
			a_S(i,j,k) = D_s 
			a_T(i,j,k) = D_t 
			a_B(i,j,k) = D_b 
			a_P(i,j,k) = a_E(i,j,k) + a_W(i,j,k) + a_N(i,j,k) & 
				+ a_S(i,j,k) + a_T(i,j,k) + a_B(i,j,k)

		end do
	end do
end do

write(6,*) "coefficients calculated"

!!!!!!!!!!!!!! Define Boundary Conditions !!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input BC concentrations in x

mean_tot = 0
do j = 1,Ny
	C(1,j,:) = y(j)
	mean_tot = mean_tot + C(1,j,1)*deltay(j)
end do

! mean_tot = 0
! do j = 1,Ny
! 	if (y(j) <= w/2) then
! 		C(1,j,:) = 0
! 	else if (y(j) > w/2) then
! 		C(1,j,:) = 1
! 	end if
! 	mean_tot = mean_tot + C(1,j,1)*deltay(j)
! end do

! mean_tot = 0
! do j = 1,Ny    
! 	if (y(j) <= w/5) then
! 		C(1,j,:) = 0
! 	else if (y(j) > w/5 .and. y(j) <= 2*w/5) then
! 		C(1,j,:) = .25
! 	else if (y(j) > 2*w/5 .and. y(j) <= 3*w/5) then
! 		C(1,j,:) = .5
! 	else if (y(j) > 3*w/5 .and. y(j) <= 4*w/5) then
! 		C(1,j,:) = .75
! 	else if (y(j) > 4*w/5) then
! 		C(1,j,:) = 1
! 	end if
! 	mean_tot = mean_tot + C(1,j,1)*deltay(j)
! end do

C(2:Nx,:,:) = mean_tot/w    !Start iterations with avg conc.

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'mean_tot')
WRITE (11, *) mean_tot
CLOSE (11)

write(6,*) "BCs defined"


!!! Calculate initial residual !!!!

errorsum = 0
do i = 2,Nx-1
	do j = 1,Ny
		do k = 1,Nz
			if (j == Ny .and. k == Nz) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
					- a_B(i,j,k)*C(i,j,k-1)
			else if (j == Ny .and. k == 1) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
					- a_T(i,j,k)*C(i,j,k+1)
			else if (j == 1 .and. k == Nz) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_B(i,j,k)*C(i,j,k-1)
			else if (j == 1 .and. k == 1) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_T(i,j,k)*C(i,j,k+1)
			else if (j == Ny) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
					- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
			else if (j == 1) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
			else if (k == Nz) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_S(i,j,k)*C(i,j-1,k) - a_B(i,j,k)*C(i,j,k-1)
			else if (k == 1) then
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1)
			else
				r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
					- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
					- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1) &
					- a_B(i,j,k)*C(i,j,k-1)
			end if
			errorsum = errorsum + r**2
		end do
	end do
end do
error = SqRt(errorsum)
print *, error
OPEN (UNIT=12, ACCESS="SEQUENTIAL", FILE=path//'r')
WRITE (12,*) '0', error
CLOSE (12)
error_orig = error

!!!!!!!!!! Initialize sweeping iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Sweep in planes from upstream to downstream until total solution is converged

do ni = 1,100000
	do i = 2,Nx-1
		!write(6,*) " "
		!disp(['plane ' num2str(i) ' of ' num2str(Nx) ', iteration ' num2str(ni)])

		!!!! Iterate line-by-line along a plane until the plane is converged
		do nj = 1,1000
			do j = 1,Ny
				! Calculate b, d, and Q coefficients
				Q(1)= 0
				a(1) = a_P(i,j,1)/alpha
				b(1) = a_T(i,j,1)
				c_it(1) = a_B(i,j,1)
				if ((j /= 1) .and. (j /= Ny)) then
					P(1) = a(1)/b(1)
				else
					P(1) = 0
				end if
				do k = 2,Nz
					a(k) = a_P(i,j,k)/alpha
					b(k) = a_T(i,j,k)
                	c_it(k) = a_B(i,j,k)
                	if (k == Nz) then
                		P(k) = 0
                	else
						P(k) = b(k)/(a(k)-c_it(k)*P(k-1))
                	end if
					if (j == 1) then
						b_P(k) = a_E(i,j,k)*C(i+1,j,k) + a_W(i,j,k)*C(i-1,j,k) &
							+ a_N(i,j,k)*C(i,j+1,k)
					else if (j == Ny) then
						b_P(k) = a_E(i,j,k)*C(i+1,j,k) + a_W(i,j,k)*C(i-1,j,k) &
							+ a_S(i,j,k)*C(i,j-1,k)
					else
						b_P(k) = a_E(i,j,k)*C(i+1,j,k) + a_W(i,j,k)*C(i-1,j,k) &
							+ a_N(i,j,k)*C(i,j+1,k) + a_S(i,j,k)*C(i,j-1,k)
					end if
					d(k) = b_P(k)+(1-alpha)*(a_P(i,j,k)/alpha)*C(i,j,k)
					if ((j == 1 .or. j == Ny) .and. k == Nz) then
						Q(k) = 0
					else
						Q(k) = (d(k)+c_it(k)*Q(k-1))/(a(k)-c_it(k)*P(k-1))
					end if
				end do

				! Solve for concentrations along a given line in z using TDMA
				C(i,j,Nz) = Q(Nz)
				do k = Nz-1,1,-1
					C(i,j,k) = P(k)*C(i,j,k+1)+Q(k)
				end do
			end do
            
			! Account for concentration at corners of a plane, which 
			! technically have no concentration flux in or out
			C(i,Ny,Nz) = C(i,Ny,Nz-1)
			C(i,Ny,1) = C(i,Ny,2)
			C(i,1,Nz) = C(i,1,Nz-1)
			C(i,1,1) = C(i,1,2)
            
            ! Calculate plane residual to determine convergence of plane
			errorsum_2D = 0
			do j = 1,Ny
				do k = 1,Nz
					if (j == Ny .and. k == Nz) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
							- a_B(i,j,k)*C(i,j,k-1)
					else if (j == Ny .and. k == 1) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
							- a_T(i,j,k)*C(i,j,k+1)
					else if (j == 1 .and. k == Nz) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_B(i,j,k)*C(i,j,k-1)
					else if (j == 1 .and. k == 1) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_T(i,j,k)*C(i,j,k+1)
					else if (j == Ny) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
							- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
					else if (j == 1) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
					else if (k == Nz) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_S(i,j,k)*C(i,j-1,k) - a_B(i,j,k)*C(i,j,k-1)
					else if (k == 1) then
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1)
					else
						r_2D = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
							- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
							- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1) &
							- a_B(i,j,k)*C(i,j,k-1)
					end if
					errorsum_2D = errorsum_2D + (r_2D)**2
				end do
			end do
			error_2D = SqRt(errorsum_2D)
			! fprintf('!.5g \n',error_2D)
			if (error_2D <= residual_drop_2D * error) then 
				exit
			end if
		end do
	end do
    
!!!!! Calculate full 3D residual to determine convergence !!!!!!!!!!!!!

	errorsum = 0
	do i = 2,Nx-1
		do j = 1,Ny
			do k = 1,Nz
				if (j == Ny .and. k == Nz) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
						- a_B(i,j,k)*C(i,j,k-1)
				else if (j == Ny .and. k == 1) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) & 
						- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) &
						- a_T(i,j,k)*C(i,j,k+1)
				else if (j == 1 .and. k == Nz) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_B(i,j,k)*C(i,j,k-1)
				else if (j == 1 .and. k == 1) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_T(i,j,k)*C(i,j,k+1)
				else if (j == Ny) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_S(i,j,k)*C(i,j-1,k) & 
						- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
				else if (j == 1) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_T(i,j,k)*C(i,j,k+1) - a_B(i,j,k)*C(i,j,k-1)
				else if (k == Nz) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_S(i,j,k)*C(i,j-1,k) - a_B(i,j,k)*C(i,j,k-1)
				else if (k == 1) then
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1)
				else
					r = a_P(i,j,k)*C(i,j,k) - a_E(i,j,k)*C(i+1,j,k) &
						- a_W(i,j,k)*C(i-1,j,k) - a_N(i,j,k)*C(i,j+1,k) &
						- a_S(i,j,k)*C(i,j-1,k) - a_T(i,j,k)*C(i,j,k+1) &
						- a_B(i,j,k)*C(i,j,k-1)
				end if
					errorsum = errorsum + (r)**2
			end do
		end do
	end do
	error = SqRt(errorsum)
	print *, error
!     OPEN (UNIT=12, ACCESS="SEQUENTIAL", FILE=path//'r', POSITION="APPEND")
!     WRITE (12, *) ni, error
!     CLOSE (12)
	if (error < residual_drop * error_orig) then
		exit
	end if
end do

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C')
WRITE (11, *) C
CLOSE (11)

!!! Compute transverse and laterally averaged concentration
!!! profiles

do i = 1,Nx-1
	do j = 1,Ny
		sum_ave_z = 0;
		do k = 1,Nz
			sum_ave_z = sum_ave_z + deltaz(k)*C(i,j,k);
		end do
			C_ave_z(i,j) = sum_ave_z/h;
	end do
end do

do i = 1,Nx-1
	do k = 1,Nz
		sum_ave_y = 0;
		do j = 1,Ny
			sum_ave_y = sum_ave_y + deltay(j)*C(i,j,k);
		end do
		C_ave_y(i,k) = sum_ave_y/w;
	end do
end do

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_ave_y')
WRITE (11, *) C_ave_y
CLOSE (11)

OPEN (UNIT=11, ACCESS="SEQUENTIAL", FILE=path//'C_ave_z')
WRITE (11, *) C_ave_z
CLOSE (11)

end program model
