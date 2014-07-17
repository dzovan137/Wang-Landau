program WangLand
implicit none
integer,dimension (:), allocatable :: qlist
integer,dimension(:,:),allocatable :: spin
integer :: i,j,L,q,n_bins,mc_steps,nb_f,w,u,o
double precision :: ran1,ran2,ran3,ran4,bin_width
double precision, dimension(:), allocatable ::  histo,g_e,abscisse,abscissem,histom,g_em! this is logaritmic input here inside g_e
!double precision, dimension(:,:), allocatable :: histo,g_e ! 2D array for 
integer, dimension(:), allocatable :: tampon
double precision :: energy,f,e1,e2,magnetization,m1,m2,ener_min,ener_max
integer :: b,s,t,k,r,mcs,g,h,cc1,cc2,tampoon
double precision :: invert,prob
integer :: histogram_check
character(len=40) :: filename
character(len=8) :: fmt,x1,x2
real :: qq





! RANDOM STUFF FROM THE CLOCK 
INTEGER ::  gj, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = gj)
ALLOCATE(seed(gj))        
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, gj) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)


! initial paramters*****************************
L = 8 !size of the system
q = 3 !how many potts states
n_bins = 32 !number of bins for the energy range -2 and 0
bin_width = 4./(L**2) ! size of a bin
mc_steps = 20000000 !number of Monte Carlo steps
f = 1. ! initial value for the weighting factor f
f = dexp(f)
nb_f = 27 ! how many times we decrease f
w = 4 ! WHICH RANGE WE WANT TO EXPLORE W=1,2,3,4



!energy range
ener_min = -2 + 0.5*(w-1)
ener_max = -1.4 + 0.5*(w-1)

if (w .eq. 4) then
ener_max = -0.05
end if



! initial allocation of size of the arrays*************
allocate(histo(n_bins),abscisse(n_bins),g_e(n_bins),abscissem(n_bins),histom(n_bins),g_em(n_bins))
allocate(spin(1:L,1:L)) ! spin 2D array
allocate(qlist(1:q)) !for different potts states
allocate(tampon(1:n_bins))
spin(:,:) = 0
qlist(:) = 0
abscisse(:) = 0
histo(:) = 0
g_e(:) = 1
abscissem(:) = 0
histom(:) = 0
g_em(:) = 0
tampon(:) = 0


! setting up the q list*************
do i=1,q
  qlist(i) = i
end do


!initial lattice***************
! the setup is a bit different for every range just to speed up
 cc2 = 0
do i=1,L
 do j=1,L 

	if (w .eq. 1 .or. w .eq. 2) then
	spin(i,j) = 1
	end if

	if (w .eq. 3) then 
	call random_number(ran1)
	b = int(ran1*q)+1
	spin(i,j) = qlist(b)
	end if

	if (w .eq. 4) then
	cc2 = cc2 + 1
	tampoon = MOD(cc2,3) +1
	spin(i,j) = tampoon
	end if

 end do
end do

!write(*,*) 'Energy:',energy(spin,L)
!write(*,*) 'Magnetization:', magnetization(spin,L,q) 

!abscisse division energy and magnetization*****************!
do i=1,n_bins
abscissem(i) = -1 + i*bin_width/2.
end do
do i=1,n_bins
	abscisse(i) = -2 + i*bin_width
	!write(*,*) abscisse(i)
end do

write(*,*) 'LATTICE',L,'x',L
write (*,*) ener_min,'<range<',ener_max


! get us in range
call ener_in_range(ener_min,ener_max,spin,L,q,qlist)



!weighting factor loop*********************
do i=1,nb_f
write(*,*) i
! Monte Carlo loop***********************
do mcs=1,mc_steps

! one sweep
do j=1,L*L

! computing the energy of a configuration**********
e1 = energy(spin,L)
m1 = magnetization(spin,L,q)
! selecting a random spin to flip to a random value from qlist*******

! look for a suitable change of spins
o = 0
do while(o < 1)
call random_number(ran1)
call random_number(ran2)
call random_number(ran3)
s = int(ran1*L)+1
t = int(ran2*L)+1
b = int(ran3*q)+1

e2 = invert(spin,L,q,qlist,s,t,b)
  if (e2 .ge. ener_min .and. e2 .le. ener_max) then
  o = 10

  end if

end do



m2 = magnetization(spin,L,q)
!finding the position in the energy range for e1, e2 (energy before flip, energy after flip respectively)****************
k = int((2 + e1)/(real(bin_width)))
k = min(k, n_bins)
k = max(k,1)

r = int((2 + e2)/(real(bin_width)))
r = min(r, n_bins)
r = max(r,1)
! finding the position for the magnetization
g = int((1 - m1)/(real(bin_width/2.)))
g = min(g, n_bins)
g = max(g,1)

h = int((1 - m2)/(real(bin_width/2.)))
h = min(h, n_bins)
h = max(h,1)





!write(*,*) k,r,g_e(k),g_e(r)
!write(*,*) e1,e2,m1,m2,g,h


! probability of acceptance*************************
if (g_e(r) .le. g_e(k)) then
	prob = 1.
else
	prob = dexp(g_e(k) - g_e(r))
end if

call random_number(ran4)
if (ran4 .le. prob) then
	
	! energy
	histo(r) = histo(r) + 1
        g_e(r) = g_e(r) + LOG(f)
        !now we accept the new configuration
	!write(*,*) 'ACCEPTED!'
	spin(s,t) = qlist(b)

 	! magnetization
	histom(h) = histo(h) + 1
	g_em(h) = g_em(h) + m2
 
else
	
	! energy
	histo(k) = histo(k) + 1
	g_e(k) = g_e(k) + LOG(f)
	
	! magnetization
	histom(g) = histo(g) + 1
	g_em(g) = g_em(g) + m1


end if



end do
u = 10
if (MOD(mcs,10000) == 0) then
write(*,*) mcs
u = histogram_check(histo,n_bins,mcs)
end if

if (u < 5) then
write(*,*) "Histogram flat at:",mcs
EXIT
end if






end do
!********************************


f = SQRT(f)
if (i .eq. nb_f) then

else
histo(:)=0
end if
end do
!*******************************







! writing in the files***************************

fmt = '(I5.5)'
write (x1,fmt) L
write (x2,fmt) w
filename='data/'//trim(x1)//'histogramq10'//trim(x2)//'.dat'
open (11,file=filename,status='unknown')

filename='data/'//trim(x1)//'density_of_stateq10'//trim(x2)//'.dat'
open (14,file=filename,status='unknown')

filename='data/'//trim(x1)//'DOSmagq10'//trim(x2)//'.dat'
open (18,file=filename,status='unknown')


 cc1 = 0
qq = q
do i=1,n_bins
if (histo(i) > 0) then
 cc1 = cc1 + 1 
 tampon(cc1) = i
! energy
write (14,*) abscisse(i),g_e(i)-g_e(tampon(1)) + LOG(qq)
write (11,*) abscisse(i),histo(i) ! the last histogram recorded, just to show the 'flatness'
! magnetization
write (18,*) abscissem(i),abs(g_em(i))/histom(i)/(L*L)
end if
end do




write(*,*) 'SIMULATION ENDED!'

end program WangLand











!***********************************************************!
!function to calculate the energy of a configuration		
function energy(spin,L)
implicit none
integer :: L,i,j
integer,dimension(L,L) :: spin
real:: ener
integer,dimension (:), allocatable :: ip,im
double precision:: energy

allocate (ip(1:L),im(1:L)) ! for periodic bc

! Periodic Boundary Conditions Array
do i=1,L
ip(i) = i+1
im(i) = i-1
end do
ip(L)=1
im(1)=L

energy = 0.
do i=1,L
  do j=1,L
	ener = 0
		if (spin(ip(i),j) .eq. spin(i,j)) then
			ener = ener - 1
    		end if
		if (spin(im(i),j) .eq. spin(i,j)) then
			ener = ener - 1
    		end if
    		if (spin(i,ip(j)) .eq. spin(i,j)) then
			ener = ener - 1
    		end if
    		if (spin(i,im(j)) .eq. spin(i,j)) then
			ener = ener - 1
    		end if

	energy = energy + ener/2

  end do
end do
energy = energy/(L*L)

return
end function energy
!*********************************************************!



!*********************************************************!
! function to invert one spin in the given configuration
function invert(spin,L,q,qlist,s,t,b)
integer :: L,q,s,t,b
double precision :: invert,energy
integer,dimension(L,L) :: spin,aux_spin
integer,dimension(q) :: qlist


aux_spin = spin
aux_spin(s,t) = qlist(b)
invert = energy(aux_spin,L)


return
end function invert
!********************************************************!

!********************************************************!
! function that checks if histogram is flat
function histogram_check(histo,n_bins,mcs)
implicit none
integer :: n_bins,mcs
double precision, dimension(n_bins) :: histo
integer :: histogram_check,u,count1,cc1
real :: summ,optimal,proc_optimal


! checking if the historgram is flat

summ = 0
 cc1 = 0
do u=1,n_bins
if (histo(u) > 0) then
  cc1 = cc1 + 1
 summ = summ + histo(u) 
end if
end do
write(*,*) mcs


! check the deviation for 80% flatness
optimal = summ/(real(cc1))
proc_optimal = optimal/5

count1 = 0
do u=1,n_bins
if (histo(u) > 0) then
if((optimal-proc_optimal > histo(u)) .OR. (histo(u) - optimal < proc_optimal)) then
count1 = count1 + 1
  end if
end if
end do

if (cc1 .eq. count1) then
	histogram_check = 1
else
	histogram_check = 10
end if



return
end function histogram_check
!***********************************************************!



!**********************************************************!
!calculates the magnetization for the Potts model with a specific formula
function magnetization(spin,L,q)
implicit none
integer :: L,i,j,q,p
real :: summ
integer, dimension(L,L) :: spin
double precision :: magnetization

summ = 0
do i=1,L
 do j=1,L
	
	if(spin(i,j) .eq. 1) then
	 	p = 1
	else 
		p = 0
	end if
	summ = summ + q*p-1
 end do
end do

magnetization = summ/(L*L*(q-1))
!write(*,*) magnetization
return
end function magnetization
!********************************************************!



!********************************************************!
! getting the initial configuration in desired range
subroutine ener_in_range(ener_min,ener_max,spin,L,q,qlist)
implicit none
integer :: L,i,j,t,q,b
integer,dimension(L,L) :: spin
integer, dimension(q) :: qlist
real :: ran1,ran2,ran3
real :: t1,t2
double precision :: energy,ener_max,ener_min,e1
! RANDOM STUFF FROM THE CLOCK 

INTEGER ::  gj, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = gj)
ALLOCATE(seed(gj))        
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, gj) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)



e1 = energy(spin,L)

t = 0
call cpu_time(t1)
do while(t < 1)

  if ((e1 .ge. ener_min .and. e1 .le. ener_max) .and. (((ener_min+ener_max)/2-e1) .le. 0.05 )) then

t = 5
write(*,*) 'Found it, we are in range!',e1


  else

 call random_number(ran1)
 call random_number(ran2)
 call random_number(ran3)
i = int(L*ran1) + 1
j = int(L*ran2) + 1
b = int(ran3*q)+1
spin(i,j) = qlist(b)

e1 = energy(spin,L)
write(*,*) e1,b,i,j
end if
  end do
call cpu_time(t2)
write(*,*) 'How much time to be in range:',t2-t1

end subroutine ener_in_range






