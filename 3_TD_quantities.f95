program Calc
implicit none
double precision, dimension(:), allocatable :: abscisse,g_e,p_e,z,tempt,mag,abscissem
integer :: i,L,n_bins,j
real :: temp
double precision :: lambda,part,c_v,u1,u2,f,s,p,summ,mag1,mag2,susc,cumulant,mag4
character(len=40) :: filename
character(len=8) :: fmt,x1

! initial paramters
L = 16
n_bins = 49

allocate (p_e(n_bins),abscisse(n_bins),g_e(n_bins),z(n_bins),tempt(n_bins),mag(n_bins),abscissem(n_bins))
abscisse(:) = 0
p_e(:) = 0
g_e(:) = 0
z(:) = 0
tempt(:) = 0
mag(:) = 0
abscissem(:) = 0

! this loading in the magnetization 
open(12,file='data/00008DOSq8magnetization00001.dat')
do i=1,n_bins
 read(12,*) mag(i)
 end do
 close(12)


! loading in the DOS
filename='data/000016density_of_stateq800001.dat'
write(*,*) filename
open(11,file=filename,status='unknown')
do i=1,n_bins
	read(11,*) abscisse(i),tempt(i)
end do
close(11)


filename='data/everything1.dat'
open (15,file=filename,status='unknown')
filename='data/everything2.dat'
open (16,file=filename,status='unknown')

! initial temperature
temp = 0.4
! loop over temperature
do j=1,100


do i=1,n_bins
g_e(i) = tempt(i) - (abscisse(i)*L*L)/real(temp)
  end do

lambda = MAXVAL(g_e)
!write(*,*) lambda
part = 0
do i=1,n_bins
part = part + dexp(g_e(i) - lambda)
  end do
  
!calculation of internal energy
u1=0
do i=1,n_bins
u1 = u1 + abscisse(i)*dexp(g_e(i) - lambda)
  end do

u1 = u1/part

u2 = 0
do i=1,n_bins
u2 = u2 + (abscisse(i)*abscisse(i))*dexp(g_e(i) - lambda)
  end do


u2 = u2/part
!Calculation of specific heat
 c_v = (u2 - u1*u1)/(temp**2)

!Calculation of free energy
part = 0
do i=1,n_bins
part = part + dexp(g_e(i))
  end do

f = -temp*LOG(part)


! Calculation of entropy
s = (u1-f/(L*L))/real(temp)

! Calculation of magnetization and susceptibility
part = 0
do i=1,n_bins
part = part + dexp(g_e(i) - lambda)
  end do
mag1 = 0
	do i=1,n_bins
		write(*,*) mag(i)
		mag1 = mag1 + abs(mag(i))*dexp(g_e(i) - lambda)
	end do

mag1 = mag1/part

mag2 = 0
	do i=1,n_bins
		mag2 = mag2 + abs(mag(i))*abs(mag(i))*dexp(g_e(i) - lambda)
	end do

susc = (mag2 - mag1*mag1)/temp 


! binder cumulant 

mag4 = 0
do i=1,n_bins
	mag4 = mag4 + abs(mag(i))*abs(mag(i))*abs(mag(i))*abs(mag(i))*dexp(g_e(i) - lambda)
end do
write (*,*) 'Cumulant:',mag4,mag2
 cumulant = 1 - mag4/(3*mag2*mag2)

    

write(15,*) temp,u1,c_v*L*L,mag1,cumulant
write(16,*) temp,f/(L*L),s,susc

temp = temp + 0.01
end do
  close(15)
  close(16)

!calculation of probability
do j=1,10
summ = 0
temp = 0.756 + (j-1)*0.00005

do i=1,n_bins
g_e(i) = tempt(i) - (abscisse(i)*L*L)/real(temp)
  end do



fmt = '(I5.5)'
write (x1,fmt) j
filename='data/probability'//trim(x1)//'.dat'
open (17,file=filename,status='unknown')
part = 0
do i=1,n_bins
part = part + dexp(g_e(i))
  end do
!write(*,*) part
do i=1,n_bins
p = dexp(g_e(i))/part
summ = summ + p
write(17,*) abscisse(i),p
end do

 close(17)

write(*,*) summ
end do




end program Calc

