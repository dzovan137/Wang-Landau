program shuff
implicit none
integer :: n_bins,i
double precision, dimension(:), allocatable :: abscisse, tempt
character(len=40) :: filename

n_bins = 9

allocate(abscisse(n_bins),tempt(n_bins))

filename='test.dat'
write(*,*) filename
open(11,file=filename,status='unknown')
filename='4.dat'
open (15,file=filename,status='unknown')
do i=1,n_bins
	read(11,*) abscisse(i),tempt(i)
	
end do



do i=3,n_bins
write(15,*) abscisse(i),tempt(1) - tempt(2) +tempt(i)
write(*,*) i,abscisse(i), tempt(i)
end do
 close(15)



end program shuff
