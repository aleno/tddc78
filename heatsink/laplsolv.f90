program laplsolv
use omp_lib
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem 
! on a square using the Jacobi method. 
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  integer, parameter                  :: n=1000, maxiter=1000
  double precision,parameter          :: tol=1.0E-3
  double precision,dimension(0:n+1,0:n+1) :: T
  double precision,dimension(n)       :: tmp1,tmp2,tmp3,Terr
  double precision                    :: error,x, err
  double precision                    :: t1,t0
  integer                             :: i,j,k, Tl, id, last
  character(len=20)                   :: str

  integer                             :: nthreads, istart, iend
  
  ! Set boundary conditions and initial values for the unknowns
  T=0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0

  ! Solve the linear system of equations using the Jacobi method

  t0 = omp_get_wtime();

  !$omp parallel private(id, istart, iend, tmp1, tmp2, tmp3, k, j)

  !$omp single
  nthreads = omp_get_num_threads()
  lines_per_thread = n / nthreads + 1
  !$omp end single

  id = omp_get_thread_num();
  istart = id * lines_per_thread + 1
  iend = min((id +1) * lines_per_thread, n)

  !$omp critical
  write (*,*) 'Hello from', id, istart, iend
  !$omp end critical

  do k=1,maxiter
    
     tmp1=T(1:n,0)
     error=0.0D0
     
!     do j=1,n
!        tmp2=T(1:n,j) ! Hämta kollumn j med rad 1:n.
!        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1)/4.0D0
!        error=max(error,maxval(abs(tmp2-T(1:n,j))))
!        tmp1=tmp2
!     end do

     tmp1 = T(1:n,istart - 1)
     tmp3 = T(1:n,iend + 1)

     !$omp barrier
     !$omp flush

     do j = istart, iend - 1
        tmp2 = T(1:n,j)
        T(1:n,j) = (T(0:n-1,j) + T(2:n+1,j) + T(1:n,j+1) + tmp1)/4.0D0

        !$omp critical
        error = max(error, maxval(abs(tmp2 - T(1:n,j))))
        !$omp end critical

        tmp1 = tmp2
     end do

     ! Specialfall för gränsgrejer
     tmp2 = T(1:n, iend)
     T(1:n,iend) = (T(0:n-1,iend) + T(2:n+1,iend) + tmp3 + tmp1)/4.0D0
     !$omp critical
     error = max(error, maxval(abs(tmp2 - T(1:n,iend))))
     !$omp end critical

     !$omp barrier
     !$omp flush

     if (error<tol) then
        exit
     end if
     !$omp barrier
  end do
  !$omp end parallel  

  t1 = omp_get_wtime();

  write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',k
  write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

  ! Uncomment the next part if you want to write the whole solution
  ! to a file. Useful for plotting. 
  
  !open(unit=7,action='write',file='result.dat',status='unknown')
  !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
  !do i=0,n+1
  !   write (unit=7,fmt=str) T(i,0:n+1)  
  !end do
  !close(unit=7)
  
end program laplsolv
