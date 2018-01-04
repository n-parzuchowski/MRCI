program print_system
  use mrci_aux
  implicit none

  character(200) :: input,ref_file,sp_file,int_file,den_file,law_file
  character(200),allocatable,dimension(:) :: obs_files
  real(8),allocatable,dimension(:) :: e,s
  real(8),allocatable,dimension(:,:) :: v
  integer :: Nstates,dim, II 
  character(2) :: Nstr
  real(8) :: spin 

  call getarg(1,input) 
  
  if (trim(input) == '') STOP "NEED INPUT FILE"
  
  call read_input_file(input,sp_file,int_file,den_file,ref_file,law_file,obs_files)

  
  call read_results(e,s,v)


  Nstates= size(e)
!!  write(Nstr,"(I2)") Nstates

  print*
  write(*, "(A)") "================================================="
  write(*, "(A)") "     SPIN      ENERGY       EX ENERGY   " 
  write(*, "(A)") "================================================="

  do II = 1, Nstates
     spin = (sqrt(s(II)*4+1.d0 )-1)/2.d0  
     write(*,"(f9.1,2(f12.4))") spin,e(ii),e(ii)-e(1)
  end do


end program print_system

     
  
  
  
  
  
