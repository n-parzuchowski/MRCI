program mrci_main
  use mrci_solve
  implicit none

  integer :: ii,jj,AA,num_refs,Aprot,Aneut ,q,t,lj
  character(200) :: input,ref_file,sp_file,int_file,den_file
  integer,allocatable,dimension(:,:) :: ref,basis
  real(8) :: x,sm
  real(8) :: omp_get_Wtime
  integer :: i,j,k,l,Jtot
  
  time_Zero = omp_get_wtime()
  
  call print_header
  call getarg(1,input)
  if (trim(input) == '') STOP "NEED INPUT FILE" 

  call read_input_file(input,sp_file,int_file,ref_file,AA,Aprot,Aneut) 
  mbas%Abody = AA
  mbas%Aprot = Aprot
  mbas%Aneut = Aneut
  
  call init_sp_basis(sp_file)

  write(*,"(A,I8)") "Number of sp states: ", mbas%ntot
  
  call read_Ref_file(ref_file,num_refs,REF,AA,den_file)

  call generate_basis(ref,basis,0,6)

  call generate_tp_basis

  ! get interaction
  call read_me2b(me2b,int_file)
  call read_me1b(me0b,me1b,int_file)

 
  ! get density matrix
  call read_me2b(lambda2b,den_file)
  call read_me1b(x,lambda1b,den_file)

  call traces
  write(*,"(A,f12.4)") "Normal-ordered E0: ", ME0B
  call unnormal_order(me0b,me1b,me2b)

  call diagonalize(me0b,me1b,me2b,basis)
end program mrci_main
  
subroutine print_header
  print*; print*; print*;
  write(*,'(A)') "============================================="
  write(*,'(A)') "MULTIREFERENCE CONFIGURATION INTERACTION CODE"
  write(*,'(A)') "============================================="
  print*; print*; print*;
end subroutine
