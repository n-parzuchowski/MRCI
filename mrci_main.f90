program mrci_main
  use mrci_solve
  implicit none

  integer :: ii,jj,AA,num_refs,Aprot,Aneut ,q,t,lj
  character(200) :: input,ref_file,sp_file,int_file,den_file
  integer,allocatable,dimension(:,:) :: ref,basis
  real(8) :: x!,!AM(1040,1040)
  
  call getarg(1,input)
  if (trim(input) == '') STOP "NEED INPUT FILE" 

  call read_input_file(input,sp_file,int_file,ref_file,AA,Aprot,Aneut) 
  mbas%Abody = AA
  mbas%Aprot = Aprot
  mbas%Aneut = Aneut
  
  call init_sp_basis(sp_file)
  
  ! do ii = 1, 50
  !    write(*,'(6(I5))') ii, mbas%nn(ii),mbas%ll(ii),mbas%jj(ii),mbas%mm(ii),mbas%tz(ii)
  !    end do
  
  
  call read_Ref_file(ref_file,num_refs,REF,AA,den_file)

  call generate_basis(ref,basis)
  call generate_tp_basis

  ! get interaction
  call read_me2b(me2b,int_file)
  call read_me1b(me0b,me1b,int_file)
  
  ! get density matrix
  call read_me2b(lambda2b,den_file)
  call read_me1b(x,lambda1b,den_file)

  print*, "E0 ground state: ", ME0B
  call unnormal_order(me0b,me1b,me2b)
  print*, "Unnormal ordered offset: ",ME0B

  call diagonalize(me0b,me1b,me2b,basis)
end program mrci_main
  
