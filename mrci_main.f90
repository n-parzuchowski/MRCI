program mrci_main
  use mrci_interactions
  implicit none

  type(tpd) :: tp_basis
  integer :: ii,AA,num_refs,Aprot,Aneut ,q,t,lj
  character(200) :: input,ref_file,sp_file,int_file,den_file
  integer,allocatable,dimension(:,:) :: ref,basis
  real(8) :: x
  
  call getarg(1,input)
  if (trim(input) == '') STOP "NEED INPUT FILE" 

  call read_input_file(input,sp_file,int_file,ref_file,AA,Aprot,Aneut) 

  call init_sp_basis(sp_file)
  
  mbas%Abody = AA
  mbas%Aprot = Aprot
  mbas%Aneut = Aneut

  call read_Ref_file(ref_file,num_refs,REF,AA,den_file)

  call generate_basis(ref,basis)
  call generate_tp_basis(tp_basis)
  call read_me2b(me2b,int_file,tp_basis)
  call read_me1b(me0b,me1b,int_file)

  ! get density matrix
  call read_me2b(lambda2b,den_file,tp_basis)
  call read_me1b(x,lambda1b,den_file)


  do ii = 1, mbas%ntot
     print*,  ii, mbas%nn(ii),mbas%ll(ii),mbas%jj(ii),mbas%mm(ii),mbas%tz(ii)
  end do
  
  print*, ME2B(1)%X(1:12)

  print*, get_ME1B(24,88,ME1B)
 
  x=get_ME2B(1,2,1,2,ME2B,tp_Basis) 
  print*, x
  x=get_JME2B(1,1,1,1,0,ME2B,tp_Basis) 
  print*, x

  print*, ME0B

  call unnormal_order(me0b,me1b,me2b,tp_basis)

  print*, ME2B(1)%X(1:12)

  print*, get_ME1B(24,88,ME1B)
 
  x=get_ME2B(1,2,1,2,ME2B,tp_Basis) 
  print*, x
  x=get_JME2B(1,1,1,1,0,ME2B,tp_Basis) 
  print*, x

  print*, ME0B
  
end program mrci_main
  
