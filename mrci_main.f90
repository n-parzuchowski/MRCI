program mrci_main
  use mrci_interactions
  implicit none

  type(spd) :: mbas,jbas
  type(tpd) :: tp_basis
  integer :: ii,AA,num_refs,Aprot,Aneut ,q
  character(200) :: input,ref_file,sp_file,int_file
  integer,allocatable,dimension(:,:) :: ref,basis

  call getarg(1,input)
  if (trim(input) == '') STOP "NEED INPUT FILE" 

  call read_input_file(input,sp_file,int_file,ref_file,AA,Aprot,Aneut) 

  call init_sp_basis(mbas,jbas,sp_file)
  
  mbas%Abody = AA
  mbas%Aprot = Aprot
  mbas%Aneut = Aneut

  call read_Ref_file(ref_file,num_refs,REF,AA)

  call generate_basis(mbas,ref,basis)
  call generate_tp_basis(tp_basis,jbas)
  call read_me2b(int_file,tp_basis)

  print*, ME2B(1)%X(1:12)
  
end program mrci_main
  
