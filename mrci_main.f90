program mrci_main
  use mrci_basis
  implicit none

  type(spd) :: jbas
  integer :: ii,AA,num_refs,Aprot,Aneut 
  character(200) :: fname
  integer,allocatable,dimension(:,:) :: ref,basis

  
  print*, 'Enter SP file: '
  read*, fname 
  call jbas%build(fname)

  print*, 'Enter number of nucleons'
  read*, AA,Aprot,Aneut

  jbas%Abody = AA
  jbas%Aprot = Aprot
  jbas%Aneut = Aneut  
  print*, 'Enter number of references'
  read*, num_refs

  allocate(ref(num_refs,AA))

  write(*,'(A,I3,A)') 'Enter the ',num_refs, 'reference states, one by one'
  do ii = 1, num_refs
     read*, ref(ii,:)
  end do
  
  call generate_basis(jbas,ref,basis)

  ! do ii = 1,jbas%Ntot
  !    print*, ii, jbas%nn(ii),jbas%ll(ii),jbas%jj(ii),jbas%mm(ii),jbas%tz(ii)
  ! end do

!
  
end program mrci_main
  
