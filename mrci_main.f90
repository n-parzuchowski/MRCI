program mrci_main
  use mrci_solve
  implicit none

  integer :: ii,jj,AA,num_refs,Aprot,Aneut ,q,t,lj
  character(200) :: input,ref_file,sp_file,int_file,den_file,law_file
  character(200),allocatable,dimension(:) :: obs_files
  integer,allocatable,dimension(:,:) :: ref,basis
  real(8) :: x,sm, Eimsrg
  real(8) :: omp_get_Wtime
  integer :: i,j,k,l,Jtot
  
  time_Zero = omp_get_wtime()
  
  call print_header
  call getarg(1,input)
  if (trim(input) == '') STOP "NEED INPUT FILE" 

  call read_input_file(input,sp_file,int_file,den_file,ref_file,law_file,obs_files) 
  call init_sp_basis(sp_file)

  write(*,"(A,I8)") "Number of sp states: ", mbas%ntot
  
  call read_Ref_file(ref_file,num_refs,REF)
  call generate_basis(ref,basis)
  call generate_tp_basis

  
  ! get interaction
  call read_me2b(me2b,int_file)
  call read_me1b(me0b,me1b,int_file)

  ! get lawson part 
  if (non_negligible(mbas%law_beta)) then
     call read_me2b(Hcm2b,law_file)
     call read_me1b(Hcm0b,Hcm1b,law_file)
     call lawsonize(mbas%law_beta)
  end if
  
  ! get density matrix
  if (num_refs > 1) then 
     call read_me2b(lambda2b,den_file)
     call read_me1b(x,lambda1b,den_file)
  else
     call generate_2b_density(lambda2b)
     call generate_1b_density(lambda1b,REF(1,:))
  end if
  
  call traces
  write(*,"(A,f12.4)") "Normal-ordered E0: ", ME0B
  Eimsrg = ME0B
  call unnormal_order(me0b,me1b,me2b)

  call diagonalize(me0b,me1b,me2b,basis,Eimsrg,10) 

  !!! handle any observables
  do q = 1, mbas%num_obs
     call read_me2b(Op2b,obs_files(q))
     call read_me1b(Op0b,Op1b,obs_files(q))

     write(*,*)
     write(*,"(A)") "Computing observables from "//trim(obs_files(q))

     write(*,"(A,f12.4)") "Normal-ordered Op0: ", Op0b
     call unnormal_order(Op0b,Op1b,Op2b)

     call observable(Op0b,Op1b,Op2b,basis,10)            
  end do
  
end program mrci_main
  
subroutine print_header
  print*; print*; print*;
  write(*,'(A)') "============================================="
  write(*,'(A)') "MULTIREFERENCE CONFIGURATION INTERACTION CODE"
  write(*,'(A)') "============================================="
  print*; print*; print*;
end subroutine
