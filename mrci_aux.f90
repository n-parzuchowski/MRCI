module mrci_aux
!!! here we keep the auxillary functions/subroutines for this code
!!! This is the most fundamental building block for the code. 
!!!
!!! character(6) function spec_not(t,n,j,l) spectral notation printer  
!!!
implicit none
!!!
TYPE :: int_vec
   integer,allocatable,dimension(:) :: Z
END TYPE int_vec


!!!===========================================================
!!!===========================================================
TYPE :: spd    ! single particle discriptor
   integer :: Ntot, Jmax, Lmax,Njjcoup,nMax
   integer :: n,l,j,t,m,Abody,Aneut,Aprot
   integer :: Ptarg,Mtarg,dTz,num_obs
   integer :: Atarg, Ntarg , Ztarg
   real(8) :: e,law_beta
   integer, allocatable, dimension(:) :: nn, ll, jj, mm,tz,jlab
   type(int_vec),allocatable,dimension(:) :: amap,qmap
   type(int_vec),allocatable,dimension(:,:) :: tlj_to_ab 
   ! for clarity:  nn, ll, nshell are all the true value
   ! jj is j+1/2 (so it's an integer) 
   ! likewise itzp is 2*tz
END TYPE spd
!!!===========================================================
!!!===========================================================
TYPE :: block_descript
   integer :: aMax,J,Tz,PAR
   integer,allocatable,dimension(:,:) :: qnums
END type block_descript
!!!===========================================================
!!!===========================================================
TYPE :: tpd
   integer :: bMax,aMaxMax
   type(block_descript),allocatable,dimension(:) :: block 
END type tpd
!!!===========================================================
!!!===========================================================
  
type(spd),public :: jbas,mbas
type(tpd),public :: tp_basis 
character(500) :: ME_DIR,SP_DIR,INI_DIR,OUTPUT_DIR
character(200) :: prefix
integer :: hw
real(8) :: tot_memory,time_Zero
contains

  subroutine print_memory(mem)
    implicit none

    real(8),intent(in) :: mem   
    
    if (mem > 1024.**3 ) then
       write(*,'(A,f6.2,A)') "Memory: ",mem/(1024.**3)," GB"
    else if (mem > 1024.**2 ) then   
       write(*,'(A,f6.2,A)') "Memory: ",mem/(1024.**2)," MB"
    else if (mem > 1024. ) then   
       write(*,'(A,f6.2,A)') "Memory: ",mem/(1024.)," KB"
    else
       write(*,'(A,f6.2,A)') "Memory: ",mem," B"
    end if
           
  end subroutine print_memory
    
  subroutine print_total_memory
    implicit none
    
    if (tot_memory > 1024.**3 ) then
       write(*,'(A,f6.2,A)') "Aggregate memory: ",tot_memory/(1024.**3)," GB"
    else if (tot_memory > 1024.**2 ) then   
       write(*,'(A,f6.2,A)') "Aggregate memory: ",tot_memory/(1024.**2)," MB"
    else if (tot_memory > 1024. ) then   
       write(*,'(A,f6.2,A)') "Aggregate memory: ",tot_memory/(1024.)," KB"
    end if
    print*
  end subroutine print_total_memory
!!!===========================================================
!!!===========================================================
  subroutine read_input_file(finput,spfile,intfile,denfile,reffile,lawfile,obsfiles)
    implicit none

    character(200) :: spfile,intfile,finput,reffile,denfile,lawfile
    character(200),allocatable,dimension(:) :: obsfiles
    integer :: AA,Aprot,Aneut,ii

    call getenv("MRCI_SP_FILES",SP_DIR)
    SP_DIR=adjustl(SP_DIR)
    call getenv("MRCI_INIFILES",INI_DIR)
    INI_DIR=adjustl(INI_DIR)
    call getenv("MRCI_OUTPUT",OUTPUT_DIR)
    OUTPUT_DIR=adjustl(OUTPUT_DIR)
    call getenv("MRCI_ME_FILES",ME_DIR)
    ME_DIR=adjustl(ME_DIR)

   
    call check_for_file(INI_DIR,finput)
    finput = adjustl(finput)    
    if (finput(1:12) == "../inifiles/") then
       finput = finput(13:200) 
       finput = adjustl(finput)
    end if
    
    open(unit=45,file=trim(INI_DIR)//trim(finput)) 

    read(45,*) !!!Enter number of nucleons (Z,N) 
    read(45,*) mbas%Aprot,mbas%Aneut
    mbas%Abody = mbas%Aneut + mbas%Aprot

    read(45,*) !!!Enter target state nucleons (Z,N) 
    read(45,*) mbas%ztarg , mbas%ntarg
    mbas%Atarg = mbas%Ztarg + mbas%Ntarg

    read(45,*) !!!Enter Parity (0,1) 
    read(45,*) mbas%ptarg

    read(45,*) !!! Enter SP file
    read(45,*) spfile
    call check_for_file(SP_DIR,spfile)
    
    read(45,*) !!! Enter INT file
    read(45,*) intfile
    call check_for_file(ME_DIR,intfile)

    read(45,*) !!! Enter INT file
    read(45,*) denfile
    call check_for_file(ME_DIR,denfile)
    
    read(45,*) !Enter REF file
    read(45,*) reffile 
    call check_for_file(INI_DIR,reffile)
    
    read(45,*) !Enter lawson param
    read(45,*) mbas%law_beta

    read(45,*) !Enter lawson file
    read(45,*) lawfile

    if (non_negligible(mbas%law_beta)) then
       call check_for_file(ME_DIR,lawfile)
    end if
       
    read(45,*) !Enter number of observables
    read(45,*) mbas%num_obs

    allocate(obsfiles(mbas%num_obs))
    read(45,*) !Enter observable files
    
    do ii = 1, mbas%num_obs
       read(45,*)  obsfiles(ii)
       call check_for_file(ME_DIR,obsfiles(ii))
    end do

    
    ii = 1
    do while (.true.) 
       if (finput(ii:ii+3) == ".ini") exit
       if (ii==100) STOP "ini file is confusing me"
       ii = ii +1
    end do
    prefix = finput(1:ii-1)

    ii = 1
    do while (.true.) 
       if (intfile(ii:ii+4) == "hwHO0") exit
       if (ii==100) STOP "ini file is confusing me"
       ii = ii +1
    end do

    read(intfile(ii+5:ii+6) ,"(I2)") hw

    close(45)

    mbas%mtarg = mod(mbas%Atarg,2)    

    call dcgi00()

    write(*,"(A)") "Target Nucleus: "//nucleus_name(mbas%Ntarg,mbas%Ztarg)
    write(*,"(A)") "Reference Nucleus: "//nucleus_name(mbas%Aneut,mbas%Aprot)
    
  end subroutine read_input_file
!!!===========================================================
!!!===========================================================
  subroutine read_ref_file(reffile,num_refs,REF)
    implicit none

    character(200) :: reffile
    integer :: ist,num_refs,AA,ii
    integer,allocatable,dimension(:,:) :: ref

    AA = mbas%Abody

    reffile = adjustl(reffile)
    open(unit=45,file=trim(INI_DIR)//trim(reffile))     
    
    read(45,*) !! Enter number of refs
    read(45,*) num_refs
    
    allocate(ref(num_refs,AA))

    
    read(45,*) !!Enter the reference states, one by one
    do ii = 1, num_refs
       read(45,*,iostat=ist) ref(ii,:)
       if (ist < 0) STOP "REF FILE TOO SHORT"
    end do
    
    close(45)
  end subroutine read_ref_file
!!!===========================================================
!!!===========================================================
  character(5) function nucleus_name(N,Z) 
!!! print out the spectral notation for an sp state
    implicit none
    
    integer,intent(in) :: N,Z
    character(2),dimension(0:50) :: names 
    character(2) :: Astr
    
    names = (/"n ","H ","He","Li","Be","B ","C ","N ","O ","F ",&
         "Ne","Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca",&
         "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",&
         "Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo",&
         "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn"/)
        
    write(Astr,"(I2)") N+Z 
    nucleus_name =  adjustl(trim(adjustl(names(Z)))//trim(adjustl(Astr)))    
    
  end function nucleus_name
!!!===========================================================
!!!===========================================================
  character(6) function spec_not(t,n,j,l) 
!!! print out the spectral notation for an sp state
    implicit none
    
    integer,intent(in) :: t,j,l,n
    character(1),dimension(6) :: l_let 
    character(1) :: j_str,n_str,l_str
    
    l_let = (/'s','p','d','f','g','h'/) 
    
    write(n_str,'(I1)') n
    write(j_str,'(I1)') j 
    
    if (l > 5) then 
       l_str = 'x' 
    else
       l_str = l_let(l+1)
    end if

    if (t==-1) then
       spec_not = 'p'//n_str//l_str//j_str//'/2'
    else
       spec_not = 'n'//n_str//l_str//j_str//'/2' 
    end if
  end function spec_not
!!!===========================================================
!!!===========================================================    
  subroutine parity_M_and_T(state,mbas,PAR,M,T)
    implicit none
    
    integer,dimension(:) :: state
    integer :: ii,smp,smm,PAR,M,T,smt
    type(spd) :: mbas
    
    smp = 0
    smm = 0 
    smt = 0
    do ii = 1,size(state)
       smp = smp + mbas%ll(state(ii))
       smm = smm + mbas%mm(state(ii))
       smt = smt + mbas%tz(state(ii))
    end do

    PAR = mod(smp,2)
    M = smm
    T = smt
  end subroutine parity_M_and_T
!!!===========================================================
!!!===========================================================
  integer function digets(num)
    implicit none

    integer,intent(in) :: num


    if (num<=1) then
       digets = 1
    else
       digets = ceiling(log10(float(num)))
    end if
  end function digets

  character(10) function fmtlen(num)
    integer,intent(in) :: num
    character(10) :: fmt

    write(fmt,'(I10)') num
    fmtlen=adjustl(fmt)

  end function fmtlen

  subroutine print_matrix(matrix)
    implicit none 

    integer :: i,m1,m2 
    real(8),dimension(:,:) :: matrix
    character(1) :: y
    character(10) :: fmt2

    m1=size(matrix(1,:))
    m2=size(matrix(:,1)) 

    write(y,'(i1)') m1

    fmt2= '('//y//'(f14.8))'	

    print*
    do i=1,m2
       write(*,fmt2) matrix(i,:)
    end do
    print* 

  end subroutine print_matrix
!=====================================================
integer function bosonic_tp_index(i,j,n)
  ! n is total number of sp states
  ! assume i <= j 
  implicit none 
  
  integer :: i,j,n
  
  bosonic_tp_index = n*(i-1) + (3*i-i*i)/2 + j - i 
  
end function 
!========================================================
logical function instr(phrase,input)
  implicit none

  character(*),intent(in) :: input,phrase
  integer :: len_i,len_p,q 

  len_i = len(input)
  len_p = len(phrase)

  do q = 1 , len_i
     if (input(q:q+len_p-1) == phrase ) then
        instr = .true.
        return
     end if
  end do
  instr=.false. 
  
end function instr
!========================================================
logical function non_negligible(X)
  implicit none

  real(8),intent(in) :: X 

  if (abs(X) > 1e-8 ) then
     non_negligible = .true.
  else
     non_negligible = .false.
  end if

end function non_negligible
!=========================================================
subroutine check_for_file(DIR,file)
  implicit none

  character(*) :: DIR,file
  logical :: ex

  inquire( file = trim(adjustl(DIR))//trim(adjustl(file)), exist=ex)

  if (.not. ex) then
     write(*,"(A)") "Could not find file: "// &
          trim(adjustl(DIR))//trim(adjustl(file)) 
     stop
  end if
end subroutine check_for_file
!==================================================================
!==================================================================
subroutine write_results(e,s,v)
  implicit none

  real(8),dimension(:) :: e,s
  real(8),dimension(:,:) :: v
  integer :: Nstates,dim

  Nstates = size(e)
  dim = size(v(1,:))
  print* 
  write(*,"(A)") "Writing Eigensystem to "//"../output/"//trim(adjustl(prefix))//&
       "_eigensystem.dat"   
  open(unit=39,file="../output/"//trim(adjustl(prefix))//&
       "_eigensystem.dat",form="unformatted")

  write(39) Nstates, dim
  write(39) e
  write(39) s
  write(39) v
  close(39)
  
end subroutine write_results
!==================================================================
!==================================================================
subroutine read_results(e,s,v)
  implicit none

  real(8),allocatable,dimension(:) :: e,s
  real(8),allocatable,dimension(:,:) :: v
  integer :: Nstates,dim



  print* 
  write(*,"(A)") "Reading Eigensystem from "//"../output/"//trim(adjustl(prefix))//&
       "_eigensystem.dat"   
  open(unit=39,file="../output/"//trim(adjustl(prefix))//&
       "_eigensystem.dat",form="unformatted")

  read(39) Nstates, dim

  allocate(e(Nstates),s(Nstates))
  allocate(v(Nstates,dim))

  
  read(39) e
  read(39) s
  read(39) v
  close(39)
  
end subroutine read_results

end module mrci_aux
  
  
