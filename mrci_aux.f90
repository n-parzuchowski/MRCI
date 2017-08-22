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
   real(8) :: e
   integer, allocatable, dimension(:) :: nn, ll, jj, mm,tz,jlab
   type(int_vec),allocatable,dimension(:) :: amap,qmap
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
character(500) :: ME_DIR,SP_DIR,INI_DIR,OUTPUT_DIR
     
contains
!!!===========================================================
!!!===========================================================
  subroutine read_input_file(finput,spfile,intfile,reffile,AA,Aprot,Aneut)
    implicit none

    character(200) :: spfile,intfile,finput,reffile 
    integer :: AA,Aprot,Aneut

    call getenv("MRCI_SP_FILES",SP_DIR)
    SP_DIR=adjustl(SP_DIR)
    call getenv("MRCI_INIFILES",INI_DIR)
    INI_DIR=adjustl(INI_DIR)
    call getenv("MRCI_OUTPUT",OUTPUT_DIR)
    OUTPUT_DIR=adjustl(OUTPUT_DIR)
    call getenv("MRCI_ME_FILES",ME_DIR)
    ME_DIR=adjustl(ME_DIR)
    
    finput = adjustl(finput)
    open(unit=45,file=trim(INI_DIR)//trim(finput)) 

    read(45,*) !!! Enter SP filee
    read(45,*) spfile
    
    read(45,*) !!! Enter INT file
    read(45,*) intfile

    read(45,*) !Enter REF file
    read(45,*) reffile 

    read(45,*) !!!Enter number of nucleons (A,Z,N) 
    read(45,*) AA,Aprot,Aneut
   
    close(45)
    
    call dcgi00()

  end subroutine read_input_file
!!!===========================================================
!!!===========================================================
  subroutine read_ref_file(reffile,num_refs,REF,AA,denfile)
    implicit none

    character(200) :: reffile,denfile
    integer :: ist,num_refs,AA,ii
    integer,allocatable,dimension(:,:) :: ref

    reffile = adjustl(reffile)
    open(unit=45,file=trim(INI_DIR)//trim(reffile)) 

    read(45,*) !! Enter number of refs
    read(45,*) denfile

    read(45,*) !! Enter number of refs
    read(45,*) num_refs
    
    allocate(ref(num_refs,AA))

    
    read(45,*) !!Enter the reference states, one by one
    do ii = 1, num_refs
       read(45,*) ref(ii,:)
    end do
    
    close(45)
  end subroutine read_ref_file
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

  
  end module mrci_aux
