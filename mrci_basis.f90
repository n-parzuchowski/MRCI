module mrci_basis
  implicit none

  TYPE :: spd    ! single particle discriptor
     integer :: Ntot, Jmax, Lmax,Njjcoup,nMax
     integer :: n,l,j,t,m,Abody,Aneut,Aprot
     real(8) :: e
     integer, allocatable, dimension(:) :: nn, ll, jj, mm,tz,jlab
     ! for clarity:  nn, ll, nshell are all the true value
     ! jj is j+1/2 (so it's an integer) 
     ! likewise itzp is 2*tz
   contains
     procedure :: build => init_sp_basis 
  END TYPE spd

  
contains

  subroutine init_sp_basis(this,fname)
    implicit none
    class(spd),intent(inout) :: this
    character(200) :: fname
    integer :: a,ii,ist,n,l,j,t,m
    real(8) :: e


    fname = adjustl(fname)

    open(unit=77,file=trim(fname))

    this%Ntot = 0
    do
       read(77,*,iostat = ist) a,n,l,j,t,e

       if ( ist > 0 ) stop "Input Error: .sps file"
       if (ist < 0 ) exit

       this%Ntot = this%Ntot + j + 1       
       
    end do

    rewind(77) 

    this%Njjcoup = a 
    allocate(this%nn(this%Ntot))
    allocate(this%ll(this%Ntot))
    allocate(this%jj(this%Ntot))
    allocate(this%tz(this%Ntot))
    allocate(this%mm(this%Ntot))
    allocate(this%jlab(this%Ntot))
    ii = 1 
    do 
       read(77,*,iostat = ist) a,n,l,j,t,e

       if ( ist > 0 ) stop "Input Error: .sps file"
       if (ist < 0 ) exit

       this%nn(ii:ii+j) = n 
       this%ll(ii:ii+j) = l
       this%jj(ii:ii+j) = j
       this%tz(ii:ii+j) = t
       this%jlab(ii:ii+j) = a

       do  m = -j,j,2
          this%mm(ii+(m+j)/2) = m
       end do
        ii = ii + j + 1
    end do
    
    this%jmax = maxval(this%jj)
    this%lmax = maxval(this%ll)
    this%nmax = maxval(this%nn)
    
    close(77)
  end subroutine init_sp_basis



  subroutine generate_basis(jbas,REF,BASIS)
    !! this subroutine generates the reference states
    !! based on the one-body density matrix
    implicit none

    type(spd) :: jbas
    integer,allocatable,dimension(:,:) :: SD_BASIS,BASIS
    integer,dimension(:,:) :: REF
    integer :: Abody,ix,jx,kx,lx,num_refs,PAR,M,jin,lout,tout
    integer :: q,mout,nin,test,BigT
    integer,allocatable,dimension(:) :: valid
    integer,dimension(jbas%Abody) :: newSD 
    logical :: present

    Abody = jbas%Abody

    num_Refs = size(REF(:,1)) 
    allocate(SD_BASIS(14000,Abody)) 

    do ix = 1, num_Refs
       call parity_M_and_T(REF(ix,:),jbas,PAR,M,BigT ) 
       if( ( PAR .ne. 0 ) .or. (M .ne. 0) .or.&
            (BigT .ne. (jbas%Aneut-jbas%Aprot)))  then
          print*, ix,BigT,jbas%Aneut-jbas%Aprot
          STOP "REFERENCES HAVE BAD QUANTUM NUMBERS."
       end if
       call sort_SD(REF(ix,:))
    end do

    !!!  NOW Systematically go through the references and construct all possible MR-CIS excitations
    SD_BASIS(1:num_refs,:) = REF 

    q = num_refs+1
    do ix = 1, num_refs
       do jx = 1, Abody 

          lout = jbas%ll(REF(ix,jx))
          mout = jbas%mm(REF(ix,jx))
          tout = jbas%tz(REF(ix,jx))
          
          do kx = 1, jbas%Ntot

             !! first check if this state has the right  
             if (jbas%mm(kx) .ne. mout) cycle
             if (mod(jbas%ll(kx)+lout,2).ne.0) cycle
             if (jbas%tz(kx) .ne. tout) cycle


             !! IF we've made it here, this state is a candidate for a 1p1h excitation.
             !! Now we see if it's in the current SD
             nin = jbas%nn(kx)
             jin = jbas%jj(kx)

             
             present = .false. 
             do lx = 1, Abody
                if (jbas%mm(lx) .ne. mout) cycle
                if (jbas%ll(lx).ne. lout) cycle
                if (jbas%tz(lx) .ne. tout) cycle
                if (jbas%jj(lx).ne. jin) cycle
                if (jbas%nn(lx) .ne. nin) cycle
                !if we're here, this means that this state is already in the SD
                present = .true.
                exit
             end do

             if (present) cycle ! don't count this state if it's already in the SD 

             !! if we're here, that means we've found a particle state that isn't in the current SD                                                                  
             newSD = REF(ix,:)
             newSD(jx) = kx

             call sort_SD(newSD) 
             

             ! check that this SD isn't already in SD_Basis

             present = .false. 
             do lx = 1, q-1                
                test = sum( abs( newSD - SD_Basis(lx,:)) )
                if (test == 0 ) then                
                   present = .true.
                   exit
                end if
             end do

             if (present) cycle

             !! if we've made it here, we have a shiny new slater determinant

             SD_Basis(q,: )= newSD
             q = q +1
             

          end do
       end do
    end do

    !! pack 'er up
    q= q-1
    allocate(BASIS(q,Abody))
    BASIS = SD_BASIS(1:q,:)

    deallocate(SD_BASIS)

    PRINT*, 'SLATER DETERMINANT BASIS GENERATED'
    PRINT*, q, 'BASIS VECTORS'
    
    !! check that nothing is wrong
    do ix = 1,q
       call parity_M_and_T(BASIS(ix,:),jbas,PAR,M,BigT ) 
       if( ( PAR .ne. 0 ) .or. (M .ne. 0) .or. &
            (BigT .ne. (jbas%Aneut-jbas%Aprot)))  then
          STOP "BASIS HAS BAD QUANTUM NUMBERS."
       end if
    end do

    
  end subroutine generate_basis

  subroutine sort_SD(ar)
    implicit none
    
    integer,dimension(:) :: ar
    integer :: n,newn,x,y ,i,xx
    n =  size(ar)

    do while (n > 1)
       newn = 1
       do i = 2, n
          if (ar(i-1) > ar(i) ) then
             x = ar(i)
             y = ar(i-1)
             ar(i-1) = x
             ar(i) = y
             newn = i
          end if
       end do
       n = newn
    end do
  end subroutine sort_SD
    
    
    
  subroutine parity_M_and_T(state,jbas,PAR,M,T)
    implicit none
    
    integer,dimension(:) :: state
    integer :: ii,smp,smm,PAR,M,T,smt
    type(spd) :: jbas
    
    smp = 0
    smm = 0 
    smt = 0
    do ii = 1,size(state)
       smp = smp + jbas%ll(state(ii))
       smm = smm + jbas%mm(state(ii))
       smt = smt + jbas%tz(state(ii))
    end do

    PAR = mod(smp,2)
    M = smm
    T = smt
  end subroutine parity_M_and_T
    
    
    
  character(6) function spec_not(t,n,j,l) 
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

 
          
          
          

    

    
    

  

  
end module 
