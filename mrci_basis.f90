module mrci_basis
  use mrci_aux
  implicit none
!!!===========================================================
!!!===========================================================

contains
!!!===========================================================
!!!===========================================================
  subroutine init_sp_basis(fname)
    !! basis descriptors for both the m-scheme and j-scheme
    !! bases. The latter is needed because that's how we store
    !! the interaction.    
    implicit none

    character(200) :: fname
    integer :: a,ii,ist,n,l,j,t,m
    real(8) :: e


    fname = adjustl(fname)

    open(unit=77,file=trim(SP_DIR)//trim(fname))

    mbas%Ntot = 0
    do
       read(77,*,iostat = ist) a,n,l,j,t,e

       if ( ist > 0 ) stop "Input Error: .sps file"
       if (ist < 0 ) exit

       mbas%Ntot = mbas%Ntot + j + 1       
       
    end do

    rewind(77) 

    mbas%Njjcoup = a 
    allocate(mbas%nn(mbas%Ntot))
    allocate(mbas%ll(mbas%Ntot))
    allocate(mbas%jj(mbas%Ntot))
    allocate(mbas%tz(mbas%Ntot))
    allocate(mbas%mm(mbas%Ntot))
    allocate(mbas%jlab(mbas%Ntot))

    jbas%Ntot=a
    allocate(jbas%nn(a))
    allocate(jbas%ll(a))
    allocate(jbas%jj(a))
    allocate(jbas%tz(a))
    allocate(jbas%jlab(a))

    ii = 1 
    do 
       read(77,*,iostat = ist) a,n,l,j,t,e

       if ( ist > 0 ) stop "Input Error: .sps file"
       if (ist < 0 ) exit

       mbas%nn(ii:ii+j) = n 
       mbas%ll(ii:ii+j) = l
       mbas%jj(ii:ii+j) = j
       mbas%tz(ii:ii+j) = t
       mbas%jlab(ii:ii+j) = a

       do  m = -j,j,2
          mbas%mm(ii+(m+j)/2) = m
       end do
       ii = ii + j + 1

       jbas%nn(a) = n 
       jbas%ll(a) = l
       jbas%jj(a) = j
       jbas%tz(a) = t
    end do
    
    mbas%jmax = maxval(mbas%jj)
    mbas%lmax = maxval(mbas%ll)
    mbas%nmax = maxval(mbas%nn)

    jbas%jmax = maxval(jbas%jj)
    jbas%lmax = maxval(jbas%ll)
    jbas%nmax = maxval(jbas%nn)
    jbas%Abody = mbas%Abody
    jbas%Aneut = mbas%Aneut
    jbas%Aprot = mbas%Aprot
    jbas%Ptarg = mbas%Ptarg
    jbas%Mtarg = mbas%Mtarg
    jbas%dTz = mbas%dtz
    jbas%Atarg = mbas%Atarg
    jbas%Ztarg = mbas%Ztarg
    jbas%Ntarg = mbas%Ntarg

    close(77)
  end subroutine init_sp_basis
!!!===========================================================
!!!===========================================================
  subroutine generate_basis(REF,BASIS)
    !! this subroutine generates the reference states
    !! based on the one-body density matrix
    implicit none

    integer,allocatable,dimension(:,:) :: SD_BASIS,BASIS
    integer,dimension(:,:) :: REF
    integer :: Abody,ix,jx,kx,lx,num_refs,PAR,M,jin,lout,tout
    integer :: q,mout,nin,test,BigT
    integer :: PTarg,MTarg,dTz
    integer,allocatable,dimension(:) :: valid
    integer,dimension(mbas%Abody) :: newSD 
    real(8) :: t1,t2,omp_get_wtime
    logical :: present

    t1=omp_get_wtime()
    write(*,'(A)')  '===================================================='
    write(*,'(A)')  'GENERATING SD BASIS'
    write(*,*)

    
    Abody = mbas%Abody
    Mtarg = mbas%Mtarg
    Ptarg = mbas%Ptarg
    dTz = mbas%dTz
    
    num_Refs = size(REF(:,1)) 
    allocate(SD_BASIS(100000,Abody)) 
    write(*,'(A)') "Allocated generic SD basis placeholder"
    call print_memory(100000*Abody*4.d0)
    tot_memory = tot_memory + 100000*Abody*4.d0
    call print_total_memory
    
    
    do ix = 1, num_Refs
       call parity_M_and_T(REF(ix,:),mbas,PAR,M,BigT ) 
       if( ( PAR .ne. 0 ) .or. (M .ne. proj) .or.&
            (BigT .ne. (mbas%Aneut-mbas%Aprot)))  then
          print*, ix,BigT,mbas%Aneut-mbas%Aprot,PAR,M
          STOP "REFERENCES HAVE BAD QUANTUM NUMBERS."
       end if
       call sort_SD(REF(ix,:))
    end do

    !!!  NOW Systematically go through the references and construct all possible MR-CIS excitations
    if((MTarg == 0).and.(PTarg==0).and.(dTz==0)) then
       !!! same quantum numbers as ground state
       SD_BASIS(1:num_refs,:) = REF 
       q = num_refs+1
    else
       q = 1
    end if

    do ix = 1, num_refs
       do jx = 1, Abody 

          lout = mbas%ll(REF(ix,jx))
          mout = mbas%mm(REF(ix,jx))
          tout = mbas%tz(REF(ix,jx))
          
          do kx = 1, mbas%Ntot

             !! first check if this state has the right symmetry 
             if (mbas%mm(kx) -  mout .ne. MTarg) cycle
             if (mod(mbas%ll(kx)+lout,2).ne. PTarg) cycle
             if (mbas%tz(kx) - tout .ne. dTz) cycle


             !! IF we've made it here, the (kx,-jx) state is a candidate for a 1p1h excitation.
             !! make sure kx isn't already in the current reference
             nin = mbas%nn(kx)
             jin = mbas%jj(kx)

             
             present = .false. 
             do lx = 1, Abody
       
                if (ref(ix,lx) .ne. kx) cycle

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
    write(*,'(A)')  "Allocated unique SD basis descriptor"
    call print_memory(q*Abody*4.d0)
    tot_memory = tot_memory + q*Abody*4.d0
    call print_total_memory

    BASIS = SD_BASIS(1:q,:)

    deallocate(SD_BASIS)
    write(*,'(A)')  "Deallocated generic SD basis placeholder"
    tot_memory = tot_memory - 100000*Abody*4.d0
    call print_total_memory

    write(*,'(A)')  '===================================================='
    write(*,'(A)')  'SLATER DETERMINANT BASIS GENERATED'
    write(*,'(I5,A)')   q, ' BASIS VECTORS'

    write(*,'(A)')  "In order to store a matrix in this basis:"
    call print_memory(dfloat(q)*dfloat(q+1)*4.d0)
    
    
    !! check that nothing is wrong
    do ix = 1,q
       call parity_M_and_T(BASIS(ix,:),mbas,PAR,M,BigT ) 
       if( ( PAR .ne. PTarg ) .or. (M .ne. MTarg) .or. &
            (BigT .ne. (mbas%Aneut-mbas%Aprot+dTz)))  then
          STOP "BASIS HAS BAD QUANTUM NUMBERS."
       end if
    end do
    t2=omp_get_wtime()
    write(*,"(A,f6.1,A,f10.1)") "Time: ", t2-t1, " Total: ", t2-time_Zero 
  end subroutine generate_basis
!!!===========================================================
!!!===========================================================
  subroutine sort_SD(ar)
    ! dumb bubble sort
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
!!!===========================================================
!!!===========================================================                  
  subroutine generate_tp_basis 
    implicit none

    integer :: Ntot,Lmax,eMax,lj,twol,twoj,j_min,j_max,numJ
    integer :: ljMax,q,idx,idxx,i,j,a,Tz,Pi,JT,x
    integer :: n1,j1,l1,t1,n2,j2,l2,t2,lj1,lj2,n
    integer,allocatable,dimension(:,:) :: SPBljs
    integer,allocatable,dimension(:) :: nMax_lj
    integer,dimension(500) :: a_list



    a_list = 0
!!! single particle ordering scheme from Heiko's code 
    Ntot = jbas%Ntot   

    Lmax = jbas%lmax
    eMax = 2*jbas%nMax
    ! populate lj array
    lj = 0
    do twol = 0, 2 * Lmax , 2
       do  twoj = abs(twol - 1) , twol+1 , 2
          lj=lj+1
       end do
    end do
    ljMax = lj 
    allocate(SPBljs(lj,2)) 
    allocate(nMax_lj(lj))

    allocate(jbas%tlj_to_ab(2,ljmax))
    
    lj = 0
    do twol = 0, 2 * Lmax , 2
       do  twoj = abs(twol - 1) , twol+1 , 2
          lj=lj+1
          SPBljs(lj,1) = twol
          sPBljs(lj,2) = twoj
          nMax_lj(lj) = (eMax - twol/2)/2

          allocate(jbas%tlj_to_ab(1,lj)%Z(nMax_lj(lj)+1)) 
          allocate(jbas%tlj_to_ab(2,lj)%Z(nMax_lj(lj)+1))          

          do n = 0 , nMax_lj(lj)

             ! now search for sp labels
             do i = 1, jbas%Ntot  !!! proton 
                if ( jbas%jj(i) .ne. twoj ) cycle
                if ( jbas%nn(i) .ne. n ) cycle
                if ( 2*jbas%ll(i) .ne. twol ) cycle
                if ( jbas%tz(i) .ne. -1 ) cycle                     
                ! i is the jlabel for this state
                jbas%tlj_to_ab(2,lj)%Z(n+1) = i 
                exit
             end do
             
             do i = 1, jbas%Ntot  !!! proton 
                if ( jbas%jj(i) .ne. twoj ) cycle
                if ( jbas%nn(i) .ne. n ) cycle
                if ( 2*jbas%ll(i) .ne. twol ) cycle
                if ( jbas%tz(i) .ne. 1 ) cycle                     
                ! i is the jlabel for this state
                jbas%tlj_to_ab(1,lj)%Z(n+1) = i 
                exit
             end do             
             
          end do

          
       end do
    end do

    
    allocate(jbas%amap(Ntot*(Ntot+1)/2))
    allocate(jbas%qmap(Ntot*(Ntot+1)/2)) 
    ! v_elem to find matrix elements
    do i = 1,Ntot
       do j = i,Ntot
          
          j_min = abs(jbas%jj(i)-jbas%jj(j))
          j_max = jbas%jj(i) + jbas%jj(j) 
          
          numJ = (j_max - j_min)/2 + 2
          
          x = bosonic_tp_index(i,j,Ntot) 
          allocate(jbas%amap(x)%Z(numJ))
          allocate(jbas%qmap(x)%Z(numJ)) 
          jbas%amap(x)%Z = 0
          jbas%amap(x)%Z(1) = j_min
          jbas%qmap(x)%Z = 0
       end do
    end do
        
    

    q = 0
    ! count the blocks
   
    do Tz = 1 , -1, -1  
       do Pi = 0,1
          do JT = 0, 2*jbas%Jmax,2 
                    
             ! count states in the block
             a = 0
             
             do lj1 = 1, ljMax
                do lj2 = 1, ljMax
                   
                   j1 = SPBljs(lj1,2) 
                   j2 = SPBljs(lj2,2)
                   l1 = SPBljs(lj1,1)/2
                   l2 = SPBljs(lj2,1)/2
                   
                   if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
                   if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 
                                      
                   do n1 = 0,nMax_lj(lj1)
                      idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                      do n2 = 0,nMax_lj(lj2) 
                         idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                         
                         if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                         if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                              (n1==n2) .and. (Tz .ne. 0) ) cycle
                         
                         a = a + 1
                         
                      end do
                   end do
                end do
             end do

             if (a ==0 ) cycle
             ! we've got a block here
             q = q + 1
             a_list(q) = a

          end do
       end do
    end do

    tp_basis%bMax = q
    tp_basis%aMaxMax = maxval(a_list)    

    allocate(tp_basis%block(q))   

    do q = 1, tp_basis%bMax
       tp_basis%block(q)%aMax = a_list(q) 
       allocate(tp_basis%block(q)%qnums(a_list(q),2))
    end do
    ! okay let's restart

    q = 0
    ! heiko's code calls protons 1 and neutrons 0
    
    
    do Tz = 1 , -1, -1  
       do Pi = 0,1
          do JT = 0, 2*jbas%Jmax,2 

             q = q + 1
             select case ( Tz)
             case ( -1 ) 
                t1 = -1
                t2 = -1
             case ( 0 ) 
                t1 = 1 
                t2 = -1
             case ( 1 ) 
                t1 = 1
                t2 = 1 
             end select

             a = 0

             do lj1 = 1, ljMax
                do lj2 = 1, ljMax

                   j1 = SPBljs(lj1,2) 
                   j2 = SPBljs(lj2,2)
                   l1 = SPBljs(lj1,1)/2
                   l2 = SPBljs(lj2,1)/2

                   if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
                   if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 

                   do n1 = 0,nMax_lj(lj1)
                      idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                      do n2 = 0,nMax_lj(lj2) 
                         idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                         
                         if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                         if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                              (n1==n2) .and. (Tz .ne. 0) ) cycle

                         ! now search for sp labels
                         do i = 1, jbas%Ntot 
                            if ( jbas%jj(i) .ne. j1 ) cycle
                            if ( jbas%nn(i) .ne. n1 ) cycle
                            if ( jbas%ll(i) .ne. l1 ) cycle
                            if ( jbas%tz(i) .ne. t1 ) cycle                     
                            exit
                         end do

                         do j = 1, jbas%Ntot 
                            if ( jbas%jj(j) .ne. j2 ) cycle
                            if ( jbas%nn(j) .ne. n2 ) cycle
                            if ( jbas%ll(j) .ne. l2 ) cycle
                            if ( jbas%tz(j) .ne. t2 ) cycle                     
                            exit
                         end do


                         a = a + 1
                         
                         tp_basis%block(q)%qnums(a,1) = i 
                         tp_basis%block(q)%qnums(a,2) = j

                         if (j < i ) then
                            x = bosonic_tp_index(j,i,Ntot) 
                            j_min = jbas%amap(x)%Z(1)
                            jbas%amap(x)%Z((JT-j_min)/2+2) = a
                            !!! figure out how this fucking basis is working.
                         else
                            x = bosonic_tp_index(i,j,Ntot) 
                            j_min = jbas%amap(x)%Z(1)
                            jbas%amap(x)%Z((JT-j_min)/2+2) = a
                         end if
                         
                         jbas%qmap(x)%Z((JT-j_min)/2+2) = q                         
                      end do
                   end do
                end do
             end do

             if  (a == 0 ) then
                ! this wasn't a real block...
                q = q -1
             else
                tp_basis%block(q)%J = JT
                tp_basis%block(q)%PAR = Pi
                tp_basis%block(q)%Tz = Tz
             end if
                
             
          end do
       end do
    end do
           
    
  end subroutine generate_tp_basis


  integer function get_tp_block_index(i,j,JT)
    implicit none

    integer :: i,j,JT,j_min,Ntot,x

    Ntot = jbas%Ntot    
    x = bosonic_tp_index(i,j,Ntot)     
    j_min = jbas%amap(x)%Z(1) 
    
    get_tp_block_index=jbas%qmap(x)%Z((JT-j_min)/2+2) 

  end function get_tp_block_index

  integer function TP_index(i,j,JT)
    implicit none
    
    integer :: i,j,JT,j_min,Ntot,x
    
    Ntot = jbas%Ntot
    x = bosonic_tp_index(i,j,Ntot) 
    j_min = jbas%amap(x)%Z(1)

    TP_index= jbas%amap(x)%Z((JT-j_min)/2+2) 
  
  end function TP_index
    
end module mrci_basis
