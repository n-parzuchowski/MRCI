module mrci_solve
  use mrci_interactions
  use mrci_operators
  implicit none

  real(8),allocatable,dimension(:) :: HAM,JTOT_MAT
contains

  subroutine diagonalize(z0,z1,z2,basis,Eimsrg)
    implicit none
    
    type(block_mat_full),allocatable,dimension(:,:) :: z1    
    type(block_mat),allocatable,dimension(:) :: z2
    integer,dimension(:,:) :: basis
    real(8) :: z0,t3,t4,t5,t6,Eimsrg
    real(8),allocatable,dimension(:) :: workl,DX,QX,resid,work,workD
    real(8),allocatable,dimension(:,:) :: V,Z
    integer :: lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),dm,x
    integer :: ishift,mxiter,nconv,mode,lworkl,ldz,nev,inc,ii,jj,aa
    real(8) :: tol,sigma,sm,Egs,t1,t2,omp_get_wtime,amp1,amp2,dcgi,spin
    character(1) :: BMAT,HOWMNY 
    character(2) :: which
    logical :: rvec
    logical,allocatable,dimension(:) :: selct

    t1 = omp_get_wtime()

    print* 
    write(*,"(A)") '=====================================' 
    write(*,"(A)") 'SOLVING MRCI EIGENVALUE PROBLEM'
    print*
    
    dm = size(basis(:,1))
    allocate(HAM(dm*(dm+1)/2))
    HAM=0.d0

    write(*,"(A)") "Allocated Hamiltonian storage" 
    call print_memory(dm*(dm+1)*4.d0)
    tot_memory = tot_memory + dm*(dm+1)*4.d0

    print*
    allocate(Jtot_mat(dm*(dm+1)/2))
    Jtot_mat=0.d0

    write(*,"(A)") "Allocated J storage" 
    call print_memory(dm*(dm+1)*4.d0)
    tot_memory = tot_memory + dm*(dm+1)*4.d0

    call print_total_memory
    
    nev = 10 ! I only care about the ground state right now. 
    ido = 0  ! status integer is 0 at start
    BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
    which = 'SA' ! compute smallest eigenvalues (algebraic) ('SM') is magnitude.
    tol = 1.0E-10 ! error tolerance? (wtf zero?) 
    info = 0
    ncv = 5*nev ! number of lanczos vectors I guess
    lworkl = ncv*(ncv+8) 
    allocate(V(dm,NCV),workl(lworkl))
    LDV = dm  
    ishift = 1
    mxiter = 500 
    mode = 1
    
    allocate(resid(dm),work(10*dm),workD(3*dm)) 


    write(*,"(A)") "Allocated workspace" 
    call print_memory((14*dm+ dm*NCV + lworkl)*8.d0)
    tot_memory = tot_memory + (14*dm+ dm*NCV + lworkl)*8.d0

    call print_total_memory

    
    iparam(1) = ishift
    iparam(3) = mxiter
    iparam(7) = mode
    ii = 0

    t3 = omp_get_Wtime()
    do 
       ! so V is the krylov subspace matrix that is being diagonalized
       ! it does not need to be initialized, so long as you have the other 
       ! stuff declared right, the code should know this. 
       call dsaupd ( ido, bmat, dm, which, nev, tol, resid, &
            ncv, v, ldv, iparam, ipntr, workd, workl, &
            lworkl, info )
       ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
       
       call progress_bar( ii )
       ii=ii+1
       
       
       if ( ido /= -1 .and. ido /= 1 ) then
          exit
       end if

       
       if (ii ==1 )then
          call first_mat_vec_prod(z0,z1,z2,workd(ipntr(1)),workd(ipntr(2)),basis,dm)
       else
          call mat_vec_prod(z0,z1,z2,workd(ipntr(1)),workd(ipntr(2)),basis,dm)
       end if
       
    end do
    t4 = omp_get_Wtime()
    write(6,*) 
    write(*,"(A,I5,A,f10.1,A)")  "converged after", ii, " iterations and "&
         ,t4-t3," seconds" 
    print*
    
    ! the ritz values are out of order right now. Need to do post
    ! processing to fix this, and get the eigenvectors
    rvec= .true. 
    howmny = 'A'
    
    allocate(selct(NCV)) 
    allocate(DX(NEV)) 
    allocate(Z(dm,NEV)) 
    ldz = dm  
    call dseupd( rvec, howmny, selct, DX, Z, ldv, sigma, &
         bmat, dm, which, nev, tol, resid, ncv, v, ldv, &
         iparam, ipntr, workd, workl, lworkl, info )
    
    
    Egs = DX(1)


    write(*, "(A)") "Computing J matrix" 
    do II = 1,dm
       do JJ = II,dm
          x = bosonic_tp_index(II,JJ,dm) 
          Jtot_MAT(x) =  Jtot_elem(II,JJ,basis)          
       end do
    end do

    print*
    write(*, "(A)") "================================================="
    write(*, "(A)") "  ENERGY       EX ENERGY     J(J+1)      SPIN    " 
    write(*, "(A)") "================================================="
    do AA = 1,10
       sm = 0.d0 
       do II = 1, dm
          amp1 = Z(ii,AA)
          do JJ = II, dm
             amp2 = Z(jj,AA)
             x = bosonic_tp_index(II,JJ,dm)
             sm = sm + amp1 * Jtot_mat(x)* amp2
          end do
       end do

       do II = 1, dm
          amp1 = Z(ii,AA)
          do JJ = 1, II-1
             amp2 = Z(jj,AA)
             x = bosonic_tp_index(JJ,II,dm)
             sm = sm + amp1 * Jtot_mat(x)* amp2
          end do
       end do


       
       spin = (sqrt(nint(sm)*4+1.d0 )-1)/2.d0 
       
       write(*,"(4(f12.4))") DX(AA),DX(AA)-EIMSRG,sm,spin
       write(66,"(3(e25.14))") DX(AA),DX(AA)-EIMSRG,sm
       print*
    end do
    t2 = omp_get_wtime()    
    print*, "TIME: ", t2-t1
    deallocate(HAM)
    t2 = omp_get_Wtime() 
    write(*,"(A,f10.1,A,f10.1)") "Time: ", t2-t1, " Total: ", t2-time_Zero 
    
  end subroutine diagonalize


  
  subroutine first_mat_vec_prod(z0,z1,z2,v,w,basis,N)
    implicit none

    type(block_mat_full),allocatable,dimension(:,:) :: z1    
    type(block_mat),allocatable,dimension(:) :: z2
    integer,dimension(:,:) :: basis
    real(8) :: z0,t5,t6,omp_get_wtime
    integer :: N,II,JJ,XX,Imax,Imin,nthr,omp_get_num_threads,q
    real(8),dimension(N) :: v,w
   
    !$OMP PARALLEL
    nthr=omp_get_num_threads()
    !$OMP END PARALLEL
    w = 0.d0 
    !$OMP PARALLEL DO PRIVATE(XX,II,JJ,q),SHARED(N,z0,z1,z2,basis,HAM,nthr,v,w)    
    do q = 1 , nthr   !!! each thread works on one "q" 
       do II = q,N,nthr  ! matrix is triangle, so make sure one thread doesn't do all the work. 
          do JJ = II,N
             XX = bosonic_tp_index(II,JJ,N)
             HAM(XX) =  mat_elem(II,JJ,basis,z1,z2)
             w(II) = w(II) +  v(JJ) * HAM(XX)
          end do
          w(II) = w(II) + v(II) * z0   ! zero body offset
          XX = bosonic_tp_index(II,II,N)
          HAM(XX) = HAM(XX) + z0 
       end do
    end do
    !$OMP END PARALLEL DO

    ! exploit hermiticity
    !$OMP PARALLEL DO PRIVATE(XX,II,JJ,q),SHARED(N,HAM,nthr,v,w)    
    do q = 1 , nthr  
       do II = q,N,nthr  
          do JJ = 1,II-1
             XX = bosonic_tp_index(JJ,II,N)
             w(II) = w(II) +  v(JJ) * HAM(XX) !* hermitian factor
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine first_mat_vec_prod


  subroutine mat_vec_prod(z0,z1,z2,v,w,basis,N)
    implicit none
    
    type(block_mat_full),allocatable,dimension(:,:) :: z1    
    type(block_mat),allocatable,dimension(:) :: z2
    integer,dimension(:,:) :: basis
    real(8) :: z0
    integer :: N,II,JJ,XX,nthr,q,omp_get_num_threads
    real(8),dimension(N) :: v,w

    !$OMP PARALLEL
    nthr=omp_get_num_threads()
    !$OMP END PARALLEL
    
    w = 0.d0
    !$OMP PARALLEL DO PRIVATE(XX,II,JJ,q),SHARED(N,HAM,nthr,v,w)    
    do q = 1 , nthr   !!! each thread works on one "q" 
       do II = q,N,nthr  
          do JJ = II,N
             XX = bosonic_tp_index(II,JJ,N)
             w(II) = w(II) +  v(JJ) * HAM(XX)
          end do

          do JJ = 1,II-1
             XX = bosonic_tp_index(JJ,II,N)
             w(II) = w(II) +  v(JJ) * HAM(XX)! * hermitian factor
          end do
       end do
    end do
    !$OMP END PARALLEL DO 
  end subroutine mat_vec_prod
  
  real(8) function mat_elem(II,JJ,basis,z1,z2)  
!    <II|H|JJ>   in SD basis 
    implicit none

    
    integer :: II ,JJ,q, rank,out,bra_a,bra_b,ket_c,ket_d,phase,sm_i
    integer :: braout,ketout,q1,sm_j
    type(block_mat_full),allocatable,dimension(:,:) :: z1    
    type(block_mat),allocatable,dimension(:) :: z2
    integer,dimension(:,:) :: basis
    integer,dimension(jbas%Abody) ::  bra,ket
    integer,dimension(jbas%Abody-1) ::  smmed
    integer,dimension(mbas%Ntot) :: brabit,ketbit,diffbit
    real(8) :: sm 
    
    !! compare occupancies    
    !! DIAGONAL
    smmed=0
    if (II == JJ ) then
       bra = basis(II,:)       
       sm = 0.d0 
       do q = 1, jbas%Abody
          sm_i = bra(q) 
          sm = sm + get_me1b(sm_i,sm_i,z1)
       end do

       do q = 1, jbas%Abody
          sm_i = bra(q) 
          do q1 = q+1,jbas%Abody
             sm_j = bra(q1)
             sm = sm + get_me2b(sm_i,sm_j,sm_i,sm_j,z2)
          end do
       end do

       mat_elem = sm

    else

       brabit = 0
       ketbit = 0 
       bra = basis(II,:)
       ket = basis(JJ,:)

       do q = 1, jbas%Abody
          brabit(bra(q)) = 1
          ketbit(ket(q)) = 1
          diffbit = brabit - ketbit
          rank = sum(abs(diffbit))/2          
       end do


       if (rank == 1) then
          !! find the indices of the 1body matrix element 
          out = 0 
          do  q = 1, mbas%Ntot
             if (diffbit(q) ==1 ) then                
                out = out + 1
                bra_a = q
             else if (diffbit(q) == -1)then                             
                out = out+1
                ket_c = q
             end if
             if (out ==2) exit
          end do

          !! find phase of ME
          phase = 1
          do q = 1, jbas%Abody
             if (bra(q) == bra_a)then
                smmed(1:q-1) = bra(1:q-1)
                smmed(q:jbas%Abody-1) = bra(q+1:jbas%Abody) 
                phase = phase *(-1)**(q-1)
                exit
             end if
          end do 
          do q = 1, jbas%Abody
             if (ket(q) == ket_c)then
                phase = phase *(-1)**(q-1)
                exit
             end if
          end do

          sm = 0.d0
          do q = 1, jbas%Abody-1
             sm_i = smmed(q) 
             sm = sm + get_me2b(bra_a,sm_i,ket_c,sm_i,z2)
          end do

          mat_elem = phase * (get_me1b(bra_a,ket_c,z1)+sm)
                    
       else if (rank == 2 ) then
          !! find the indices of the 2body matrix element 
          out = 0 
          braout = 0
          ketout = 0
          phase = 1
          do  q = 1, mbas%Ntot
             
             if (diffbit(q) ==1 ) then                
                out = out + 1
                braout = braout+1                
                if (braout ==1 ) then
                   bra_a = q
                else
                   bra_b = q
                end if
             else if (diffbit(q) == -1)then                             
                out = out+1
                ketout = ketout+1                
                if (ketout ==1 ) then
                   ket_c = q
                else
                   ket_d = q
                end if

             end if
             if (out ==4) exit
          end do


          !! get the phase
          do q = 1, jbas%Abody
             if (ket(q) == ket_c)then
                phase = phase *(-1)**(q-1)
                exit
             end if
          end do

          do q = 1, jbas%Abody
             if (ket(q) == ket_d)then
                phase = phase *(-1)**(q)
                exit
             end if
          end do

          if (ket_c > ket_d) then
             phase = phase*(-1)
          end if
          

          do q = 1, jbas%Abody
             if (bra(q) == bra_a)then
                phase = phase *(-1)**(q-1)
                exit
             end if
          end do

          do q = 1, jbas%Abody
             if (bra(q) == bra_b)then
                phase = phase *(-1)**(q)
                exit
             end if
          end do

          if (bra_a > bra_b) then
             phase = phase*(-1)
          end if
        
          mat_elem=get_me2b(bra_a,bra_b,ket_c,ket_d,z2)*phase


          
       else
          mat_elem = 0.d0
          return
       end if
       
    end if
    
  end function mat_elem

  
  subroutine progress_bar( step  )  
    implicit none

    integer :: step

    if ( step .ne. 0) then 
       ! hold the backspace key down
       write(6,'(A)',advance='no') char(8)//char(8)//char(8)//char(8)//char(8) &
            //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)&
            //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8) &
            //char(8)//char(8)//char(8)//char(8)
       flush 6
    end if

    write(6,'((A19),(I5))',advance='no') 'Lanczos iteration: ',step
    flush 6 

  end subroutine progress_bar

  end module mrci_solve
