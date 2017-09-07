module mrci_interactions
  use mrci_basis
  use gzipmod
  implicit none

  type :: block_mat
     real(8),allocatable,dimension(:) :: X 
  end type block_mat

    type :: block_mat_full
     real(8),allocatable,dimension(:,:) :: XX
  end type block_mat_full


  type(block_mat),allocatable,dimension(:),public :: ME2B,Lambda2b
  type(block_mat_full),allocatable,dimension(:,:),public :: ME1B,Lambda1b
  real(8),public :: ME0B

  
contains

  subroutine read_me1b(z0,z1,intfile)
    implicit none

    type(block_mat_full),allocatable,dimension(:,:) :: z1 
    real(8) :: z0 
    real(8) :: me
    integer :: q,a,totme,buflen,ist,bMax,b,endpos
    integer :: a1,a2,a1len,a2len,aa,menpos,lenme,ii
    integer ::  lj,t,l,twoj,nMax,eMax,Lmax
    character(200) :: intfile,me1bfile 
    character(20) :: mes,memstr
    character(3) :: units
    character(10) :: fm,fma2,fma1,fme
    type(c_ptr) :: buf
    integer(c_int) :: hndle,sz
    character(kind=C_CHAR,len=200) :: buffer
    
     me1bfile = intfile(1:len(trim(intfile))-5)//'1b.gz'
  !  me1bfile = intfile(1:len(trim(intfile))-5)//'1b'
    
     hndle=gzOpen(trim(ME_DIR)//trim(adjustl(me1bfile))//achar(0),"r"//achar(0))
     !open(unit=45,file= trim(ME_DIR)//trim(adjustl(me1bfile)) ) 
     Lmax = mbas%lmax
     eMax = 2*mbas%nmax
     sz = 200

     allocate(z1(2,2*Lmax+1))
     
     do t = 0,1 !neutrons are 0, protons 1
        lj = 0
        do l = 0, Lmax
           do  twoj = abs(2*l - 1) , 2*l+1 , 2
              lj=lj+1             
              nMax = (eMax - l)/2
              allocate(z1(t+1,lj)%XX(nMax+1,nMax+1))
              z1(t+1,lj)%XX=0.d0
           end do
        end do
     end do


    ! read verion line, and then some integer
    buf=gzGets(hndle,buffer,sz) 
!    read(45,*) 
    ! the integer probably has to do with the file size
    buf=gzGets(hndle,buffer,sz) 
!    read(45,*) bMax
      
    read(buffer(1:4),'(I4)',iostat=ist) bMax 
    if (ist .ne. 0 ) then 
       read(buffer(1:3),'(I3)',iostat=ist) bMax 
       if (ist .ne. 0 ) then 
          read(buffer(1:2),'(I2)',iostat=ist) bMax
          if (ist .ne. 0 ) then 
             read(buffer(1:1),'(I1)',iostat=ist) bMax
          end if
       end if
    end if
    
    ! me0b 
    buf=gzGets(hndle,buffer,sz) 
!    read(45,*) z0 
    
    read(buffer(1:7),'(f7.1)')  me
        
    if (buffer(1:1)=='-') then
       lenme = 1+digets(floor(abs(me)))+7
    else
       lenme = digets(floor(abs(me)))+7
    end if

    mes = fmtlen(lenme)

    read(buffer(1:lenme),'(f'//trim(mes)//'.6)')  z0 

    !! the rest of the file is "t lj  a  aa  me"
    do   ii = 1, bMax 

 !      read(45,*) t,lj,a1,a2,me
       
       buf=gzGets(hndle,buffer,sz)
!       z1(t+1,lj+1)%XX(a1+1,a2+1)  = me
       
       read(buffer(1:2),'(I2)') t
       read(buffer(3:5),'(I3)') lj
       read(buffer(6:9),'(I4)') a1
       read(buffer(10:12),'(I3)') a2
       read(buffer(13:24),'(f11.6)') z1(t+1,lj+1)%XX(a1+1,a2+1) 
       
    end do

!    close(45) 
    sz=gzclose(hndle)
  end subroutine read_me1b
    
   
  subroutine read_me2b(z2,intfile)
    implicit none

    type(block_mat),allocatable,dimension(:) :: z2
    integer :: q,a,totme,buflen,ist,bMax,b,endpos,i,j
    integer :: a1,a2,a1len,a2len,aa,menpos,lenme,ii,spot
    real(8) :: mem,me,x
    real(8) :: t1,t2,omp_get_wtime
    character(200) :: intfile
    character(20) :: mes,memstr
    character(3) :: units
    character(10) :: fm,fma2,fma1,fme
    logical :: do_read,get_out,found_a1,found_a2
    type(c_ptr) :: buf
    integer(c_int) :: hndle,sz
    character(kind=C_CHAR,len=200) :: buffer
    
    t1 = omp_get_Wtime()
    write(*,*)
    write(*,"(A)") "======================================="
    write(*,"(A)") "READING OPERATOR FILE"    
    write(*,*)
    
    allocate(z2(tp_basis%bMax))

    totme = 0
    do q = 1, tp_basis%bMax
       a = tp_basis%block(q)%aMax
       allocate(z2(q)%X(a*(a+1)/2))
       z2(q)%X = 0.d0
       totme = totme + a*(a+1)/2 
    end do

    write(mes,'(I20)') totme
    mes = adjustl(mes) 
    write(*,"(A)") "Allocated space for "//trim(mes)//" matrix elements." 

    tot_memory = tot_memory +  totme*8.d0
    mem = totme*8.d0/1024.d0/1024.d0
    units = ' MB'
    if (mem > 1024.d0 ) then
       mem = totme*8.d0/1024.d0/1024.d0/2024.d0
       units=' GB'
    end if
    write(memstr,'(f20.3)') mem       
    memstr = adjustl(memstr)

    
    
    
    write(*,"(A)") "Memory: "//trim(memstr)//units
    call print_total_memory 
    write(*,*)
    intfile = adjustl(intfile)
    write(*,"(A)") "Reading interaction from "//trim(intfile)//"..."


    hndle = gzOpen(trim(ME_DIR)//trim(intfile)//achar(0),"r"//achar(0))
    sz = 200

    buf=gzGets(hndle,buffer,sz)  !! comments in file
    buf=gzGets(hndle,buffer,sz)  !! comments in file
    buf=gzGets(hndle,buffer,sz) 
    !! buffer now contains "bMax=XXXXXX"
    !! we would like to know what the length of that string is a priori
    !! I have a good idea because I've already computed bMax 
    buflen = digets(tp_basis%bMax-1)
    fm=fmtlen(buflen)
    
    read(buffer(6:5+buflen),'(I'//trim(fm)//')') bMax 

    if (bMax+1 .ne. tp_basis%bMax) then ! add 1 because c vs. fortran
       print*, "basis bMax does not match INT bMax"
       print*, "basis: ",tp_Basis%bMax
       print*, "INT: ", bMax+1
       stop
    end if

    bMax = bMax + 1 
    sz = 20 !! lines are shorter now I guess
    
    do q = 1,bMax

       buf=gzGets(hndle,buffer,sz) 

       read(buffer(1:7),'(I7)')  b
       read(buffer(10:15),'(I6)')  a

       if ((b+1) .ne. q) then
          print*, 'somehow we lost a block'
          print*, 'basis:',q
          print*, 'INT: ',b+1 
          stop
       end if

       if ((a+1) .ne. tp_Basis%block(q)%aMax) then
          print*, 'block aMax do not match'
          print*, 'block: ',q 
          print*, 'basis:',tp_Basis%block(q)%aMax
          print*, 'INT: ',a+1 
          stop
       end if
    end do

    sz= 200
    get_out = .false.
    do_read=.true. 
    do q = 1, tp_basis%bMax
       if (.not. get_out) buf=gzGets(hndle,buffer,sz) !! empty line
       buf=gzGets(hndle,buffer,sz) !! comment
       get_out = .false. 
       ii = 1
       do a1 = 1 , tp_basis%block(q)%aMax
          do a2 = a1 , tp_basis%block(q)%aMax

             if (do_read) buf=gzGets(hndle,buffer,sz)           
             do_read = .true.

             !! break down the string
             found_a1 = .false. 
             found_a2 = .false. 
             spot = 1 
             do while (.true. )
                if (found_a1) then
                   !! find a2
                   if ( buffer(spot:spot) == " ") then
                      a2len = spot-3-a1len 
                      found_a2 = .true.
                      exit
                   end if
                else
                   !! find a1 
                   if ( buffer(spot:spot) == " ") then
                      a1len = spot-1
                      found_a1 = .true.
                      spot = spot + 1
                   end if
                end if
                spot = spot + 1
             end do
             fma1 = fmtlen(a1len+2)                   
             fma2 = fmtlen(a2len+2)
             
             ! all of this is because gzgets sucks in fortran
             menpos = a1len+a2len+5
             
             read( buffer(menpos:menpos+9),'(f10.2)',iostat=ist) me
             if (ist .ne. 0) then
                
                ! we've reached the end of the block
                get_out = .true.
                exit
             end if
             
             
             
             lenme = digets(floor(abs(me))) 
             
             if (buffer(menpos:menpos) == '-') then
                fme = fmtlen(lenme+10)
                endpos = menpos+9+lenme
             else
                fme = fmtlen(lenme+9)
                endpos = menpos+8+lenme
             end if

             read(buffer(1:2+a1len),'(I'//trim(fma1)//')',iostat=ist)  a

          
             
             if(ist .ne. 0) then
                !! the correct index was not read
                do_read = .false.
                ii = ii + 1
                cycle
             end if
             
             read(buffer(3+a1len:4+a1len+a2len),&
                  '(I'//trim(fma2)//')',iostat=ist) aa

           
             
             if(ist .ne. 0) then
                !! the correct index was not read
                do_read = .false.
                ii = ii +1 
                cycle                
             end if
             
             
             if ((a+1) .ne. a1 ) then
                do_read = .false. 
                ii = ii + 1
                cycle 
             end if

             if ((aa+1) .ne. a2 ) then
                 do_read = .false. 
                 ii = ii + 1
                 cycle
             end if
             
             !! here, we have presumably matched the q, a1 and a2 indices to
             !! b, a and aa
             
             read(buffer(menpos:endpos) ,'(f'//trim(fme)//'.8)') me


             i = tp_Basis%block(q)%qnums(a1,1)
             j = tp_Basis%block(q)%qnums(a1,2) 

             if (i > j ) then
                ! we actually store the opposite orientation
                me = me * (-1)**(( jbas%jj(i) - jbas%jj(j) +   tp_basis%block(q)%J)/2) 
             end if
             
             i = tp_Basis%block(q)%qnums(a2,1)
             j = tp_Basis%block(q)%qnums(a2,2) 
             
             if (i > j ) then
                ! we actually store the opposite orientation
                me = me * (-1)**(( jbas%jj(i) - jbas%jj(j) +   tp_basis%block(q)%J)/2) 
             end if
             
             z2(q)%X(ii) = me
             ii = ii +1
             
          end do
          if (get_out) exit
       end do
    end do
    sz= gzClose(hndle) 
    t2 = omp_get_Wtime() 
    write(*,"(A,f6.1,A,f10.1)") "Time: ", t2-t1, " Total: ", t2-time_Zero 
  end subroutine read_me2b
    


  real(8) function get_me1b(a,b,z1)
    ! a and b are m-scheme indices
    implicit none

    type(block_mat_full),dimension(:,:) :: z1
    integer :: a,b,lj,t
    
    if (mbas%mm(a) .ne. mbas%mm(b) ) then
       get_me1b = 0.d0
       return
    end if

    if (mbas%ll(a) .ne. mbas%ll(b) ) then
       get_me1b = 0.d0
       return
    end if

    if (mbas%jj(a) .ne. mbas%jj(b) ) then
       get_me1b = 0.d0
       return
    end if

    if (mbas%tz(a) .ne. mbas%tz(b) ) then
       get_me1b = 0.d0
       return       
    end if

    t = (1-mbas%tz(a))/2

    if (mbas%ll(a) == 0 )then
       lj = 1
    else
       lj = 1+mbas%ll(a) + (mbas%jj(a)-1)/2 
    end if

    get_me1b = z1(t+1,lj)%XX(mbas%nn(a)+1,mbas%nn(b)+1 )


  end function get_me1b

    real(8) function get_Jme1b(a,b,z1)
    ! a and b are j-scheme indices
    implicit none

    type(block_mat_full),dimension(:,:) :: z1
    integer :: a,b,lj,t
    

    if (jbas%ll(a) .ne. jbas%ll(b) ) then
       get_Jme1b = 0.d0
       return
    end if

    if (jbas%jj(a) .ne. jbas%jj(b) ) then
       get_Jme1b = 0.d0
       return
    end if

    if (jbas%tz(a) .ne. jbas%tz(b) ) then
       get_Jme1b = 0.d0
       return       
    end if

    t = (1-jbas%tz(a))/2

    if (jbas%ll(a) == 0 )then
       lj = 1
    else
       lj = 1+jbas%ll(a) + (jbas%jj(a)-1)/2 
    end if

    get_Jme1b = z1(t+1,lj)%XX(jbas%nn(a)+1,jbas%nn(b)+1 )


  end function get_Jme1b

  real(8) function get_Jme2b(ay,by,cy,dy,JT,z2)
    !! UNNORMALIZED V^2_{ABCD}  jbas indices
    implicit none

    type(block_mat),dimension(:) :: z2
    integer :: JT,j_min,j_start,j_end
    integer :: a,b,c,d,ja,jb,jc,jd,MT,q
    integer :: ax,bx,cx,dx,ay,by,cy,dy
    integer :: ma,mb,mc,md ,A1,A2,BB,Ntot,Amin,Amax,pre
    real(8) :: me,dcgi
    
    
    pre = 1

    ja = jbas%jj(ay)
    jb = jbas%jj(by)
    jc = jbas%jj(cy)
    jd = jbas%jj(dy) 

    if ((JT> ja+jb).or.(JT > jc+jd)) then 
       get_Jme2b = 0.d0
       return
    end if
       
    if ((JT < abs(ja-jb)).or.(JT < abs(jc-jd))) then
       get_Jme2b = 0.d0
       return
    end if
    
    a = ay
    b = by
    if ( ay > by ) then
       b = ay
       a = by       
       pre = pre*(-1)**((ja-jb+JT)/2) 
    end if
    
    c = cy
    d = dy
    if ( cy > dy ) then
       d = cy
       c = dy
       pre = pre*(-1)**((jc-jd+JT)/2) 
    end if

    if ((a==b).and.(mod(JT/2,2)==1)) then
       get_Jme2b = 0.d0
       return
    end if

    if ((c==d).and.(mod(JT/2,2)==1)) then
       get_Jme2b = 0.d0
       return
    end if

    ja = jbas%jj(a)
    jb = jbas%jj(b)
    jc = jbas%jj(c)
    jd = jbas%jj(d) 

    
    if ((jbas%tz(a)+jbas%tz(b) ).ne.(jbas%tz(c)+jbas%tz(d))) then
       get_Jme2b = 0.d0
       return
    end if

    if (mod(jbas%ll(a)+jbas%ll(b)+jbas%ll(c)+jbas%ll(d),2).ne.0) then
       get_Jme2b = 0.d0
       return
    end if

    me = 0.d0

    q = get_tp_block_index(a,b,JT)
    if (q==0) then
       get_Jme2b=0.d0
       return
    end if
       
    Ntot = tp_Basis%block(q)%aMax
    A1 = TP_index(a,b,JT)       
    A2 = TP_index(c,d,JT)
    
    Amin = min(A1,A2)
    Amax = max(A1,A2) 

    me = z2(q)%X(bosonic_Tp_index(Amin,Amax,Ntot))*pre

    get_Jme2b = me
  end function get_Jme2b

  
    
  real(8) function get_me2b(ay,by,cy,dy,z2)
    !! V_{abcd} m-scheme   (m-scheme indices)
    implicit none

    type(block_mat),dimension(:) :: z2
    integer :: JT,j_min,j_start,j_end
    integer :: a,b,c,d,ja,jb,jc,jd,MT,q
    integer :: ax,bx,cx,dx,ay,by,cy,dy
    integer :: ma,mb,mc,md ,A1,A2,BB,Ntot,Amin,Amax,pre
    real(8) :: me,dcgi

    ax = mbas%jlab(ay)
    bx = mbas%jlab(by)
    cx = mbas%jlab(cy)
    dx = mbas%jlab(dy)

    pre = 1
    a = ay
    b = by 
    if (ax > bx ) then
       a = by
       b = ay
       pre = pre*(-1)
    end if

    c = cy
    d = dy 
    if (cx > dx ) then
       c = dy
       d = cy
       pre = pre*(-1)
    end if

    
    ja = mbas%jj(a)
    jb = mbas%jj(b)
    jc = mbas%jj(c)
    jd = mbas%jj(d) 

    ma = mbas%mm(a)
    mb = mbas%mm(b)
    mc = mbas%mm(c)
    md = mbas%mm(d) 
    
    if ((ma+mb).ne.(mc+md) ) then
       get_me2b = 0.d0
       return
    end if
    
    if ((mbas%tz(a)+mbas%tz(b) ).ne.(mbas%tz(c)+mbas%tz(d))) then
       get_me2b = 0.d0
       return
    end if

    if (mod(mbas%ll(a)+mbas%ll(b)+mbas%ll(c)+mbas%ll(d),2).ne.0) then
       get_me2b = 0.d0
       return
    end if
    
    j_start = max(abs(ja-jb),abs(jc-jd))
    j_end = min(ja+jb,jc+jd)

    MT = ma+mb

    me = 0.d0

    ax = mbas%jlab(a)
    bx = mbas%jlab(b)
    cx = mbas%jlab(c)
    dx = mbas%jlab(d)

    do JT = j_start,j_end,2
       q = get_tp_block_index(ax,bx,JT)
       if (q==0) cycle
       Ntot = tp_Basis%block(q)%aMax
       A1 = TP_index(ax,bx,JT)       
       A2 = TP_index(cx,dx,JT)

       if (A1 == 0) cycle
       if (A2 == 0) cycle 
       Amin = min(A1,A2)
       Amax = max(A1,A2) 

       
       me = me + z2(q)%X(bosonic_Tp_index(Amin,Amax,Ntot)) &
            *dcgi(ja,ma,jb,mb,JT,MT)*dcgi(jc,mc,jd,md,JT,MT)*pre
    end do

    get_me2b = me
  end function get_me2b
    
  subroutine unnormal_order(z0,z1,z2)
    implicit none

    real(8) :: z0
    type(block_mat_full),dimension(:,:) :: z1
    type(block_mat),dimension(:) :: z2
    integer :: t,lj,nMax,ljMax,n1,n2,a,b,i,j,ji
    integer :: J_min,J_max,JT,jj,ja,jb,A1,A2,q
    real(8) :: sm,dir,ex,t1,t2,omp_get_Wtime

    t1 = omp_get_Wtime() 
    print* 
    write(*,"(A)") '=====================================' 
    write(*,"(A)") 'UNNORMAL ORDERING INTERACTION'
    print*
    
    ljMax = size(jbas%tlj_to_ab(1,:))

!!! zero body peice

    do i = 1, jbas%Ntot
       ji = jbas%jj(i)
       do j = 1, jbas%Ntot           
          z0 = z0 - (ji+1.d0) * get_jme1b(i,j,z1) * get_jme1b(j,i,lambda1b)           
       end do
    end do

    do i = 1, jbas%Ntot
       ji = jbas%jj(i)
       do j = 1,jbas%Ntot
          jj = jbas%jj(j)

          do a = 1, jbas%Ntot
             ja = jbas%jj(a)
             do b = 1, jbas%Ntot
                jb = jbas%jj(b)
                
                J_min = max(abs(ji-jj),abs(ja-jb))
                J_max = min(ji+jj,ja+jb)


                dir = get_Jme1b(i,a,lambda1b) * get_Jme1b(j,b,lambda1b) 
                ex =  get_Jme1b(i,b,lambda1b) * get_Jme1b(j,a,lambda1b)

                
                do JT = J_min,J_max,2
                   z0 = z0 + (JT+1.d0) * get_Jme2b(i,j,a,b,JT,z2) * &
                        ( dir  - (get_Jme2b(i,j,a,b,JT,lambda2b) + dir &
                        - (-1)**((ji+jj+JT)/2)*ex)*0.25d0)
                end do

                
             end do
          end do
       end do
    end do
    
    
!!! one body piece
    do t = 1,2
       do lj=1, ljMax

          ! all of these n qnums start at 1 because fortran, when
          ! they actually start at zero. 
          nMax = size(jbas%tlj_to_ab(t,lj)%Z)

          do n1 = 1, nMax
             a = jbas%tlj_to_ab(t,lj)%Z(n1)
             ja = jbas%jj(a) 
             do n2 = 1,nMax
                b = jbas%tlj_to_ab(t,lj)%Z(n2)
                jb = jbas%jj(b) 

                sm = 0.d0 
                do i = 1, jbas%Ntot
                   ji = jbas%jj(i) 
                   do j = 1, jbas%Ntot
                      jj = jbas%jj(j) 

                      J_min = max(abs(ji-ja),abs(jj-jb))
                      J_max = min(ji+ja,jj+jb)

                      do JT = J_min,J_max,2
                         z1(t,lj)%XX(n1,n2) = z1(t,lj)%XX(n1,n2)  - (JT+1.d0)/(ja+1.d0)* &
                              get_Jme2b(i,a,j,b,JT,z2) *get_jme1b(i,j,lambda1b)                   
                      end do
                   end do
                end do

             end do
          end do
       end do
    end do

    write(*,"(A)") "Unnormal ordering successful" 
    write(*,"(A,f12.4)") "Unnormal-ordered offset: ",z0
    t2 = omp_get_Wtime() 
    write(*,"(A,f8.1,A,f10.1)") "Time: ", t2-t1, " Total: ", t2-time_Zero 

  end subroutine unnormal_order

  subroutine traces
    implicit none
  
    integer :: ii,jj,Jtot,A 
    real(8) :: sm

    A = mbas%Abody

    write(*,*)
    
    sm = 0.d0
    do ii = 1, jbas%ntot
       sm = sm + (jbas%jj(ii)+1.d0)*get_jme1b(ii,ii,lambda1b)
    end do

    write(*,"(A,f9.2)") "Tr[rho1b]:",sm
    if ( abs(sm - A ) > 1e-3  ) then
       STOP "one-body density matrix has incorrect trace!"
    end if

    sm=0.d0
    do ii = 1, jbas%ntot
       do jj = 1,jbas%ntot
          
          do Jtot = 0,18,2 
             
             sm = sm + (Jtot+1.d0)*get_Jme2b(ii,jj,ii,jj,Jtot,lambda2b) 
          end do
          
          sm = sm   + get_jme1b(ii,ii,lambda1b)*get_jme1b(jj,jj,lambda1b)&
               *(jbas%jj(ii)+1.d0)*(jbas%jj(jj)+1.d0) &
               - get_jme1b(ii,jj,lambda1b)*get_jme1b(jj,ii,lambda1b)*(jbas%jj(ii)+1.d0) 
       end do
    end do

    write(*,"(A,f9.2)") "Tr[rho2b]:",sm
    if ( abs(sm - A*(A-1) ) > 1e-2  ) then
       STOP "two-body density matrix has incorrect trace!"
    end if

  end subroutine traces
  
    
 end module
