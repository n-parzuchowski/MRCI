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


  type(block_mat),allocatable,dimension(:),public :: ME2B
  type(block_mat_full),allocatable,dimension(:,:),public :: ME1B
  real(8),public :: ME0B

  
contains

  subroutine read_me1b(intfile)
    implicit none

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
    
    hndle=gzOpen(trim(ME_DIR)//trim(adjustl(me1bfile))//achar(0),"r"//achar(0))

    Lmax = mbas%lmax
    eMax = 2*mbas%nmax
    sz = 200

    allocate(me1b(2,2*Lmax+1))
    
    do t = 0,1 !neutrons are 0, protons 1
       lj = 0
       do l = 0, Lmax
          do  twoj = abs(2*l - 1) , 2*l+1 , 2
             lj=lj+1             
             nMax = (eMax - l)/2
             allocate(me1b(t+1,lj)%XX(nMax+1,nMax+1))
             me1b(t+1,lj)%XX=0.d0
          end do
       end do
    end do


    ! read verion line, and then some integer
    buf=gzGets(hndle,buffer,sz) 
    ! the integer probably has to do with the file size
    buf=gzGets(hndle,buffer,sz) 

      
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

    read(buffer(1:7),'(f7.1)')  me
    
    
    
    if (buffer(1:1)=='-') then
       lenme = 1+digets(floor(abs(me)))+7
    else
       lenme = digets(floor(abs(me)))+7
    end if

    mes = fmtlen(lenme)

    read(buffer(1:lenme),'(f'//trim(mes)//'.6)')  ME0B

    !! the rest of the file is "t lj  a  aa  me"
    do   ii = 1, bMax 
       
       buf=gzGets(hndle,buffer,sz)
       
       read(buffer(1:2),'(I2)') t
       read(buffer(3:5),'(I3)') lj
       read(buffer(6:9),'(I4)') a1
       read(buffer(10:12),'(I3)') a2
       read(buffer(13:24),'(f11.6)') me1b(t+1,lj+1)%XX(a1+1,a2+1) 
       
    end do
    
    sz=gzclose(hndle)
  end subroutine read_me1b
    
   
  subroutine read_me2b(intfile,tp_basis)
    implicit none

    type(tpd) :: tp_basis
    integer :: q,a,totme,buflen,ist,bMax,b,endpos
    integer :: a1,a2,a1len,a2len,aa,menpos,lenme,ii
    real(8) :: mem,me
    character(200) :: intfile
    character(20) :: mes,memstr
    character(3) :: units
    character(10) :: fm,fma2,fma1,fme
    type(c_ptr) :: buf
    integer(c_int) :: hndle,sz
    character(kind=C_CHAR,len=200) :: buffer
  
    
    allocate(ME2B(tp_basis%bMax))

    totme = 0
    do q = 1, tp_basis%bMax
       a = tp_basis%block(q)%aMax
       allocate(ME2B(q)%X(a*(a+1)/2))
       totme = totme + a*(a+1)/2 
    end do

    write(mes,'(I20)') totme
    mes = adjustl(mes) 
    write(*,*) "Allocated space for "//trim(mes)//" matrix elements." 

    mem = totme*8.d0/1024.d0/1024.d0
    units = ' MB'
    if (mem > 1024.d0 ) then
       mem = totme*8.d0/1024.d0/1024.d0/2024.d0
       units=' GB'
    end if
    write(memstr,'(f20.3)') mem       
    memstr = adjustl(memstr) 

    
    write(*,*) "Total Memory: "//trim(memstr)//units
    write(*,*)
    intfile = adjustl(intfile)
    write(*,*) "Reading interaction from "//trim(intfile)//"..."


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
     !  print*, q,buffer(1:7),buffer(10:16)
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

    do q = 1, tp_basis%bMax
       buf=gzGets(hndle,buffer,sz) !! empty line
       buf=gzGets(hndle,buffer,sz) !! comment

       ii = 1
       do a1 = 1 , tp_basis%block(q)%aMax
          a1len = digets(a1-1)          
          fma1 = fmtlen(a1len+2)
          do a2 = a1 , tp_basis%block(q)%aMax
             a2len = digets(a2-1) 
             fma2 = fmtlen(a2len+2)

             buf=gzGets(hndle,buffer,sz)

             !!! all of this is because gzgets sucks in fortran
             menpos = a1len+a2len+5
             read( buffer(menpos:menpos+9),'(f10.2)') me

             lenme = digets(floor(abs(me))) 
             
             if (buffer(menpos:menpos) == '-') then
                fme = fmtlen(lenme+10)
                endpos = menpos+9+lenme
             else
                fme = fmtlen(lenme+9)
                endpos = menpos+8+lenme
             end if

             read(buffer(1:2+a1len),'(I'//trim(fma1)//')')  a
             read(buffer(3+a1len:4+a1len+a2len),'(I'//trim(fma2)//')') aa      
             read(buffer(menpos:endpos) ,'(f'//trim(fme)//'.8)') me

             if ((a+1) .ne. a1 ) then
                print*, "first index doesn't match"
                print*, "block ",q
                print*, "basis: ",a1
                print*, "INT:", a+1
                stop
             end if

             if ((aa+1) .ne. a2 ) then
                print*, "second index doesn't match"
                print*, "block ",q
                print*, "basis: ",a2
                print*, "INT:", aa+1
                stop
             end if

             ME2B(q)%X(ii) = me
             ii = ii +1
             
          end do
       end do
    end do
    sz= gzClose(hndle) 
 
        
  end subroutine read_me2b
    


  real(8) function get_me1b(a,b)
    ! a and b are m-scheme indeces
    implicit none

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

    get_me1b = me1b(t+1,lj)%XX(mbas%nn(a)+1,mbas%nn(b)+1 )


  end function get_me1b
    
    
  real(8) function get_me2b(a,b,c,d,tp_Basis)
    implicit none

    type(tpd) :: tp_basis
    integer :: JT,j_min,j_start,j_end
    integer :: a,b,c,d,ja,jb,jc,jd,MT,q
    integer :: ax,bx,cx,dx
    integer :: ma,mb,mc,md ,A1,A2,BB,Ntot,Amin,Amax,pre
    real(8) :: me,dcgi
    
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

    print*, ax,bx,cx,dx
    do JT = j_start,j_end,2
       q = get_tp_block_index(ax,bx,JT)
       if (q==0) cycle
       Ntot = tp_Basis%block(q)%aMax
       A1 = TP_index(ax,bx,JT)       
       A2 = TP_index(cx,dx,JT)

       Amin = min(abs(A1),abs(A2))
       Amax = max(abs(A1),abs(A2)) 
       pre = sign(1,A1)*sign(1,A2)

       
       me = me + ME2B(q)%X(bosonic_Tp_index(Amin,Amax,Ntot)) &
            *dcgi(ja,ma,jb,mb,JT,MT)*dcgi(jc,mc,jd,md,JT,MT)*pre
       print*, JT, ME2B(q)%X(bosonic_Tp_index(Amin,Amax,Ntot)),q
    end do

    get_me2b = me
  end function get_me2b
    

    
    
 end module
