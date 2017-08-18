module mrci_interactions
  use mrci_basis
  use gzipmod
  implicit none

  type :: block_mat
     real(8),allocatable,dimension(:) :: X 
  end type block_mat

  type(block_mat),allocatable,dimension(:),public :: ME2B
contains

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


    hndle = gzOpen(trim(intfile)//achar(0),"r"//achar(0))
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
    
    
    
 end module
