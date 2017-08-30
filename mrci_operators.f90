module mrci_operators
  use mrci_basis
  implicit none

contains

  real(8) function get_Jtot_1b(a,b)
    implicit none

    integer :: a,b,ja,jb

    ja= mbas%jj(a)
    jb= mbas%jj(b)

    if (mbas%mm(a) .ne. mbas%mm(b)) then
       get_Jtot_1b = 0.d0
       return
    end if

    if (mbas%tz(a) .ne. mbas%tz(b)) then
       get_Jtot_1b = 0.d0
       return
    end if

    if (mod(mbas%ll(a) + mbas%ll(b),2)==1) then
       get_Jtot_1b = 0.d0
       return
    end if
    
    if (ja==jb)then
       get_Jtot_1b = float(ja)/2.d0*(float(ja)/2.d0+1.d0)       
    else
       get_Jtot_1b = 0.d0
    end if

  end function get_Jtot_1b

  real(8) function get_Jtot_2b(a,b,c,d)
    implicit none

    integer :: a,b,c,d
    
    get_Jtot_2b =(jplus(a,c) * jminus(b,d) + &
          jminus(a,c) * jplus(b,d) + 2*jz(a,c)*jz(b,d) &
          -(jplus(a,d) * jminus(b,c) + &
            jminus(a,d) * jplus(b,c) + 2*jz(a,d)*jz(b,c)))

    
  end function get_Jtot_2b

  real(8) function jplus(a,b)
    implicit none

    integer :: a,b,ja,jb,ma,mb

    ja = mbas%jj(a)
    ma = mbas%mm(a) 
    jb = mbas%jj(b)
    mb = mbas%mm(b)

    if (mbas%tz(a) .ne. mbas%tz(b)) then
       jplus = 0.d0
       return
    end if

    if (mod(mbas%ll(a) + mbas%ll(b),2)==1) then
       jplus = 0.d0
       return
    end if
    
    if (ja .ne. jb) then
       jplus = 0.d0
       return
    end if

    if (ma .ne.  mb + 2) then
       jplus = 0.d0
       return
    end if

    jplus = sqrt((float(jb + mb)/2 +1.d0) *(float(jb-mb)/2.d0))
  end function jplus

   real(8) function jminus(a,b)
    implicit none

    integer :: a,b,ja,jb,ma,mb

    ja = mbas%jj(a)
    ma = mbas%mm(a) 
    jb = mbas%jj(b)
    mb = mbas%mm(b)

    if (mbas%tz(a) .ne. mbas%tz(b)) then
       jminus = 0.d0
       return
    end if

    if (mod(mbas%ll(a) + mbas%ll(b),2)==1) then
       jminus = 0.d0
       return
    end if
    
    if (ja .ne. jb) then
       jminus = 0.d0
       return
    end if
    
    if (ma .ne.  mb - 2) then
       jminus = 0.d0
       return
    end if

    jminus = sqrt(( float(jb - mb)/2.d0 +1.d0) *(float(jb+mb)/2.d0))
  end function jminus
 

   real(8) function jz(a,b)
    implicit none

    integer :: a,b,ja,jb,ma,mb

    ja = mbas%jj(a)
    ma = mbas%mm(a) 
    jb = mbas%jj(b)
    mb = mbas%mm(b)

    if (mbas%mm(a) .ne. mbas%mm(b)) then
       jz = 0.d0
       return
    end if

    if (mbas%tz(a) .ne. mbas%tz(b)) then
       jz = 0.d0
       return
    end if

    if (mod(mbas%ll(a) + mbas%ll(b),2)==1) then
       jz = 0.d0
       return
    end if
    
    if (ja .ne. jb) then
       jz = 0.d0
       return
    end if

    if (ma .ne.  mb) then
       jz = 0.d0
       return
    end if

    jz = float(mb)/2.d0
  end function jz


    real(8) function Jtot_elem(II,JJ,basis)  
!    <II|H|JJ>   in SD basis 
    implicit none

    
    integer :: II ,JJ,q, rank,out,bra_a,bra_b,ket_c,ket_d,phase,sm_i
    integer :: braout,ketout,q1,sm_j
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
          sm = sm + get_Jtot_1b(sm_i,sm_i)
       end do

       do q = 1, jbas%Abody
          sm_i = bra(q) 
          do q1 = q+1,jbas%Abody
             sm_j = bra(q1)
             sm = sm + get_Jtot_2b(sm_i,sm_j,sm_i,sm_j)
          end do
       end do

       Jtot_elem = sm

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
             sm = sm + get_Jtot_2b(bra_a,sm_i,ket_c,sm_i)
          end do

          Jtot_elem = phase * (get_Jtot_1b(bra_a,ket_c)+sm)
                    
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
        
          Jtot_elem=get_Jtot_2b(bra_a,bra_b,ket_c,ket_d)*phase


          
       else
          Jtot_elem = 0.d0
          return
       end if
       
    end if
    
  end function Jtot_elem

end module
