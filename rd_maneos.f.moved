 
c----------------------------------------------------------------------
      implicit none

      integer nt,nw,n,nr, lun, i, j, ir, jt, jst

      parameter (nt=233,nr=339)

      double precision a(NR*NT*3 + NR + NT + 2),rhov(nr),tv(nt),
     *                 rho(nr,nt),t(nr,nt),p(nr,nt),e(nr,nt),h(nr,nt)
c----------------------------------------------------------------------- 



      lun = 112


      open (lun,file='MANEOS_SIO2.txt')

      nw = nr*nt*3 + nr + nt

      read (lun,*) (a(i), i = 1, nw)
     
      n = 2

      do i = 1, nr 

         n = n + 1
         rhov(i) = a(n)*1d3 

      end do 

      do i = 1, nt

         n = n + 1
         tv(i) = a(n) 

      end do 

      do i = 1, nr 
         do j = 1, nt

            rho(i,j) = rhov(i)

            if (rhov(i).ge.2600d0) then 
               ir = i
               exit
            end if 

         end do 
      end do 

      do i = 1, ir 
         do j = 1, nt

            if (tv(j).lt.2001d0) jst = j

            t(i,j) = tv(j)

            if (tv(j).ge.7d3) then 
               jt = j
               exit
            end if 

         end do 
      end do 

      do i = 1, nr 
         do j = 1, nt

            n = n + 1

            p(i,j) = a(n)*1d4

         end do 
      end do 

      do i = 1, nr 
         do j = 1, nt

            n = n + 1

            e(i,j) = a(n)*1d3

         end do 
      end do 

      do i = 1, nr 
         do j = 1, nt

            n = n + 1

            h(i,j) = a(n)*1d3

         end do 
      end do 

      close (lun)

      open (lun,file='MANEOS_SIO2.tab')

      write (lun,'(a)') '|6.6.6'
      write (lun,'(a)') 'output_from_rd_MANEOS'
      write (lun,*) 2

      write (lun,'(a9)') 't(K)       '
      write (lun,'(f14.7)') tv(1)
      write (lun,'(g14.7)') 0.
      write (lun,'(i4)') jt - jst + 1

      write (lun,'(a9)') 'rho_kg/m'
      write (lun,'(f14.7)') rhov(1)
      write (lun,'(g14.7)') 0.
      write (lun,'(i4)') ir

      write (lun,*) 5
      write (lun,'(9a)') 'T(K) rho(kg/m3) p(bar) e(kJ/kg) a(kJ/kg)'

      do i = 1, ir
         do j = jst, jt
            write (lun,'(120(g14.7,1x))') 
     *            t(i,j), rho(i,j), p(i,j), e(i,j), h(i,j)
         end do 
      end do 

      close (lun)



      end 
    
   