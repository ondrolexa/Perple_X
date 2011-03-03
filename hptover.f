c   program to read the original thermodynamic datafile of H & P 1990

      program trans

      open (9, file='hp96.dat', status= 'old')
      open (10,file='junk.dat')

      call rdin
      end

      subroutine rdin

c   this subroutine reads part of the H&P data

      character     text(132)*1, name*8, cnum*80, gnum*80, twod*80 
      character     snum*80, vnum*80
      real          rnum, comp(16), rgib, g, reas, s, reav
      real          v, a, b, c, e, newa, k298, l1, l2, l3
      integer       index, igst, igen, isst, isen, ivst, iven
      integer       ihptov(16)
c     hp ind      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 
      data ihptov/4,7,3,9,2,8,6,1,5,13,14,15,12,16,10,11/

c component indices in hp vs ver:

c                    hp          ver
c na2o               4            1
c mgo                7            2
c al2o3              3            3
c sio2               9            4
c k2o                2            5
c cao                8            6
c tio2               6            7
c mno                1            8
c feo                5            9
c nio                15           10
c zro2               16           11
c cl2                13           12
c o2                 10           13
c h2o                11           14
c co2                12           15
c electron           14           16

1000  format(132a1)
2000  format(a8,4(i2),' H= ',g13.7)
2002  format(6(g13.7,1x))
2001  format(16(f5.2,1x))
c                          entropies of the elements at 298, 1bar
          icomp = 15

          sk = 64.68
          sna = 51.30
          sca = 41.63
          sc  = 5.74
          sti = 30.63
          sal = 28.35
          ssi = 18.81
          smg = 32.68
          smn = 32.01
          sfe = 27.28 
c                           values from Robie
          so = 205.20 /2.
          sh = 130.70 /2.
          ctoj = 4.184
          scl = ctoj*53.288/2.
          sni = ctoj*7.14
          szr = ctoj*9.32
       
5       read(9,1000,end=900) text

c   inst and inen specify the positions of initial and terminal 
c   character of the phase name in the array called text
c   the phase name is then written to the variable name

      inst = 0
      inen = 0

      do i = 1, 132
            if (text(i) .ne. ' ') then
                  inst = i
                  goto 10
            end if
      end do

10      do i = inst, 132
            if (text(i) .eq. ' ') then
                  inen = i-1
                  goto 20
            end if
      end do
c                          blank card, assume end of
c                          file.
      goto 900 

20      write(name,1000) (text(i), i = inst, inen)

c      reading indices and real volues of components from the array called text
c      indices are read and then written into the variable index as integers.
c      Then the value for the component corresponding to the actual index 
c      is read and written ("internal write" statement) to the varaible cnum.
c      cnum is written to rnum which is then assigned to the corresponding 
c      element of the array compon.

      do i = 1, 16
            nst = 15 + (i-1)*11
            if (text(nst).eq.'0'.and.text(nst-1).eq.' ') goto 50
           write(twod,1000) text(nst-1), text(nst)
           read(twod,*) index
            write(cnum,1000) (text(j), j = nst+1,nst+7)
            read(cnum,*) rnum
            comp(ihptov(index)) = rnum
      end do


c   finding, reading and writing the values for g

50      do i=nst+1, 132
                if(text(i).ne.' ') then
                        igst = i
                        goto 60
                end if
        end do
60      do i=igst+1, 132
                if(text(i).eq.' ') then
                        igen = i-1
                        goto 70
                end if
        end do

70      write(gnum,1000)(text(j),j=igst, igen)
        read(gnum,*) rgib
        g = rgib

c   finding, reading and writing values for s

        do i = igen+1, 132
                if(text(i).ne.' ') then
                        isst = i
                        goto 80
                end if
        end do    
80      do i = isst+1, 132
                if(text(i).eq.' ') then
                        isen = i-1
                        goto 90
                end if
        end do
90      write(snum,1000)(text(j), j = isst, isen)
        read(snum,*) reas
        s = reas

c   finding, reading and writing the values for v

        do i = isen+1, 132
                if(text(i).ne.' ') then
                        ivst = i
                        goto 100
                end if
        end do    
100     do i = ivst+1, 132
                if(text(i).eq.' ') then
                        iven = i-1
                        goto 110
                end if
        end do
110     write(vnum,1000)(text(j), j = ivst, iven)
        read(vnum,*) reav
        v = reav


c   end of first line in h&p data file

c   as the second line of the data file has a fixed format, it will be read in a standard
c   fortran way

      read(9,*) a, b, c, e, newa, crap1, k298, crap2, kt, lam

      if (crap1.ne.1d1.or.crap2.ne.4d0) then 
         write (*,*) ' mayday mayday ',name
      end if 
      if(lam.ne.0) then
         backspace(9)
         read(9,*) a, b, c, e, newa, crap1, k298, crap2, kt, lam, 
     *             l1, l2, l3
      end if

       ilam = lam * 10

       g = g * 1000.
       s = s * 1000.
       a = a * 1000.
       b = b * 1000.
       c = c * 1000.
       d = d * 1000.
       e = e * 1000.
      tr =298.15d0
       k298 = k298*1000.
       b1 = newa
       b5 = -10.d0*newa
       b6 = k298 * (1 + 1.5d-4 * tr);
       b7 = - k298*1.5d-4
       b8 = 4.d0

       write (10,2000) name,1,0,ilam,0,g
c                                    convert hp H to G
       dsf = s        - comp(1) * sna - comp(2) * smg
     *                - comp(3) * sal - comp(4) * ssi
     *                - comp(5) * sk  - comp(6) * sca
     *                - comp(7) * sti - comp(8) * smn
     *                - comp(9) * sfe - comp(10) * sni
     *                - comp(11) * szr - comp(12) * scl
     *                - comp(13) * so - comp(14) * sh
     *                - comp(15) * sc 

       g = g - 298.15 * dsf 
c                                   convert to oxide stoichiometry
       comp(1) = comp(1)/2.
       comp(3) = comp(3)/2.
       comp(5) = comp(5)/2.
       comp(14) = comp(14)/2.
       comp(13) = comp(13) - comp(1) - comp(2) - 3. * comp(3)
     *             - 2. * comp(4) - comp(5) - comp(6) 
     *             - 2. * comp(7) - comp(8) - comp(9)
     *             - comp(10) - 2. * comp(11)
     *             - comp(14) - 2. * comp(15)

       comp(13) = comp(13) / 2.
       comp(12) = comp(12) / 2.

       write (10,2001) (comp(j),j=1,icomp)
       write (10,2002) g,s,v,a,b,c,0.,e,0.,0.,b1,0.,0.,0.,b5,b6,b7,b8
      if (ilam.ne.0) then
         if (l3.ne.0.d0) then
            write (*,*) name, l2*1000/l3
         else 
            write (*,*) name, ' noclap '
         end if 
         write (10,2002) l1,l2*1000.,l3,0.,0.,0.,0.,0.,0.,0.
      end if 

       do i = 1, icomp
          comp(i) = 0.
       end do 

       goto 5

900     end