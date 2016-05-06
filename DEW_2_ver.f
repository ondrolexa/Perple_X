c----------------------------------------------------------------------
c DEW_2_ver converts aqueous species data from the DEW spreadsheet 
c to perple_x format/units. the program does not write a formula
c for the converted data, therefore such a formula must be added
c before the data is entered into a perple_x 



c rewrite 2010 is a program to rewrite *ver.dat thermodynami data files
c from before May 2010 to the current format.

c to run this code you must temporarily modify perplex_parameters such that:

c k5 = max number of components in the data base (<=k0)
c m7 = number of parameters in a transition (was formerly m7 = 12)

c----------------------------------------------------------------------

      implicit none
 
      include 'perplex_parameters.h'

      integer i, ibeg, iend, len, ier, iscan, iscnlt, jbeg(3), jend(3),
     *        nel

      double precision nums(13)

      character rec*240, name8*8, elem(50)*2, two*2

      external iscan, iscnlt

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(240)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
c     iam = 11
c                                 the DEW data, space delimited
      open (n9, file='DEW_data.txt', status= 'old')
c                                 a case sensitive list of acceptable elements
c                                 species that contain elements not in this list
c                                 will be rejected.
      open (11,file='DEW_elements.txt',status='old')
      nel = 0
      do
         read (11,'(a)',iostat=ier) two
         if (ier.ne.0) exit
         nel = nel + 1
         elem(nel) = two
      end do 

      if (nel.eq.0) then 
         write (*,*) 'gork'
         stop
      end if 
c
      open (10,file='DEWver.dat')
c                                 DEW data consists of a name, formula
c                                 Gf, Hf, S0, v0, cp0, omega, z, a1, a2, a3, a4, c1, c2, comment
      do 

         call redcd3 (n9,len,ier)
         if (ier.ne.0) exit
c                                find end of first name
         ibeg = 1 
         iend = iscan (ibeg,len,' ') - 1

         jbeg(1) = 1
         jend(1) = iend
c                                find the second name
         ibeg = iscnlt (iend+1,len,' ')
         iend = iscan (ibeg,len,' ') - 1 

         jbeg(2) = ibeg
         jend(2) = iend
c                                find a comment, if there is one
         ibeg = iscnlt (iend+1,len,'@') 
         jbeg(3) = ibeg
         jend(3) = len
c                                write the numbers to record
         write (rec,'(240a1)') (chars(i), i = jend(2)+1, len)
c                                read the numbers
         read (rec,*) nums 

         if (chars(jbeg(1)).eq.'"') then 
            jbeg(1) = jbeg(1) + 1
            jend(1) = jend(1) - 1
         end if

         if (jbeg(3).eq.jend(3)) jend(3) = jbeg(3) - 1

         if (jend(1)-jbeg(1).gt.7) then
            write (*,*) 'truncated name from ',
     *                  (chars(i),i=jbeg(1),jend(1)),' to ',
     *                  (chars(i),i=jbeg(1),jbeg(1)+7)
            jend(1) = jbeg(1) + 7
         end if 

         write (name8,'(8a1)') (chars(i),i=jbeg(1),jend(1))
         write (*,'(a,a,f8.0,240a1)') name8,' EoS = 16 | Gf = ',nums(1),
c                                 formula
     *                                (chars(i),i=jbeg(2)-1,jend(2)+1),
c                                 comment
     *                                (chars(i),i=jbeg(3)-1,jend(3))
      end do

      

      end 

      subroutine redcd3 (nloc,len,ier)
c----------------------------------------------------------------------
c readcd - read 240 column card image from unit 9, strip out unwanted
c characters. ier = 1 no card found.
c----------------------------------------------------------------------    
      implicit none

      integer len, ier, iscan, ict, i, iscnlt, ibeg, nloc

      character card*240

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(240)
c----------------------------------------------------------------------

      ier = 0 

      ibeg = 0
  
      len = 0 

      card = ' '

      do while (ibeg.ge.len) 

         read (nloc,'(a)',end=90) card

         if (card.ne.' ') then 

            read (card,'(240a)') chars
c                                 find end of data marker '|'
            len = iscan (1,240,'|') - 1
c                                 '|' in first column
            if (len.eq.0) cycle
c                                 find a non blank character
            ibeg = iscnlt (1,len,' ')

         end if 

      end do

      ict = 1

      do i = 2, len 
c                                 strip out '+' and '*' chars
         if (chars(i).eq.'+'.or.chars(i).eq.'*') chars(i) = ' '
c                                 eliminate blanks after '/' and '-'
c                                 and double blanks
         if (chars(ict).eq.' '.and.chars(i  ).eq.' ') cycle
         ict = ict + 1
         chars(ict) = chars(i)

      end do 

      len = ict

      goto 99

90    ier = 3

99    end