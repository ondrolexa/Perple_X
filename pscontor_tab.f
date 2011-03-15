      PROGRAM PSPLOT 

      implicit none

      include 'perplex_parameters.h'

      integer ier, i

      logical ratio
 
      character y*1

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer ix,iy,mvar
      double precision z,zt 
      common/ dim   /z(nx,ny),zt(nx,ny),ix,iy,mvar

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer  iop0 
      common / basic /iop0
c----------------------------------------------------------------------
c                                 version info
      call vrsion

      ratio = .false.
c                         
      do 
c                                 get input file 
         write (*,1000) 
          
         read (*,'(a)') tfname
c                                 extract the root
         call getrt
         
         open (n4,iostat=ier,file=tfname,status='old')

         if (ier.ne.0) then
       
            write (*,1010) tfname
            read (*,'(a)') y

            if (y.eq.'Y'.or.y.eq.'y') cycle 

            stop

         end if

         exit  

      end do 

      if (jvar.eq.2) then
c                                 for 2d query for ratio plots
         write (*,1020) tfname
         read (*,'(a)') y

         if (y.eq.'Y'.or.y.eq.'y') then

            ratio = .true.
c                                 ratio plot 
            do 
c                                 get input file 
               if (i.eq.1) then 
                  write (*,1040) 'numerator'
               else
                  write (*,1040) 'denominator'         
               end if 
          
               read (*,'(a)') tfname

               if (i.eq.1) then 

                  open (n4,iostat=ier,file=tfname,status='old')

                  call getrt

               else

                  open (n5,iostat=ier,file=tfname,status='old')

               end if 

               if (ier.ne.0) then
       
                  write (*,1010) tfname
                  read (*,'(a)') y

                  if (y.eq.'Y'.or.y.eq.'y') cycle 

                  stop

               end if

               exit  

            end do  

         end if

      end if 
c                                 read plot option file, set
c                                 default transformation
      call rdopt 
c                                 open output file 
      call psopen
c                                 allow drafting options prompt
      iop0 = 0

      write (*,1030) 
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') iop0 = 1

      if (jvar.eq.2) then 
c                                 contour plotting
         call psxypl (ratio)

      else 
c                                 x-y plotting
         call psplt1

      end if 
 
      call psclos
 
      close (n4)
 
1000  format (/,'Enter the complete plot file name [e.g., ',
     *       'my_project.tab or my_project.ctr]:')
1010  format (/,'**warning ver191** cannot find file:',/,a,/,
     *       'run WERAMI/FRENDLY to generate the ',
     *       'file or try a different name (y/n)?')
1020  format (/,'Contour the ratio of values in two contour plot ',
     *       'files (y/n)?',/,'If you answer yes the data from the',
     *       'file just read will define',/,'the numerator of the '
     *       'ratio and you will be prompted next for a file',/,
     *       'containing the data for the denominator.')
1030  format (/,'Modify the default plot (y/n)?')
1040  format (/,'Enter the full name of the plot file name that ',
     *          'contains the ',a,' data',/,
     *          '[e.g., my_project1.tab or my_project1.ctr]:')

      end

c---------------------------------------------------------------------
      subroutine psxypl (ratio)
 
c psxypl - subroutine to output x-y plot.

      implicit none

      include 'perplex_parameters.h'

      logical ratio

      character y*1

      integer i,j,jmn,imn,imx,jop0,ncon,jmx,iop1,jy,jx

      double precision dx,dy,xpmn,xpmx,cmin,cmax,dcon,ypmx,ypmn,
     *                 z0min,z0max

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer  iop0 
      common / basic /iop0

      integer ix,iy,mvar
      double precision z,zt 
      common/ dim   /z(nx,ny),zt(nx,ny),ix,iy,mvar

      double precision zmin,zmax
      common/ stuff /zmax,zmin

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------------
      call redtab (n4)

      if (ratio) then 

         jx = ix
         jy = jy 
         do i = 1, ix
            do j = 1, iy
               zt(i,j) = z(i,j)
            end do
         end do 

         call redtab (n5)

         if (jx.ne.ix.or.jy.ne.iy) then 
            write (*,'(a)') 'the plots do not have consistent ',
     *                      'dimensions'
            stop
         end if 

         do i = 1, ix
            do j = 1, iy
               z(i,j) = zt(i,j)/z(i,j)
            end do
         end do

      end if 

      if (iop0.eq.1) then 

         write (*,1050) 
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then 
            do j = 1, iy
               do i = 1, ix
                  if (z(i,j).ne.0d0) z(i,j) = dlog10(dabs(z(i,j)))
               end do 
            end do  
         end if 

      end if 

      write (*,1060) 
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then 

         write (*,1070) vmx(2),vmn(2),vmx(1),vmn(1)
         read (*,*) ypmx,ypmn,xpmx,xpmn
           
         imn = int(xpmn/dvr(1)) + 1
         imx = int(xpmx/dvr(1)) + 1
         jmn = int(ypmn/dvr(2)) + 1
         jmx = int(ypmx/dvr(2)) + 1

         ix = (imx-imn+1)
         iy = (jmx-jmn+1)
         vmx(1) = xpmn + (ix-1)*dvr(1)
         vmx(2) = ypmn + (iy-1)*dvr(2)
         vmn(1) = xpmn
         vmn(2) = ypmn 
c                                      reload mini matrix:
         do i = 1, ix
            do j = 1, iy
               z(i,j) = z(i+imn-1,j+jmn-1)
            end do
         end do

      end if 
c                                 get some options and
c                                 set up transformations
      call psaxop (1,jop0,iop1)
        
      zmin = 1d9
      zmax = -1d9
      z0min = 1d30
      z0max = -1d30
c                                      set up contour intervals                                      
      do i = 1, ix
         do j = 1, iy 
            if (z(i,j).lt.zmin) zmin = z(i,j)
            if (z(i,j).gt.zmax) zmax = z(i,j)
            if (z(i,j).lt.z0min.and.z(i,j).ne.0d0) z0min = z(i,j)
            if (z(i,j).gt.z0max.and.z(i,j).ne.0d0) z0max = z(i,j)
         end do 
      end do 
c                                      set up contour intervals
      write (*,1020) zmin, zmax, z0min, z0max 
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then

         write (*,1030) 
         read (*,*) cmin, cmax, dcon
         ncon = int((cmax-cmin)/dcon) + 1

      else 

         dcon = (zmax-zmin)/11.
         cmax = zmax - 0.5d0 * dcon
         cmin = zmin + 0.5d0 * dcon
         ncon = 11

      end if 

      call pscontor (cmin,ncon,dcon)
 
      call psaxes (jop0)
 
1020  format ('Contoured variable range:',g14.6,'->',g14.6,/,
     *        'Range excluding zero values:',g14.6,'->',g14.6,/,
     *        'Modify default contour interval (y/n)?')
1030  format ('Enter min, max and interval for contours:')
1050  format ('Contour log10 of the z-value (y/n)?')
1060  format (/,'Reset plot limits (y/n)?')
1070  format (/,'Old values were: ',4(g12.4),/,'Enter new values:')

      end

