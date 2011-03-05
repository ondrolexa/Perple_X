      PROGRAM PSCNTR 

      implicit none

      include 'perplex_parameters.h'

      integer ier, i

      logical ratio
 
      character yes*1

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer  iop0 
      common / basic /iop0
c----------------------------------------------------------------------
c                                 version info
      call vrsion

      write (*,1020) 
      read (*,'(a)') yes

      if (yes.ne.'Y'.and.yes.ne.'y') then

         ratio = .false.
c                                 simple plot 
         do 
c                                 get input file 
            write (*,1000) 
          
            call readrt

            call mertxt (tfname,prject,'.ctr',0)
         
            open (n4,iostat=ier,file=tfname,status='old')

            if (ier.ne.0) then
       
               write (*,1010) tfname
               read (*,'(a)') yes

               if (yes.eq.'Y'.or.yes.eq.'y') cycle 

               stop

            end if

            exit  

         end do  

      else 

         ratio = .true.

         do i = 1, 2 
c                                 ratio plot 
            do 
c                                 get input file 
               if (i.eq.1) then 
                  write (*,1040) 
               else
                  write (*,1050)          
               end if 
          
               call readrt

               call mertxt (tfname,prject,'.ctr',0)
         
               if (i.eq.1) then 
                  open (n4,iostat=ier,file=tfname,status='old')
               else
                  open (n4,iostat=ier,file=tfname,status='old')          
               end if 
               if (ier.ne.0) then
       
                  write (*,1010) tfname
                  read (*,'(a)') yes

                  if (yes.eq.'Y'.or.yes.eq.'y') cycle 

                  stop

               end if

               exit  

            end do  

         end do 

      end if 
c                                 read plot option file, set
c                                 default transformation
      call rdopt 
c                                 open output file 
      call psopen
c                                 allow drafting options prompt
      iop0 = 0
      write (*,1030) 
      read (*,'(a)') yes
      if (yes.eq.'y'.or.yes.eq.'Y') iop0 = 1

      call psxypl (ratio)
 
      call psclos
 
      close (n4)
 
1000  format (/,'Enter the CONTOUR plot file name [',
     *       'without the .ctr suffix]:')
1010  format (/,'**warning ver191** cannot find file:',/,a,/,
     *       'run WERAMI to generate the ',
     *       'file or try a different name (y/n)?')
1020  format (/,'Contour the ratio of the values in two contour plot',
     *       'files (y/n)?')
1030  format (/,'Modify the default plot (y/n)?')
1040  format (/,'Enter the numerator CONTOUR plot file name [',
     *       'without the .ctr suffix]:')
1050  format (/,'Enter the denominator CONTOUR plot file name [',
     *       'without the .ctr suffix]:')

      end

c---------------------------------------------------------------------
      subroutine psxypl (ratio)
 
c psxypl - subroutine to output x-y plot.

      implicit none

      include 'perplex_parameters.h'

      logical ratio

      character y*1, fname*162

      integer nx,ny,i,j,iox,ioy,jmn,imn,imx,jop0,ncon,jmx,iop1,jy,jx

      double precision dx,dy,xpmn,xpmx,cmin,cmax,dcon,ypmx,ypmn,
     *                 z0min,z0max
 
      parameter (nx=500,ny=500)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer  iop0 
      common / basic /iop0

      integer ix,iy
      double precision z,zt 
      common/ dim   /z(nx,ny),zt(nx,ny),ix,iy

      double precision zmin,zmax
      common/ stuff /zmax,zmin

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------------
      read (n4,'(a)') fname
      read (n4,*) ix,iy,xmin,ymin,dx,dy
      read (n4,'(10a)') (vnm(i),i=1,2)
      if (ix.gt.nx) call error (1,dx,nx,'NX, PSXYPL')
      if (iy.gt.ny) call error (1,dx,ny,'NY, PSXYPL')

      if (ratio) then 
         read (n4,'(10a)') (vnm(i),i=1,2)
         read (n5,'(a)') fname
         read (n5,*) jx,jy,xmin,ymin,dx,dy
         read (n5,'(10a)') (vnm(i),i=1,2)
         if (jx.ne.ix.or.jy.ne.iy) then 
            write (*,'(a)') 'the contour plots do not have the same ',
     *                      'dimensions'
            stop
         end if 
      end if 

      read (n4,*) ((z(i,j), i = 1, ix), j = 1, iy)

      if (ratio) then 
         read (n5,*) ((zt(i,j), i = 1, ix), j = 1, iy)
         do i = 1, nx
            do j = 1, ny
               z(i,j) = z(i,j)/zt(i,j)
            end do
         end do
      end if 

      if (iop0.eq.1) then 
         write (*,1050) 
         read (*,'(a)') y
         if (y.eq.'y') then 
            do j = 1, ny
               do i = 1, nx
                  if (z(i,j).ne.0d0) z(i,j) = dlog10(dabs(z(i,j)))
               end do 
            end do  
         end if 
      end if 

      ypmn = ymin
      xpmn = xmin
      ypmx = ymin + (iy-1)*dy
      xpmx = xmin + (ix-1)*dx    
      ymax = ypmx
      xmax = xpmx     

      write (*,1060) 
      read (*,'(a)') y
      if (y.eq.'y') then 
         write (*,1070) ypmx,ypmn,xpmx,xpmn
         read (*,*) ypmx,ypmn,xpmx,xpmn
      end if 

      iox = ix
      ioy = iy

      if (y.eq.'y') then 
         jmn = int(ypmn/dy) + 1
         jmx = int(ypmx/dy) + 1
         imn = int(xpmn/dx) + 1
         imx = int(xpmx/dx) + 1
         iy = (jmx-jmn+1)
         ix = (imx-imn+1)
         ymax = ypmn + (iy-1)*dy
         ymin = ypmn 
         xmin = xpmn
         xmax = xpmn + (ix-1)*dx
      end if 
c                                      reload mini matrix:
      if (y.eq.'y') then
         do i = 1, ix
            do j = 1, iy
               z(i,j) = z(i+imn-1,j+jmn-1)
            end do
         end do
      end if 

      vmn(1) = xmin
      vmn(2) = ymin
      vmx(1) = xmax
      vmx(2) = ymax
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

      if (y.eq.'y') then
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

