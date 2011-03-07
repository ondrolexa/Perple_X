      program specis

      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,jgo,ipl

      character specie(nsp)*3, vname*10, yes*1

      double precision 
     *       x(nsp,1000),f(1000),y(1000),ns(1000),no(1000),nn(1000),
     *       nh(1000),nc(1000),fug(nsp,1000),z1(1000),z2(1000),
     *       x1(1000),x2(1000),x3(1000),nco2, nh2o, nch4

      integer igo, ins(nsp), isp

      double precision tentoe, fo2, fs2, xmn, xmx, xonc, fmn, fmx, 
     *                 tot, fmin, fmax, finc, xfmn, xfmx,
     *                 yh, nht, yo, yh2o, yco2, rno, coht

      double precision p,t,xo,u
      common/ cst5 /p,t,xo,u(6)

      double precision xs,g,v
      common/ cstcoh /xs(nsp),g(nsp),v(nsp)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer iam
      common/ cst4 /iam

      data tentoe, fo2, fs2, specie /2.302585093d0, 0d0, 0d0,
     *      'H2O','CO2','CO ','CH4','H2 ','H2S','O2 ',
     *      'SO2','COS','N2 ','NH3'/
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 12
c                                 version info
      call vrsion

      jgo = 0

50    igo = 0
      ipl = 1
      fs2 = -9999d0*tentoe/2d0
      elag = 0d0
c                                  initialize species fractions/volumes
      do i = 1, nsp
         xs(i) = 0d0
      end do 
c                                  get the users choice of EoS:   
      call rfluid (1, ifug)

      if ((ifug.eq.24.or.ifug.eq.7.or.ifug.eq.8).and.ibuf.ne.3) then
         write (*,*) 'you must choose buffer 3 (user specified fo2) ',
     *               ' for this EoS, try again...'
         goto 50 
      end if 

      write (*,1220) 'Plot (log) species fractions (y/n)?'
      write (*,1220) '(if you answer no, then fugacities are plotted).'
      read (*,1220) yes
      if (yes.eq.'y'.or.yes.eq.'Y') ipl = 0

      igo = 0
      k = 0
c                                  for multispecies fluids set
c                                  up species indices:
      if (ifug.gt.6.and.ifug.lt.13) then
         vname = '  X(O)  '
         if (ifug.eq.7.or.ifug.eq.8) vname = 'log(fo2)'
         isp = 5
         do 1 i = 1, 6
1           ins(i) = i
         if (ifug.eq.12) then
            isp = 8
            ins(7) = 8
            ins(8) = 9
         end if
      else if (ifug.eq.16) then
         vname = '  X(O)  '
         isp = 3
         ins(1) = 1
         ins(2) = 5
         ins(3) = 7
      else if (ifug.eq.17) then
         vname = '  X(O)  '
         isp = 4
         ins(1) = 1
         ins(2) = 5
         ins(3) = 6
         ins(4) = 8
      else if (ifug.eq.13.or.ifug.eq.15) then 
         vname = ' X(H2)  '
      else if (ifug.eq.24) then 
         write (*,*) 'lg f(O2) is log10[f(O2)]'
         vname = 'lg f(O2)'
         isp = 7 
         do i = 1, 5
            ins(i) = i 
         end do 
         ins(6) = 10
         ins(7) = 11
      end if 

      write (*,2900) 
      read (*,*) p, t

      if (ifug.eq.24) then 
         write (*,2940) 
         read (*,*) gz, elag
         elag = elag*tentoe
      end if 

77    write (*,2910) vname
      read (*,*) xmn, xmx, xonc
      if (xmx.lt.xmn) then 
         write (*,2950)
         goto 77
      end if 

78    write (*,2920) 
      read (*,*) fmn, fmx 
      if (fmx.lt.fmn) then 
         write (*,2950)
         goto 78
      end if 
      

2900  format (/,'Enter pressure (bar), temperature (K):')
2910  format (/,'Enter min, max, and increments for the',
     *        ' independent variable (',a,'):')
2920  format ('Enter min-max (log) limits for the dependent variables,'
     *       ,' in the plot file,',/,'values that exceed these limits',
     *        ' will be set to the limits established here:')
2940  format (/,'Enter molar N/C ratio and log10[a(graphite)]:',/)
2950  format (/,'The maximum value must be greater than the minimum',/)

      if (jgo.eq.0) then 

         open (10,file='species.plt')
         open (11,file='GCOH_surf.plt')
         open (12,file='GCOO_surf.plt')
         open (13,file='GCHH_surf.plt')

c                                  this is the max X(O), min X(O),
c                                  max log(Y), min log(Y) for
c                                  psvdraw
         write (10,2030) 0.,' ',' ',' ',' ',xmx,xmn,fmx,fmn,
     *                vname,' log(Y) '
         write (11,2030) 0.,' ',' ',' ',' ',1.,0.,1.,0.,
     *                ' X',' Y '
         write (12,2030) 0.,' ',' ',' ',' ',1.,0.,1.,0.,
     *                ' X ',' Y '
         write (13,2030) 0.,' ',' ',' ',' ',1.,0.,1.,0.,
     *                ' X ',' Y '

      end if 

      write (*,3000) vname,(specie(ins(i)),i=1,isp)
      write (*,3040)

      if (ipl.eq.0) then
         write (*,3030) 
      else 
         write (*,3010) 
      end if 

3000  format (/,'values which follow are ',a,
     *        'and the mole fractions of:',/,1x,10(a3,1x))
3030  format (/,'the plot file includes the log of the '  
     *         ,' C-O-H-S-N atomic proportions and ',
     *        /,'the log of the species fractions.'/)
3010  format (/,'the plot file includes the log of the '  
     *         ,' C-O-H-S-N atomic proportions and ',
     *        /,'the log of the species fugacities.'/)
3040  format ('and log10(fO2).',/)    

      xo = xmn
      if (xo.eq.0d0) xo = 0.00000001d0

      if (ifug.eq.7.or.ifug.eq.8.or.ifug.eq.24) then 
         fmin = xmn*tentoe
         finc = xonc*tentoe
         fmax = xmx*tentoe
         dlnfo2 = fmin
      end if 

c                                  call fluid routine:
10    call cfluid (fo2, fs2)

      if (ifug.eq.24.or.ifug.eq.7.or.ifug.eq.8) then 
c                                  check it didn't hit max fo2
         if (xs(2).eq.1d0) goto 20  

      end if 

      write (*,2000) xo,(xs(ins(i)),i=1,isp),fo2/tentoe

      k = k + 1

      do 90 i = 1, isp
         j = ins(i)
         if (xs(j).gt.1d0.or.xs(j).lt.0d0) xs(j) = 1d0
         fug(j,k) = dlog10(xs(j)*g(j)*p)
         x(j,k) = dlog10(xs(j))

         if (fug(j,k).lt.fmn) then
            fug(j,k) = fmn
         else if (fug(j,k).gt.fmx) then
            fug(j,k) = fmx
         end if 

         if (x(j,k).lt.fmn) then
            x(j,k) = fmn
         else if (x(j,k).gt.fmx) then
            x(j,k) = fmx
         end if 

90    continue

      f(k) = fo2/tentoe
      if (f(k).lt.fmn) f(k) = fmn
      if (f(k).gt.fmx) f(k) = fmx

      if (ifug.ne.7.and.ifug.ne.8.and.ifug.ne.24) then 
         y(k) = xo
      else 
         y(k) = dlnfo2/tentoe
      end if 

      ns(k) = xs(6) + xs(8) + xs(9)
      no(k) = xs(1) + xs(2)*2d0+ xs(3) + xs(7)*2d0 
     *               + xs(8)*2d0 + xs(9) 
      nc(k) = xs(2) + xs(3) + xs(4) + xs(9)
      nh(k) = (xs(1) + xs(5) + xs(6))*2d0 + xs(4)*4d0 + xs(11)*3d0 
      nn(k) = 2d0*xs(10) + xs(11)

      tot = ns(k) + no(k) + nh(k) + nc(k) + nn(k)    

      ns(k) = ns(k) / tot
      no(k) = no(k) / tot
      nc(k) = nc(k) / tot
      nh(k) = nh(k) / tot
      nn(k) = nn(k) / tot
c                                     calculate saturation surface
c                                     position for various coordinate 
c                                     systems:
      coht  = no(k) + nc(k) + nh(k)
      z1(k) = (no(k) + 0.5d0 * nc(k))/coht
      z2(k) = nc(k) * 0.866025d0 / coht
c                                     this is to calculate the
c                                     coordinates in O2-H2O-CO2
       nco2 = nc(k) 
       nh2o = nh(k) / 2d0
       rno = no(k) - 2d0*nc(k) - nh(k)/2d0
       tot = rno + nco2 + nh2o
       yco2 = nco2 / tot
       yh2o = nh2o / tot
       yo = rno / tot
  
       x1(k) = yco2 + 0.5d0 * yo
       x2(k) = -yo * 0.866025d0 
c                                     this is to calculate the
c                                     coordinates in H2-H2O-CH4
       nh2o = no(k) 
       nch4 = nc(k)
       nht = nh(k) - 4d0*nc(k) - 2d0*no(k)
       tot = nht + nch4 + nh2o
       yh2o = nh2o / tot
       yh = nht / tot  
       x3(k) = yh2o + 0.5d0 * yh

      if (ns(k).eq.0d0) ns(k) = 1d0
      if (no(k).eq.0d0) no(k) = 1d0 
      if (nh(k).eq.0d0) nh(k) = 1d0  
      if (nc(k).eq.0d0) nc(k) = 1d0 
      if (nn(k).eq.0d0) nn(k) = 1d0

      xfmx = 1d1**fmx
      xfmn = 1d1**fmn

      if (ns(k).lt.xfmn) ns(k) = xfmn
      if (ns(k).gt.xfmx) ns(k) = xfmx
      if (no(k).lt.xfmn) no(k) = xfmn
      if (no(k).gt.xfmx) no(k) = xfmx
      if (nh(k).lt.xfmn) nh(k) = xfmn
      if (nh(k).gt.xfmx) nh(k) = xfmx
      if (nc(k).lt.xfmn) nc(k) = xfmn
      if (nc(k).gt.xfmx) nc(k) = xfmx
      if (nn(k).lt.xfmn) nn(k) = xfmn
      if (nn(k).gt.xfmx) nn(k) = xfmx

      if (ifug.eq.7.or.ifug.eq.8.or.ifug.eq.24) then 
         dlnfo2 = dlnfo2 + finc
         if (dlnfo2.gt.fmax)  goto 20
      else 
         if (xo.ge.0.9999999d0.and.xo.le.xmx) goto 20

         xo = xo + xonc

         goto (20), igo

         if (xo.gt.xmx) then
            xo = xmx
            if (xo.eq.1.d0) xo = 0.9999999d0
            igo = 1
         end if
      end if  

      goto 10

20    do i = 1, isp
         j = ins(i)
         write (10,2010) 2*k,j,1,1,1,1,1,1,1,0.,specie(j)
         write (10,2010) 0
         if (ipl.eq.0) then 
            write (10,2020) (y(l),x(j,l),l=1,k)
         else 
            write (10,2020) (y(l),fug(j,l),l=1,k)
         end if 
      end do
c                                            O2 fugacity
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'log(fO2)'
         write (10,2010) 0
         write (10,2020) (y(l),f(l),l=1,k)
c                                            bulk fractions:
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'C' 
         write (10,2010) 0
         write (10,2020) (y(l), dlog10(nc(l)),l=1,k)
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'H'
         write (10,2010) 0
         write (10,2020) (y(l),dlog10(nh(l)),l=1,k)
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'O'
         write (10,2010) 0
         write (10,2020) (y(l),dlog10(no(l)),l=1,k)
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'S'
         write (10,2010) 0
         write (10,2020) (y(l),dlog10(ns(l)),l=1,k)
         write (10,2010) 2*k,1,1,1,1,1,1,1,1,0.,'N'
         write (10,2010) 0
         write (10,2020) (y(l),dlog10(nn(l)),l=1,k)

         write (12,2010)  2*k,1,1,1,1,1,1,1,1,0.,'sat co2-h2o-o' 
         write (12,2010) 0
         write (12,2020) (x1(l), x2(l) ,l=1,k)

         write (13,2010) 2*k,1,1,1,1,1,1,1,1,0.,'sat ch4-xh2o-h' 
         write (13,2010) 0
         write (13,2020) (x3(l) , y(l), l=1,k)

         write (11,2010) 2*k,1,1,1,1,1,1,1,1,0.,'sat C-O-H' 
         write (11,2010) 0
         write (11,2020) (z1(l) , z2(l), l=1,k)

2000  format (8(g12.6,2x))
2010  format (i5,1x,8(i3,1x),f4.1,1x,a)
2020  format (6(g12.6,1x))
1210  format (/,'Change EoS, buffer, or graphite activity (y/n)?',/)
1220  format (a)
2030  format ('1',/,'0 0 0',/,'0 0 0 0 0 0 0',/,
     *        g9.1,1x,a162,3(/,a162),/,
     *        '2 1 2 0 0',/,'0 0 0 0. 0. 0. 0. 0.',/,
     *        4(g12.6,1x),/,a,/,a)

      write (*,1210) 
      read (*,1220) yes
      if (yes.eq.'y'.or.yes.eq.'Y') then 
         jgo = 1
         goto 50
      end if 

      end 
