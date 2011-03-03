c convert lars stixrude's phase data files to vertex format.

      implicit double precision (a-z)

      character fname*5, junk*1

      integer ier

      do 

         write (*,*) 'file?'
         read (*,*) fname
         if (fname.eq.'end') exit

         open (unit=10,file=fname,iostat=ier)

         if (ier.ne.0) cycle

         read (10,*) junk
         read (10,*) n
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) a0
         read (10,*) v0
         read (10,*) k0
         read (10,*) k0p
         read (10,*) junk
         read (10,*) theta0
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) gam0
         read (10,*) q0
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) junk
         read (10,*) mu0
         read (10,*) mup
         read (10,*) mut

         n   = -n
         a0  = a0*1d3
         v0  = -v0/1d1
         k0  = k0*1d4
         mu0 = mu0*1d4

         write (99,*) fname
         write (*,1000) a0,n,v0,k0,k0p,theta0,gam0,q0,mut,mu0,mut
         write (99,1000) a0,n,v0,k0,k0p,theta0,gam0,q0,mut,mu0,mup

      end do    

1000  format (g13.7,1x,f4.0,1x,f8.4,1x,f8.0,1x,f6.3,1x,f8.3,1x,
     *        3(f6.3,1x),9('0. '),/,f8.0,1x,f6.3,1x,7('0. '))

      end 