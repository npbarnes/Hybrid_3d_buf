
      MODULE misc

      USE global
!      USE gutsf
!      USE boundary

      contains


!----------------------------------------------------------------------
      real FUNCTION pad_ranf()
! This is the random number generator that works on foo.
!----------------------------------------------------------------------
!      include 'incurv.h'

!      integer function irand
!      integer iflag, irnum
!      external irand

!      irnum = irand(iflag)

!      ranf = irnum/32767.0
!      ranf = irnum/2147483647.

      call random_number(pad_ranf)

      return
      end FUNCTION pad_ranf
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE random_initialize ( )
!----------------------------------------------------------------------
      integer, allocatable :: seed(:)
      integer :: n,istat

      call random_seed(size = n)
      allocate(seed(n))

      open(unit=999,file="/dev/urandom", access="stream",
     x form="unformatted", action="read", status="old", iostat=istat)

      if(istat == 0) then
      read(999) seed
      close(999)
      else
      write(*,*) "Random_Initialize error"
      stop
      endif


      call random_seed(put=seed)

      end SUBROUTINE random_initialize
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      SUBROUTINE debug_random_initialize ( )
!----------------------------------------------------------------------
      integer, allocatable :: seed(:)
      integer :: n,istat
      logical :: ex

      call random_seed(size = n)
      allocate(seed(n))

      inquire(file="./seeds/"//int_to_str(my_rank)//"_seed",
     x exist=ex)


      if(ex) then
      open(unit=999,file="./seeds/"//int_to_str(my_rank)//"_seed",
     x access="stream", form="unformatted", action="read",
     x status="old", iostat=istat)
      read(999) seed
      close(999)
      call random_seed(put=seed)

      else

      call random_initialize()
      call random_seed(get=seed)
      open(unit=998,file="./seeds/"//int_to_str(my_rank)//"_seed",
     x form="unformatted")
      write(998) seed
      close(998)

      endif


      end SUBROUTINE debug_random_initialize
!----------------------------------------------------------------------
      pure function int_to_str(num)
      integer,intent(in) :: num
      character(len=12) :: temp
      character(len=:), allocatable :: int_to_str

      if(num .eq. 0) then

      allocate(character(len=1) :: int_to_str)
      int_to_str = "0"

      else if(num .lt. 0) then

      allocate(character(len=int(log10(float(-num)))+2) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      else

      allocate(character(len=int(log10(float(num)))+1) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      endif
      end function int_to_str

      if(num .eq. 0) then

      allocate(character(len=1) :: int_to_str)
      int_to_str = "0"

      else if(num .lt. 0) then

      allocate(character(len=int(log10(float(-num)))+2) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      else

      allocate(character(len=int(log10(float(num)))+1) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      endif
      end function int_to_str

!----------------------------------------------------------------------
      real FUNCTION ranf()
! This is the random number generator that works on foo.
!----------------------------------------------------------------------
!      include 'incurv.h'

!      integer function irand
!      integer iflag, irnum
!      external irand

!      irnum = irand(iflag)

!      ranf = irnum/32767.0
!      ranf = irnum/2147483647.

      call random_number(ranf)

      return
      end FUNCTION ranf
!----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE Col_Freq(nu,aj)
!c----------------------------------------------------------------------
!      include 'incurv.h'
      
!      real nu(nx,ny,nz),    !collision frequency
!     x     aj(nx,ny,nz,3)    !particle density
!c     x     up(nx,ny,nz,3),  !particle bulk velocity
!c     x     nf(nx,ny,nz),    !fluid density
!c     x     uf(nx,ny,nz,3)   !fluid velocity

!c      parameter(sigma = 1)   !need a reasonable coulomb collsion xsec
!c      real ntot             !total ion density
!c      real ui(3)
!       real maxaj

!c      maxaj = 0.0
!c      do 5 i=1,nx
!c         do 5 j=1,ny
!c            do 5 k=1,nz
!c               do 5 m=1,3
!c                  if (abs(aj(i,j,k,m)) .gt. maxaj) then
!c                     maxaj = abs(aj(i,j,k,m))
!c                     endif
!c 5             continue

!c      write(*,*) 'maxaj....',maxaj
!c      if (maxaj .lt. 1.0) maxaj = 1.0
!c      write(*,*) 'maxaj....',maxaj

!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz    !Kelley p. 464
!               nu(i,j,k) = 65.0 
!c   + 100.0*abs(aj(i,j,k,1))/maxaj +
!c     x                             100.0*abs(aj(i,j,k,2))/maxaj +
!c     x                             100.0*abs(aj(i,j,k,3))/maxaj 
!c               ntot = np(i,j,k) + nf(i,j,k)
!c               do 20 m=1,3
!c                  ui(m) = (np(i,j,k)*up(i,j,k,m)/ntot) + 
!c     x                      (nf(i,j,k)*uf(i,j,k,m)/ntot))
!c 20               continue
!c                  nu(i,j,k) = sqrt(ui(1)**2 + ui(2)**2 + 
!c     x                             ui(3)**2)*sigma
! 10            continue

!c      do 20 j=1,ny
!c         do 20 k=1,nz
!c            nu(2,j,k) = 250.0
!c            nu(3,j,k) = 200.0
!c            nu(4,j,k) = 150.0
!c            nu(nx,j,k) = 250.0
!c            nu(nx-1,j,k) = 200.0
!c            nu(nx-2,j,k) = 150.0
!c 20         continue


!c      do 30 i=1,nx
!c         do 30 k=1,nz
!c            nu(i,2,k) = 250.0
!c            nu(i,3,k) = 2000.0
!c            nu(i,4,k) = 150.0
!c            nu(i,ny,k) = 250.0
!c            nu(i,ny-1,k) = 200.0
!c            nu(i,ny-2,k) = 150.0
!c 30         continue


!c      do 40 i=1,nx
!c         do 40 j=1,ny
!c            nu(i,j,2) = 250.0
!c            nu(i,j,3) = 200.0
!c            nu(i,j,4) = 150.0
!c            nu(i,j,nz) = 250.0
!c            nu(i,j,nz-1) = 200.0
!c            nu(i,j,nz-2) = 150.0
!c 40         continue

!      return
!      end
!c----------------------------------------------------------------------


!----------------------------------------------------------------------
!      SUBROUTINE periodic_bc(i,j,k,ip,im,jp,jm,kp,km)
!----------------------------------------------------------------------
!      include 'incurv.h'

!      if (i .eq. 1) then
!         im = nx
!      else
!         im = i-1
!      endif

!      if (j .eq. 1) then
!         jm = ny
!      else
!         jm = j-1
!      endif

!      if (k .eq. 1) then
!         km = nz
!      else
!         km = k-1
!      endif
 
!      if (i .eq. nx) then
!         ip = 1
!      else
!         ip = i+1
!      endif
 
!      if (j .eq. ny) then
!         jp = 1
!      else
!         jp = j+1
!      endif

!      if (k .eq. nz) then
!         kp = 1
!      else
!         kp = k+1
!      endif

         
!      return
!      end
!----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE check_np(np,nf)
!c----------------------------------------------------------------------
!      include 'incurv.h'
      
!      real np(nx,ny,nz)
!      real nf(nx,ny,nz)

!      real maxnp
!      real nfrac

!      maxnp = 0.0
      
!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz
!c               nfrac = np(i,j,k)/nf(i,j,k)
!c               if (nfrac .gt. 10.0) then
!c                  np(i,j,k) = nf(i,j,k)*10.0
!c                  endif

!               if (np(i,j,k) .gt. maxnp) then
!                  maxnp = np(i,j,k)
!                  endif
! 10            continue

!      write(*,*) 'maxnp, nf....',maxnp,nf(1,1,1)
!      write(*,*) 'np/nf....',maxnp/nf(1,1,1)
!      write(*,*)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE fix_b1z(b1)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real b1(nx,ny,nz,3)
!      real offsetx,offsety,offsetz

!      offsetx = (b1(3,3,3,1) + b1(nx-1,3,3,1) + b1(3,ny-1,3,1) + 
!     x         b1(nx-1,ny-1,3,1) + b1(3,3,nz-1,1) + b1(nx-1,3,nz-1,1) +
!     x         b1(3,ny-1,nz-1,1) + b1(nx-1,ny-1,nz-1,1))/8.0


!      offsety = (b1(3,3,3,2) + b1(nx-1,3,3,2) + b1(3,ny-1,3,2) + 
!     x         b1(nx-1,ny-1,3,2) + b1(3,3,nz-1,2) + b1(nx-1,3,nz-1,2) +
!     x         b1(3,ny-1,nz-1,2) + b1(nx-1,ny-1,nz-1,2))/8.0


!      offsetz = (b1(3,3,3,3) + b1(nx-1,3,3,3) + b1(3,ny-1,3,3) + 
!     x         b1(nx-1,ny-1,3,3) + b1(3,3,nz-1,3) + b1(nx-1,3,nz-1,3) +
!     x         b1(3,ny-1,nz-1,3) + b1(nx-1,ny-1,nz-1,3))/8.0

!      write(*,*) 'b1x offset....',offsetx
!      write(*,*) 'b1y offset....',offsety
!      write(*,*) 'b1z offset....',offsetz
!      write(*,*)

!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz
!               b1(i,j,k,1) = b1(i,j,k,1) - (0.0 + offsetx)
!               b1(i,j,k,2) = b1(i,j,k,2) - (0.0 + offsety)
!               b1(i,j,k,3) = b1(i,j,k,3) - (0.0 + offsetz)
! 10            continue

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_chex_rate(np,nn,up,m,chex_rate)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real nn(nx,ny,nz)
!      real up(nx,ny,nz,3)
!      real vn(3)

!      real cx,cy,cz          !neutral cloud center
!      real rx,ry,rz  
!      real t            !run time
!      real vr           !relative velocity between ions and neutrals
!      real sigma_chex
!      parameter (sigma_chex = 1.0e-24)   !km^2  check this
!      real vol

!      call Neut_Center(m,t,cx,cy,cz)
      
!      chex_rate = 0.0
!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz
!               rx = qx(i) - cx
!               ry = qy(j) - cy
!               rz = qz(k) - cz
!               vn(1) = vsat + rx/t
!               vn(2) = ry/t
!               vn(3) = rz/t
!               vr = sqrt((up(i,j,k,1) - vn(1))**2 + 
!     x                   (up(i,j,k,2) - vn(2))**2 +
!     x                   (up(i,j,k,3) - vn(3))**2)
!               vol = dx*dy*dz_grid(k)
!               chex_rate = chex_rate + 
!     x                     vr*sigma_chex*np(i,j,k)*nn(i,j,k)*vol
! 10            continue 


!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE charge_exchange(np,xp,vp,vp1,m,chex_rate,
!     x                           input_p)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real xp(Ni_max,3)
!      real vp(Ni_max,3)
!      real vp1(Ni_max,3)
!      real vn(3)
!      real input_p(3)

!      real cx,cy,cz          !neutral cloud center
!      real rx,ry,rz,r  
!      real t            !run time
!      real vr           !relative velocity between ions and neutrals
!      real sigma_chex
!      parameter (sigma_chex = 1.0e-24)   !km^2  check this
!      real vol
!      real dNcx
!      integer*4 nchex
!      real rnd
!      real nn           !neutral density
!      real nconst
!      real initial_E

!      nconst = vth*sqrt(pi)

!      call Neut_Center(m,t,cx,cy,cz)
      
!      initial_E = input_E
!      chex_rate = 0.0
!      nchex = 0      
!      do 10 l=1,Ni_tot

!         i=ijkp(l,1)
!         j=ijkp(l,2)
!         k=ijkp(l,3)

!         rx = qx(i) - cx
!         ry = qy(j) - cy
!         rz = qz(k) - cz
!         r = sqrt(rx**2 + ry**2 + rz**2)
!         vn(1) = vsat + rx/t
!         vn(2) = ry/t
!         vn(3) = rz/t
!         vr = sqrt((vp(l,1) - vn(1))**2 + 
!     x             (vp(l,2) - vn(2))**2 +
!     x             (vp(l,3) - vn(3))**2)

!         if (r .gt. 2.33*t) then    !2.33 km/s as distbn limit
                                    !otherwise float underflow
!            nn = 0.0
!         else
!            nn = (No/(4*pi*r*r*t*nconst)) *
!     x                exp(-(r-vo*t)**2 / (vth*t)**2)
!         endif

!            dNcx = dt*vr*sigma_chex*nn
!c            write(*,*) 'dNcx....',dNcx,dNcx/beta
!            rnd = pad_ranf()
!            if (rnd .lt. dNcx) then
!               nchex = nchex + 1
!               vol = dx*dy*dz_grid(k)
!               do 30 m=1,3             !remove neutral energy
!                  vp1(l,m) = vp(l,m)
!                  input_E = input_E -  
!     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
!                  input_p(m) = input_p(m) - mBa*vp(l,m) / beta
! 30            continue
!               vp(l,1) = vn(1)
!               vp(l,2) = vn(2)
!               vp(l,3) = vn(3)
!c               xp(l,1) = xp(l,1)
!c               xp(l,2) = xp(l,2)
!c               xp(l,3) = xp(l,3)
!               do 40 m=1,3             !add ion energy
!                  vp1(l,m) = vp(l,m)
!                  input_E = input_E +  
!     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
!                  input_p(m) = input_p(m) + mBa*vp(l,m) / beta
! 40               continue
!               endif
         
! 10      continue 


!c      write(*,*) 'nchex,chex_rate...',real(nchex)/beta,
!c     x            chex_rate/beta

!      input_chex = input_chex + (input_E - initial_E)
!      chex_rate = (real(nchex))/(dt*beta)
!      write(*,*) 'Normalized charge exchange energy gain...',
!     x            input_chex/(input_E - input_chex - input_bill),
!     x            input_E,input_chex,input_bill 
!      write(*,*) 'Charge exchange rate...',chex_rate
!

!      return
!      end SUBROUTINE charge_exchange
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_gradP(gradP,np,nf,itstep,etemp)
!c----------------------------------------------------------------------
!      include 'incurv.h'
      
!      real gradP(nx,ny,nz,3),
!     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
!     x     etemp(nx,ny,nz)

!      real etemp0
!      parameter (etemp0 = 1.0e5)
!c      parameter (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
!      real a0,ntot

!c      write(*,*) 'tstep...',itstep,dt
!c      etemp = etemp0*exp(-itstep*dt/0.5)

!      do 5 i=1,nx
!         do 5 j=1,ny
!            do 5 k=1,nz
!               etemp(i,j,k) = etemp0*np(i,j,k)/(nf(i,j,k)+np(i,j,k))
! 5             continue

!      do 10 i=2,nx-1
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1
!               ntot = np(i,j,k) + nf(i,j,k)
!               a0 = kboltz*etemp(i,j,k)/(mO*ntot)
!               a1 = kboltz*np(i,j,k)/(mO*ntot)               
!               gradP(i,j,k,1) = a0*( 0.5*(np(i+1,j,k)+np(i,j,k))
!     x                          - 0.5*(np(i,j,k)+np(i-1,j,k)) )/dx 
!     x                  + a1*( 0.5*(etemp(i+1,j,k)+etemp(i,j,k))
!     x                    - 0.5*(etemp(i,j,k)+etemp(i-1,j,k)) )/dx
!               gradP(i,j,k,2) = a0*( 0.5*(np(i,j+1,k)+np(i,j,k))
!     x                          - 0.5*(np(i,j,k)+np(i,j-1,k)) )/dy
!     x                  + a1*( 0.5*(etemp(i,j+1,k)+etemp(i,j,k))
!     x                    - 0.5*(etemp(i,j,k)+etemp(i,j-1,k)) )/dy
!               gradP(i,j,k,3) = a0*( 0.5*(np(i,j,k+1)+np(i,j,k))
!     x                      - 0.5*(np(i,j,k)+np(i,j,k-1)) )/dz_grid(k)
!     x                  + a1*( 0.5*(etemp(i,j,k+1)+etemp(i,j,k))
!     x                - 0.5*(etemp(i,j,k)+etemp(i,j,k-1)) )/dz_grid(k)
! 10         continue

!      return
!      end
!c----------------------------------------------------------------------



!c----------------------------------------------------------------------
!      SUBROUTINE get_np_at_sat(np,m,satnp)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real t
!      real cx,cy,cz
!      real Rsx,Rsy,Rsz                !offset coords of satellite
!      parameter (Rsx = -1.83)         !G1 release (km)
!      parameter (Rsy = -3.42)
!      parameter (Rsz = 1.14)
      

!      call Neut_Center(m,t,cx,cy,cz)

!      ii=nint((cx+Rsx)/dx)
!      jj=rj + nint(Rsy/dy)
!      kk=rk + nint(Rsz/dz_grid(rk)) 
      
!      write(*,*) 'sat coords........',ii,jj,kk
!      write(*,*) 'neutral center....',cx,cy,cz

!      satnp = np(ii,jj,kk)

!      return
!      end
!c----------------------------------------------------------------------

!----------------------------------------------------------------------
      SUBROUTINE get_bndry_Eflux(b1,E)
!----------------------------------------------------------------------
!      include 'incurv.h'
     
      real b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3)

      real vol
      real uf_flux
      real exb_flux

      real mO_q
      mO_q = mO/q


! Energy flux through boundary faces

!c i = 2 face 

      do 20 j=2,ny
         do 20 k=2,nz
            m=1
            i=2
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dy*dz_cell(k)*
     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !+ sign since pos is flux into domain
 20         continue

!c i = nx face

      do 30 j=2,ny
         do 30 k=2,nz
            m=1
            i=nx
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dy*dz_cell(k)*
     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
     x           km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 30         continue
!c**********************
!c j = 2 face

      do 40 i=2,nx
         do 40 k=2,nz
            m=2
            j=2
!            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dx*dz_cell(k)*
     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !+ sign since neg is flux into domain
 40         continue

!c j = ny face

      do 50 i=2,nx
         do 50 k=2,nz
            m=2
            j=ny
!            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dz_cell(k)*
     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 50         continue
!c****************
! k = 2 face

      do 60 i=2,nx
         do 60 j=2,ny
            m=3
!            k=rk-20
            k=2
!            vol = uf(i,j,k,m)*dtsub*dx*dy
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 60         continue

! k = nz face

      do 70 i=2,nx
         do 70 j=2,ny
!            m=3
!            k=rk+20
            k=nz-1
!            vol = uf(i,j,k,m)*dtsub*dx*dy
!            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 70         continue

      return
      end SUBROUTINE get_bndry_Eflux
!----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_ui(uf,nf,up,np,ui)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real uf(nx,ny,nz,3),
!     x     nf(nx,ny,nz),
!     x     up(nx,ny,nz,3),
!     x     np(nx,ny,nz),
!     x     ui(nx,ny,nz,3)     

!      real ntot(3),fnp(3),fnf(3)

!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz

!               ip = i+1
!               jp = j+1
!               kp = k+1

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
!     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
!     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
!     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

!               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
!               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
!               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

!               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
!               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
!               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

!               do 10 m=1,3
!                  ui(i,j,k,m) = fnp(m)*up(i,j,k,m)+fnf(m)*uf(i,j,k,m)
! 10               continue

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_vp_out(vprest,vp_out)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real vprest(3),
!     x     vp_out(3)

!      real theta,phi,mtheta,mphi
!      real v1(3),v2(3)
!      real vdx,vdy,vdz,vdx1,vdy1,vdz1,vm
!      real Ein,delta_E
!      real rnd
!      real anglei

!      phi = 0.0
     
!      if ((vprest(1) .gt. 0) .and. (vprest(2) .gt. 0)) then 
!         phi = atan(vprest(2)/vprest(1))
!c         write(*,*) 'First quadrant....'
!      endif
      
!      if ((vprest(1) .gt. 0) .and. (vprest(2) .le. 0)) then 
!         phi = atan(vprest(2)/vprest(1))
!c         write(*,*) 'Forth quadrant....'
!      endif

!      if ((vprest(1) .le. 0) .and. (vprest(2) .gt. 0)) then  
!         phi = pi - atan(abs(vprest(2))/abs(vprest(1)))
!c         write(*,*) 'Second quadrant....'
!      endif

!      if ((vprest(1) .le. 0) .and. (vprest(2) .le. 0)) then 
!         phi = pi + atan(abs(vprest(2))/abs(vprest(1)))
!c         write(*,*) 'Third quadrant....'
!      endif

!      Ein = vprest(1)**2 + vprest(2)**2 + vprest(3)**2

!      theta = acos(vprest(3)/sqrt(Ein))

!c      write(*,*) 'theta,phi...',theta*180.0/pi, phi*180.0/pi

!      mtheta = -theta
!      mphi = phi

!c      v1(1) = cos(mphi)*vprest(1) + sin(mphi)*vprest(2)
!c      v1(2) = -sin(mphi)*vprest(1) + cos(mphi)*vprest(2)
!c      v1(3) = vprest(3)

!c      v2(1) = cos(mtheta)*v1(1) + sin(mtheta)*v1(3)
!c      v2(2) = v1(2)
!c      v2(3) = -sin(mtheta)*v1(1) + cos(mtheta)*v1(3)

!      rnd = ranf()
!      anglei = asin(rnd)
!      delta_E = Ein*cos(anglei)*cos(anglei)

!      vm = sqrt(delta_E)

!      rnd = ranf()*2*pi
!      vdx = vm*sin(anglei)*cos(rnd)
!      vdy = vm*sin(anglei)*sin(rnd)
!      vdz = vm*cos(anglei)

!c      vdm = sqrt(vdx^2 + vdy^2 + vdz^2)

!      vdx1 = cos(theta)*vdx + sin(theta)*vdz
!      vdy1 = vdy
!      vdz1 = -sin(theta)*vdx + cos(theta)*vdz

!      vp_out(1) = cos(-phi)*vdx1 + sin(-phi)*vdy1
!      vp_out(2) = -sin(-phi)*vdx1 + cos(-phi)*vdy1
!      vp_out(3) = vdz1

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE billiard(vp,vp1,input_p,nn,timestep,bill_rate)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real vp(Ni_max,3),
!     x     vp1(Ni_max,3),
!     x     input_p(3),
!     x     nn(nx,ny,nz)


!      real cx,cy,cz          !neutral cloud center
!      real rx,ry,rz,r  
!      real t            !run time
!      real vr           !relative velocity between ions and neutrals
!      real sigma_Ba
!      real vol
!      real dNcol
!      integer*4 ncol
!      real rnd
!      real nneut           !neutral density
!      real nconst
!      real anglei
!      real vprest(3)
!      real vp_out(3)
!      real initial_E, final_E
!      real vn(3)

!      parameter (sigma_Ba = 3.0e-26)   !km^2  check this
      
!      initial_E = input_E
!      nconst = vth*sqrt(pi)

!      call Neut_Center(cx,cy,cz)
      
!      ncol = 0
!      bill_rate = 0.0      
!      do 10 l=1,Ni_tot

!         i=ijkp(l,1)
!         j=ijkp(l,2)
!         k=ijkp(l,3)

!c         rx = qx(i) - cx
!c         ry = qy(j) - cy
!c         rz = qz(k) - cz
!c         r = sqrt(rx**2 + ry**2 + rz**2)
!c         vn(1) = vsat + rx/t
!c         vn(2) = ry/t
!c         vn(3) = rz/t
!         vn(1) = 0.0
!         vn(2) = 0.0
!         vn(3) = 0.0
!         vr = sqrt((vp(l,1) - vn(1))**2 + 
!     x             (vp(l,2) - vn(2))**2 +
!     x             (vp(l,3) - vn(3))**2)

!c         write(*,*) 'r, 2.33*t....',r,2.33*t
!c         if (r .gt. 2.33*t) then    !2.33 km/s as distbn limit
!c                                    !otherwise float underflow
!c            nneut = 0.0
!c         else
!c            nneut = (5.0*No/(4*pi*r*r*t*nconst)) *
!c     x                exp(-(r-vo*t)**2 / (vth*t)**2)
!c         endif
!         nneut = nn(i,j,k)

!         dNcol = dt*vr*sigma_Ba*nneut
!c     write(*,*) 'dNcol....',dNcol,nneut
!         rnd = ranf()
!         if (rnd .lt. dNcol) then
!            ncol = ncol + 1
!c     vprest(1) = vn(1) - vp(l,1)
!c     vprest(2) = vn(2) - vp(l,2)
!c     vprest(3) = vn(3) - vp(l,3)
!            vprest(1) = -vp(l,1)
!            vprest(2) = -vp(l,2)
!            vprest(3) = -vp(l,3)
            
!            call get_vp_out(vprest,vp_out)
            
!            do 35 m=1,3         !remove ion energy
!c     vp1(l,m) = vp(l,m)
!               input_E = input_E -  
!     x              0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
!               input_p(m) = input_p(m) - mBa*vp(l,m) / beta
! 35         continue
            
!            write(*,*) 'vp in...',vp(l,1),vp(l,2),vp(l,3),
!     x           vp_out(1),vp_out(2),vp_out(3)
!            vp(l,1) = vp_out(1)+vp(l,1)
!            vp(l,2) = vp_out(2)+vp(l,2)
!            vp(l,3) = vp_out(3)+vp(l,3)
!            write(*,*) 'vp out...',vp(l,1),vp(l,2),vp(l,3)
            
!            vol = dx*dy*dz_grid(k)
!            do 40 m=1,3         !add ion energy
!c     vp1(l,m) = vp(l,m)
!               input_E = input_E +  
!     x              0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
!               input_p(m) = input_p(m) + mBa*vp(l,m) / beta
! 40         continue
!         endif
         
! 10   continue 
      
!      final_E = input_E - initial_E
!c     input_bill = input_bill + (input_E - initial_E)
!      bill_rate = (real(ncol))/(dt*beta)
!      write(*,*) 'Normalized collisional energy gain...',
!     x     input_bill/(input_E - input_bill - input_chex) 
!      write(*,*) 'Collision rate...',bill_rate
      
!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_etar(np,aj)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real aj(nx,ny,nz,3)

!      real ajperp

!      do 10 i = 1,nx 
!         do 10 j = 1,ny
!            do 10 k = 1,nz
!               do 10 m = 1,3
!                  etar(i,j,k,m) = 0.0
! 10            continue

!      do 20 i = 1,nx 
!         do 20 j = 1,ny
!            do 20 k = 1,nz
!               ajperp = sqrt(aj(i,j,k,1)**2 + aj(i,j,k,2)**2)
!               if ((np(i,j,k) .gt. 0.0) .and. (ajperp .gt. 1.0)) then
!                  etar(i,j,k,1) = eta_init
!                  etar(i,j,k,2) = eta_init
!                  etar(i,j,k,3) = 0
!               endif
! 20            continue

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE update_nu(nu,np,nf)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real nu(nx,ny,nz)
!      real np(nx,ny,nz)
!      real nf(nx,ny,nz)

!      real ntot
!      real sigma_coulb
!      parameter (sigma_coulb = 1.4e-19)   !units of km^2
!      real vthermal
!      parameter (vthermal = 210.0)  !km/s 1000 K
!      real numax,ntotmax

!      numax = 0.0
!c      ntotmax = 0.0
!      do 10 i=1,nx
!         do 10 j=1,ny
!            do 10 k=1,nz
!               ntot = np(i,j,k)+nf(i,j,k)
!               nu(i,j,k) = (melec/mO)*ntot*sigma_coulb*vthermal
!               if (nu(i,j,k) .gt. numax) then 
!                  numax = nu(i,j,k)
!                  endif
!c               if (ntotmax .gt. ntot) then
!c                  ntotmax = ntot
!c                  endif
! 10            continue

!      write(*,*) 'Nu max.........',numax
!c      write(*,*) 'np+nf max.....',ntotmax

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_nuin(nuin,nn,uf)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real nuin(nx,ny,nz)
!      real nn(nx,ny,nz)
!      real uf(nx,ny,nz,3)

!      real sigma_in
!      parameter (sigma_in = 3.0e-26)

!      do i = 1,nx
!         do j = 1,ny
!            do k = 1,nz
!               nuin(i,j,k) = sqrt(uf(i,j,k,1)**2 + uf(i,j,k,2)**2 + 
!     x                       uf(i,j,k,3)**2)*sigma_in*nn(i,j,k)
!            enddo
!         enddo
!      enddo

!      return
!      end
!c----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_beta()
!----------------------------------------------------------------------
!      include 'incurv.h'

!      real Np_tot
      real src(200)
      real r_xyz(200)
      real tau
      real Nofr(200),Np_tot

      tau = tau_photo
      
      src(2) = Qo !1.5e27   !mol/s at Pluto

      do i = 1,200 
         r_xyz(i) = i*dx 
      enddo

!      Np_tot = 0.0
!      do i = 1,S_radius
!c         Nofr(i) =  (N_o - src(i)*tau)*exp(-r_xyz(i)/(tau*vrad)) + 
!c     x                src(i)*tau         
!         Nofr(i) = Qo*dx/vrad
!         Np_tot = Np_tot + Nofr(i)*dt/tau
!c         write(*,*) 'S...',src(i)
!c         write(*,*) Nofr(i),Np_tot
!      enddo
!c      write(*,*) sum(Nofr),Np_tot

! divide particles up between procnum processors      
!      beta = (Ni_tot_sys/((nx*dx*ny*dy*nz*delz))/nf_init
      beta = (Ni_tot_sys/((qx(nx)-qx(1))*(qy(ny-1)-qy(1))*
     x                    (qz(nz-1)-qz(1))))/nf_init

!      beta = (5e6/(nx*dx*ny*dy*nz*delz))/nf_init

!      write(*,*) 'beta...',beta,Ni_tot,nx*dx*ny*dy*nz*delz,nf_init

!      dNi = (beta*Np_tot)/procnum
!      dNi = (beta*Np_tot)   ! use this for all pickup on one processor
      dNi_sw = nf_init*vsw*((ny-2)*dy*(nz-2)*delz)*(dt/2)*beta
      print *,'dNi_sw....',dNi_sw

!      write(*,*) 'dNi...',dNi
!      if (dNi .lt. 1.0) then 
!         write(*,*) 'dNi lt 1.0....'
!         stop
!      endif

!      write(*,*) 'dNi...',dNi,Np_tot
!      dNi = Ni_max
!      beta = (dNi/Np_tot)

      write(*,*) 'beta, dNi....',beta,dNi
!      stop

      return
      end SUBROUTINE get_beta
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_dipole(bdp)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real bdp(nx,ny,nz,3)
!      real Bo,phi,lam,br,blam,qdpx,qdpy,qdpz
!      real x,y,z,r,M
!      real eoverm
!      parameter(eoverm = q/mO)

!      Bo = 300e-9
!      M = Bo*(5.0)**3
      
!      x0 = dx/2
!      y0 = dy/2
!      z0 = dz_grid(nz/2)/2
      
!      qdpx = qx(ri) + x0
!      qdpy = qy(rj) + y0
!      qdpz = qz(rk) + z0
      
!      do 10 i = 1,nx 
!         do 10 j = 1,ny
!            do 10 k = 1,nz
!               x = qx(i) - qdpx
!               y = qy(j) - qdpy
!               z = qz(k) - qdpz
!               r = sqrt(x**2 + y**2 + z**2)
!               phi = atan(y/x)
!               if ((x .ge. 0) .and. (y .ge. 0)) phi = atan(y/x)
!               if ((x .le. 0) .and. (y .ge. 0)) then 
!                  phi = pi/2 + atan(abs(x/y))
!               endif
!               if ((x .le. 0) .and. (y .le. 0)) then 
!                  phi = pi + atan(abs(y/x))
!               endif
!               if ((x .ge. 0) .and. (y .le. 0)) then 
!                  phi = (3.0*pi/2.0) + atan(abs(x/y))
!               endif

!c               write(*,*) 'phi...',phi,qx(i),qdpx
               
!               lam = asin(z/r)
               
!               br = -M*sin(lam)/r**3
!               blam = M*0.5*cos(lam)/r**3
!               bdp(i,j,k,3) = br*sin(lam) + blam*cos(lam)
!               bx = br*cos(lam) + blam*sin(lam)
!               bdp(i,j,k,1) = bx*cos(phi)
!               bdp(i,j,k,2) = bx*sin(phi)
!               bdp(i,j,k,3) = bdp(i,j,k,3)
! 10         continue
!            bdp(:,:,:,1) = bdp(:,:,:,1)*eoverm
!            bdp(:,:,:,2) = bdp(:,:,:,2)*eoverm
!            bdp(:,:,:,3) = bdp(:,:,:,3)*eoverm

!            open(10,file='bdp.dat',status='unknown',
!     x           form='unformatted')
!            write(10) bdp
!            close(10)
!            return
!            end
!c----------------------------------------------------------------------


      end MODULE misc
