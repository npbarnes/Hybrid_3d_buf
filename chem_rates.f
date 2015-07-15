
      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp_dd
      USE grid_interp
      
      contains

!----------------------------------------------------------------------
      real FUNCTION neutral_density(r)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real r
      real nn0

!      nn0 = nf_init*2e6
!      nn0 = Ncol*(pwl+1)/RIo
      nn0 = 2e8*1e15
! gaussian 
!      neutral_density = nn0*exp(-(r)**2/(RIo)**2)                  

! power law profile
!      if (r .le. RIo) then
!         neutral_density = nn0*10
!      endif
!      if (r .gt. RIo) then 
!c         neutral_density =  nn0*(RIo/r)**(pwl)
!         neutral_density =  nn0*exp(-(r-RIo)/20.0)
!      endif

! Pluto isotropic escape

      neutral_density = Qo/(4*PI*r**2*vrad)

! Pluto Strobel Atm

!      neutral_density = 4e21*(RIo/r)**16 + 4e16*(RIo/r)**5.0 !+ 
!c     x     3.4e27/(4*PI*(r*1e3)**2*100.)               !m^-3
!      neutral_density = neutral_density*1e9 !km^-3


!      if (neutral_density .ge. 1e22) then 
!         neutral_density = 1e22
!      endif

!      write(*,*) 'nden...',neutral_density

      return
      end FUNCTION neutral_density
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE res_chex(xp,vp,vp1)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real cx,cy,cz,r,vrel
      real nn0,nn !,neutral_density
      real chex_tau,chex_prob

      real sigma_chex
      PARAMETER (sigma_chex = 1e-25)  !10 A^2 in units of km^2

      nn0 = Ncol*(pwl+1)/RIo

      call Neut_Center(cx,cy,cz)
      
      do l = 1,Ni_tot 
         r = sqrt((xp(l,1)-cx)**2 + (xp(l,2)-cy)**2 + 
     x            (gz(ijkp(l,3))-cz)**2) !use global coords
         vrel = sqrt(vp(l,1)**2 + vp(l,2)**2 + vp(l,3)**2)
!         if (r .ge. RIo) then 
         nn = neutral_density(r)
!            nn = nn0*(RIo/r)**(pwl)
!         else
!            nn = nn0
!         endif
         chex_tau = 1./(nn*sigma_chex*vrel)
         chex_prob = dt/chex_tau

         if (pad_ranf() .lt. chex_prob) then
            do m=1,3
               input_E = input_E - 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
!               input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
            enddo                     
            
            vp(l,:) = 0.0
            vp1(l,:) = 0.0
!            m_arr(l) = m_pu*mproton
            mrat(l) = 1./m_pu
!            beta_p(l) = beta_pu
            beta_p(l) = 1.0
!            write(*,*) 'chex...',l,chex_prob
         endif


      enddo

      return
      end SUBROUTINE res_chex
!----------------------------------------------------------------------



!----------------------------------------------------------------------
      SUBROUTINE Ionize_Io(np,vp,vp1,xp,up,ndot)
! Ionizes the neutral cloud with a 28 s time constant and fill particle
! arrays, np, vp, up (ion particle density, velocity, 
! and bulk velocity).   
!----------------------------------------------------------------------
!      include 'incurv.h'

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     up(nx,ny,nz,3),
     x     ndot(nx,ny,nz)

      real ndot_chex
!      real neutral_density
      real r,z1,y1,x1
!      real function ranf      

      integer flg           !flag for while loop
      real rnd              !random number

      integer cnt
      real delta_N
      integer dNion
      real vol

      integer recvbuf
      integer count
      count = 1

      call Neut_Center(cx,cy,cz)

      cnt = 0
      do i = 1,nx-1
         do j = 1,ny-1
            do k = 2,nz-1
               vol = dx*dy*dz_cell(k)

               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = gz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
!               ndot_chex = 0.
!               ndot_chex = nuin(i,j,k)*nf(i,j,k)
!               if (r .gt. 3.0*Rio) then
!c                  ndot_chex = nuin(i,j,k)*nf(i,j,k)
!                  ndot_chex = 0.0
!               endif

               if (r .gt. 2*RIo) then
                  delta_N = vol*beta*neutral_density(r)*dt/tau_photo
               endif
!               delta_N = vol*(ndot(i,j,k))*dt*beta



               if (delta_N .ge. 1.0) then 
!               write(*,*) 'delta_N...',my_rank,delta_N,i,j,k
!                  write(*,*) 'ndot...',delta_N
                  dNion = nint(delta_N)
                  l1 = Ni_tot + 1
                  Ni_tot = Ni_tot + dNion
                  do l = l1,Ni_tot
                     vp(l,:) = 0.0
                     
                     xp(l,1) = qx(i) + (pad_ranf())*dx_grid(i)
                     xp(l,2) = qy(j) + (pad_ranf())*dy_grid(j)
                     xp(l,3) = qz(k) + (pad_ranf())*dz_grid(k)
                     
!                     ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
!                     ijkp(l,2) = floor(xp(l,2)/dy)


                     ii=0
 17                  continue
                     ii = ii + 1
                     if (xp(l,1) .gt. qx(ii)) go to 17 !find i on non-uniform 
                     ii = ii-1
                     ijkp(l,1)= ii                     
                     
                     jj=0
 18                  continue
                     jj = jj + 1
                     if (xp(l,2) .gt. qy(jj)) go to 18 !find j on non-uniform 
                     jj = jj-1
                     ijkp(l,2)= jj

                     kk=0
 16                  continue
                     kk = kk + 1
                     if (xp(l,3) .gt. qz(kk)) go to 16 !find k on non-uniform 
                     kk = kk-1
                     ijkp(l,3)= kk

!                     kk=1
!                     do 16 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
!                        ijkp(l,3) = kk !grid
!                        kk=kk+1
! 16                  continue
!                     kk=ijkp(l,3)
!                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!                        ijkp(l,3) = kk+1
!                     endif
                     
                     mrat(l) = 1.0/m_pu
!                     m_arr(l) = mproton*m_pu
!                     Ni_tot = l
                     cnt = cnt + 1
                     do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 
     x                     0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /beta
!                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta

                     enddo                     
                  enddo

               endif
               if (delta_N .lt. 1.0) then 
                  if (delta_N .gt. pad_ranf()) then
                     l = Ni_tot + 1
                     vp(l,:) = 0.0
                     
                     xp(l,1) = qx(i) + (pad_ranf())*dx_grid(i)
                     xp(l,2) = qy(j) + (pad_ranf())*dy_grid(j)
                     xp(l,3) = qz(k) + (pad_ranf())*dz_grid(k)
                     
!                     ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
!                     ijkp(l,2) = floor(xp(l,2)/dy)


                     ii=0
 27                  continue
                     ii = ii + 1
                     if (xp(l,1) .gt. qx(ii)) go to 27 !find i on non-uniform 
                     ii = ii-1
                     ijkp(l,1)= ii                     
                     
                     jj=0
 28                  continue
                     jj = jj + 1
                     if (xp(l,2) .gt. qy(jj)) go to 28 !find j on non-uniform 
                     jj = jj-1
                     ijkp(l,2)= jj

                     kk=0
 29                  continue
                     kk = kk + 1
                     if (xp(l,3) .gt. qz(kk)) go to 29 !find k on non-uniform 
                     kk = kk-1
                     ijkp(l,3)= kk
                     
!                     kk=1
!                     do 18 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
!                        ijkp(l,3) = kk !grid
!                        kk=kk+1
! 18                  continue
!                     kk=ijkp(l,3)
!                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!                        ijkp(l,3) = kk+1
!                     endif

                     if (ijkp(l,3) .le. 0) then 
                     write(*,*) 'index error...',ijkp(l,3),qz(k),xp(l,3)
                     endif
                     
                     mrat(l) = 1.0/m_pu
!                     m_arr(l) = mproton*m_pu
                     Ni_tot = l
                     cnt = cnt + 1
                     do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 
     x                    0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /beta
!                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta
                     enddo                     
                  endif
               endif
            enddo
         enddo
      enddo


      write(*,*) 'total new ions...',cnt,Ni_tot,my_rank


!      call MPI_Barrier(MPI_COMM_WORLD,ierr)

!      call MPI_ALLREDUCE(cnt,recvbuf,count,
!     x     MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

!      write(*,*) 'total dNi....',recvbuf
!      stop

      
      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue
      

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)

      return
      end SUBROUTINE Ionize_Io
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_ndot(ndot)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl


      real recvbuf
      integer count
      count = 1


      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = gz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
!               theta = atan2(sqrt(x1**2 + y1**2),z1)
!               phi = atan2(y1,x1)
               
!               dvol = dx*dy*dz_grid(k)
! shell distribution
!               ndot(i,j,k) = dvol*exp(-(r - 1.4*RIo)**2/(0.3*RIo)**2)*
!     x                       sin(theta)
! power law distribution
!               if (r .gt. RIo) then
!                  ndot(i,j,k) = dvol*(Rio/r)**3.5*
!     x                 sin(theta)
!               endif
!               if (r .le. RIo) then 
!                  ndot(i,j,k) = 0.0
!               endif

               npofr = vol*beta*neutral_density(r)*dt/tau_photo
               ndot(i,j,k) = exp(-(r - 1.4*RIo)**2/(0.2*RIo)**2)*
     x                       sin(theta)*(cos(phi)+1)/2 !+
!     x                    0.2*dvol*exp(-(r - 1.2*RIo)**2/(0.1*RIo)**2)*
!     x                       sin(theta)
               ndot_intgl = ndot_intgl + ndot(i,j,k)
            enddo
         enddo
      enddo

      write(*,*) 'ndot_intgl....',ndot_intgl,my_rank
      

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

!      write(*,*) 'ndot_intgl...',ndot_intgl
!      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot
      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot/(dx*dy*delz)
      if (ndot_intgl .lt. 0.001) ndot = 0.0

! add a power law contribution to ndot

!      ndot_intgl = 0.0
!      do i = 1,nx
!         do j = 1,ny
!            do k = 1,nz
!               x1 = qx(i) - cx
!               y1 = qy(j) - cy
!               z1 = gz(k) - cz !global coordinate
!               r = sqrt(x1**2 + y1**2 + z1**2)
!               theta = atan2(sqrt(x1**2 + y1**2),z1)
!               phi = atan2(y1,x1)
               
!               dvol = dx*dy*dz_grid(k)
!               if (r .gt. RIo) then
!                  ndot2(i,j,k) = dvol*(Rio/r)**12*
!     x                 sin(theta)*(cos(phi+PI)+1)/2
!               endif
!               if (r .le. RIo) then 
!                  ndot2(i,j,k) = 0.0
!               endif

!               ndot_intgl = ndot_intgl + ndot2(i,j,k)*dvol
!            enddo
!         enddo
!      enddo

!      call MPI_Barrier(MPI_COMM_WORLD,ierr)

!      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
!     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

!      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

!      write(*,*) 'ndot_intgl...',ndot_intgl
!      ndot2 = (Mdot/(mO*recvbuf))*ndot2
!      if (ndot_intgl .lt. 0.001) ndot2 = 0.0

! combine shell with power law

!      ndot = ndot + ndot2

!      write(*,*) 'Mdot...',sum(ndot)*dx*dy*delz*mO,my_rank

      write(*,*) 'Mdot...',sum(ndot)*(dx*dy*delz)*m_pu*mproton,my_rank,
     x         (ndot(nx/2,ny/2,nz/2))

!      write(*,*) 'Max ndot...',maxval(ndot)
!      stop

      return
      end SUBROUTINE get_ndot
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_ndot_gauss(ndot)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl


      real recvbuf
      integer count
      count = 1


      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = qz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
               theta = atan2(sqrt(x1**2 + y1**2),z1)
               phi = atan2(y1,x1)
               dvol = dx*dy*dz_grid(k)
               ndot(i,j,k) = exp(-(r)**2/(RIo)**2)
               ndot_intgl = ndot_intgl + ndot(i,j,k)
            enddo
         enddo
      enddo

      write(*,*) 'ndot_intgl....',ndot_intgl,my_rank
      

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot/(dx*dy*delz)
      if (ndot_intgl .lt. 0.01) ndot = 0.0

      write(*,*) 'Mdot...',sum(ndot)*(dx*dy*delz)*m_pu*mproton,my_rank,
     x         (ndot(nx/2,ny/2,nz/2))

      return
      end SUBROUTINE get_ndot_gauss
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_ndot_Xianzhe(ndot,nn)
!     ndot is a density rate (cm-3 s-1)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)
      real nn(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl

!      real RIo_off !now defined as a parameter in para.h
!      parameter (RIo_off = 0.0)

      real recvbuf
      integer count
      count = 1

      if(my_rank.eq.0) then
          write(*,*)'Scaled ionization rate XIANZHE (maind/) '
          write(*,*)'----------------------------------------------'
          endif

      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = qz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
!               theta = atan2(sqrt(x1**2 + y1**2),z1)
!               phi = atan2(y1,x1)
               dvol = dx*dy*dz_grid(k)

               if ((r.ge.(Rio+RIO_off)) .and. r .lt. 3*Rio)   then  ! to compare to Linker values

! try Linker again
!                   ndot(i,j,k)= 5.e6*(r/Rio)**(-3.5)*1.e15/(25.*3600.)  ! no pwl to make it smooth

                     ndot(i,j,k) =  nn(i,j,k)/ (25. * 3600.) !  timescale for Ioniz = 25 hours
                     ndot(i,j,k) = ndot(i,j,k) *49./64. !  to scale to Xianzhe 7e27 s-1
                      endif

               if (r.lt.(Rio+RIO_off)) then
                     ndot(i,j,k)=0.
                     endif

               if( (r.ge.(Rio+RIO_off)).and.(r.le.(6.*Rio))  ) then  ! to compare to Linker values
                     ndot_intgl = ndot_intgl + ndot(i,j,k)*dvol !cm-3 integrated on a cell volume
                     endif
            enddo
         enddo
      enddo

      write(*,*) '    Proc ionioz rate XIANZHE(in 6Rio)= ',ndot_intgl,my_rank


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(my_rank.eq.0) then
      write(*,*) '   TOTAL XIANZHE (6 Rio) ndot_intgl_global ',recvbuf,my_rank
      endif

!     DOLS I remove the scaling to Mdot from para.h
! Linker give a rate and that is it (kg cm-3 s-1) but the code wants only km-3s-1
      return
      end SUBROUTINE get_ndot_Xianzhe
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_nuin(nuin,nn,uf,nf,ndot)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real nuin(nx,ny,nz)
      real nn(nx,ny,nz)
      real uf(nx,ny,nz,3)
      real nf(nx,ny,nz)
      real ndot(nx,ny,nz)
      real ufc(nx,ny,nz,3) !gather at cell center

!      real sigma_in
!      parameter (sigma_in = 3.0e-26)

!      call periodic(uf)
!      call periodic_scalar(nf)

      call face_to_center(uf,ufc)

!      call periodic(ufc)
      
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
!               uf2(1) = 0.5*(uf(i,j,k,1) + uf(i-1,j,k,1)) 
!               uf2(2) = 0.5*(uf(i,j,k,2) + uf(i,j-1,k,2)) 
!               uf2(3) = 0.5*(uf(i,j,k,3) + uf(i,j,k-1,3)) 
!               nuin(i,j,k) = sqrt(uf2(1)**2 + uf2(2)**2 + 
!     x                       uf2(3)**2)*sigma_in*nn(i,j,k)
               nuin(i,j,k) = sqrt(ufc(i,j,k,1)**2 + ufc(i,j,k,2)**2 + 
     x                       ufc(i,j,k,3)**2)*sigma_in*nn(i,j,k)
!               nuin(i,j,k) = 10.0*ndot(i,j,k)/nf(i,j,k)
            enddo
         enddo
      enddo

      call periodic_scalar(nuin)

      return
      end SUBROUTINE get_nuin
!----------------------------------------------------------------------



!----------------------------------------------------------------------
      SUBROUTINE Ionize_pluto_mp(np,np_2,vp,vp1,xp,m_tstep,input_p,up)
! Ionizes the neutral cloud with a 28 s time constant and fill particle
! arrays, np, vp, up (ion particle density, velocity, 
! and bulk velocity).   
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real np(nx,ny,nz),
     x     np_2(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)
!     x     gz(nz)

      real function ranf      
      real uptot1(3),uptot2(3)
!      integer*4 dNi         !# of new born ions created in dt
!      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
!      real r_xyz       !radial distance
!      real src(200)         !particle source distribution
!      real Nofr(200)        !number of neutrals as func of r
      real nnofr       !neutral density vs. r
      real npofr       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      real ddni
      integer cnt, l1
!      real neutral_density
      real npmax

!      integer*4 ion_cnt(nx,ny,nz)  !keeps running count of ions 
                                   !in cell for calculating the bulk
                                   !flow velocity
!      real vsum(nx,ny,nz,3) !total particle velocity in cell of the
                            !new born ions...for get bulk flow

      call Neut_Center(cx,cy,cz)

! get source density

      
      vol = dx**3
      cnt = 0
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               r = sqrt((qx(i)-cx)**2 + (qy(j)-cy)**2 + (gz(k)-cz)**2)
!               r = sqrt((xp(l,1)-cx)**2 + (xp(l,2)-cy)**2 + 
!     x              (gz(ijkp(l,3))-cz)**2) !use global coords

               !np = electron density, used for recombination

!               npmax = (-k_rec*np(i,j,k) + sqrt((k_rec*np(i,j,k))**2 + 
!     x              4*k_rec*neutral_density(r)/tau_photo))/(2*k_rec)
               
             
               npmax = sqrt(neutral_density(r)/(tau_photo*k_rec))

!               if ((np_2(i,j,k) .gt. npmax) .and.
!     x            (np_2(i,j,k) .gt. 0.0)) then
!                  write(*,*) 'limit ion production...',np_2(i,j,k),
!     x               npmax
!               endif
               if ((r .le. dx*S_radius) .and.
!               if ((r .le. dx*4) .and.
     x              (np_2(i,j,k) .lt. npmax)) then

!                   write(*,*) 'npmax..',npmax,np_2(i,j,k)

!                  write(*,*) 'r...',r/(dx*S_radius)
!                  write(*,*) 'npofr...',r
!                  nnofr = Qo/(4*pi*r**2*vrad)
!                  npofr = vol*beta*nnofr*dt/tau_photo
               if (r .le. dx*4) then
                  bpu = beta_pu
                  npofr = vol*beta*bpu*
     x                 neutral_density(r)*dt/tau_photo
               else 
                  bpu = 1.0
                  npofr = vol*beta*
     x                    neutral_density(r)*dt/tau_photo
               endif
!                  write(*,*) 'npofr...',npofr
                  if ((npofr .ge. 1) .and. (npofr+l1 .lt. Ni_max)) then
                     do ll = 1,nint(npofr)
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        

                        xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                        xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                        xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)

                        ii=0
 26                     continue
                        ii = ii + 1
                        if (xp(l,1) .gt. qx(ii)) go to 26 !find i on non-uniform 
                        ii = ii-1
                        ijkp(l,1)= ii

!                        ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
!                        ijkp(l,2) = floor(xp(l,2)/dy)
                    
                        jj=0
 18                     continue
                        jj = jj + 1
                        if (xp(l,2) .gt. qy(jj)) go to 18 !find j on non-uniform 
                        jj = jj-1
                        ijkp(l,2)= jj

                        kk=2
 15                     continue
                        kk = kk + 1
                        if (xp(l,3) .gt. qz(kk)) go to 15 !find k on non-uniform 
                        kk = kk-1
                        ijkp(l,3)= kk

!                        kk=1
!                       do 15 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find ck
!                           ijkp(l,3) = kk !grid
!                           kk=kk+1
! 15                     continue
!                        kk=ijkp(l,3)
!                        if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!                           ijkp(l,3) = kk+1
!                        endif
                        
                        mrat(l) = 1.0/m_pu
!                        m_arr(l) = mproton*m_pu
                        beta_p(l) = bpu
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                      0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x                          (beta*beta_p(l))
                    input_p(m) = input_p(m) + (mion/mrat(l))*vp(l,m)/
     x                          (beta*beta_p(l))
                        enddo                     

                     enddo
                  endif
                  if ((npofr .lt. 1).and.(npofr + l1 .lt. Ni_max)) then
                     if (npofr .gt. pad_ranf()) then
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        
                        
                        xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                        xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                        xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)


                        ii=0
 27                     continue
                        ii = ii + 1
                        if (xp(l,1) .gt. qx(ii)) go to 27 !find i on non-uniform 
                        ii = ii-1
                        ijkp(l,1)= ii

!                        ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
!                        ijkp(l,2) = floor(xp(l,2)/dy)
                        

                        jj=0
 17                     continue
                        jj = jj + 1
                        if (xp(l,2) .gt. qy(jj)) go to 17 !find i on non-uniform 
                        jj = jj-1
                        ijkp(l,2)= jj                     

                        kk=2
 16                     continue
                        kk = kk + 1
                        if (xp(l,3) .gt. qz(kk)) go to 16 !find k on non-uniform 
                        kk = kk-1
                        ijkp(l,3)= kk

!                        kk=1
!                       do 16 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find ck
!                           ijkp(l,3) = kk !grid
!                           kk=kk+1
! 16                     continue
!                        kk=ijkp(l,3)
!c                        if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!                          ijkp(l,3) = kk+1
!                        endif
                        
                        mrat(l) = 1.0/m_pu
!                        m_arr(l) = mproton*m_pu
                        beta_p(l) = bpu
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                          0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2/
     x                          (beta*beta_p(l))
                           input_p(m)=input_p(m)+(mion/mrat(l))*vp(l,m)/
     x                          (beta*beta_p(l))
                        enddo                     
                        
                        
                     endif

                  endif

!                  vol_shell = (4./3.)*pi*((r+dx/2)**3 - (r-dx/2)**3)
!                  vol_shell_min = (4./3.)*pi*(2.0*Rp)**3
!                  if (vol_shell .lt. vol_shell_min) then 
!                     vol_shell = vol_shell_min
!                  endif
!                  vol = dx**3
!c                  dNi_shell = (Qo*beta)*dx*dt/(vrad*tau_photo)
!c                  dNi_cell = dNi_shell*vol/vol_shell
!c                  write(*,*) 'dNi_shell...',cart_rank,
!c     x                  dNi_shell*vol/vol_shell
!                  rnd = pad_ranf()
!                  if (dNi_shell*vol/vol_shell .gt. rnd) then
!                     l = Ni_tot + 1
!                     vp(l,1) = 0.0
!                     vp(l,2) = 0.0
!                     vp(l,3) = 0.0
!c                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
!c                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
!c                     xp(l,3) = qz(k) + (0.5-pad_ranf())*dz_grid(k)

!                     xp(l,1) = qx(i) + (pad_ranf())*dx
!                     xp(l,2) = qy(j) + (pad_ranf())*dy
!                     xp(l,3) = qz(k) + (pad_ranf())*dz_grid(k)

!                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!                     ijkp(l,2) = nint(xp(l,2)/dy)

!c                     kk=1
!c                     do 50 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
!c                        ijkp(l,3) = kk !grid
!c                        kk=kk+1
!c 50                  continue
!c                     kk=ijkp(l,3)
!c                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!c                        ijkp(l,3) = kk+1
!c                     endif

!                     kk=1
!                     do 15 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
!                        ijkp(l,3) = kk !grid
!                        kk=kk+1
! 15                  continue
!                     kk=ijkp(l,3)
!                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
!                        ijkp(l,3) = kk+1
!                     endif

!                     mrat(l) = 1.0/m_pu
!                     m_arr(l) = mproton*m_pu
!                     Ni_tot = l
!                     cnt = cnt + 1
!                     do 45 m=1,3
!                        vp1(l,m) = vp(l,m)
!                        input_E = input_E + 
!     x                       0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
!                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta
! 45                  continue                     
!                  endif                     
               endif
            enddo
         enddo
      enddo

      
!      write(*,*) 'total new ions....',my_rank,cnt,Ni_tot         
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      stop

      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            write(*,*) 'removing ion...',my_rank,ijkp(l,:)
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

!      stop

!      write(*,*) 'Ni_tot after wake....',Ni_tot

!      call get_interp_weights(xp)
!      call update_np(np)
!      call update_up(vp,np,up)


      
      return
      end SUBROUTINE Ionize_pluto_mp
!----------------------------------------------------------------------

      end MODULE chem_rates
