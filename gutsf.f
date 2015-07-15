
      MODULE gutsf

      USE global
      USE boundary
      USE grid_interp
      
      contains


!----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(b1,b12,b1p2,bt,b0)
! loops run 1 to n since values are only being copied
!----------------------------------------------------------------------
!      include 'incurv.h'
      
      real b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(nx,ny,nz,3)
!     x     bdp(nx,ny,nz,3)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               bt(i,j,k,1) = b1p2(i,j,k,1) + b0(i,j,k,1)
               bt(i,j,k,2) = b1p2(i,j,k,2) + b0(i,j,k,2)
               bt(i,j,k,3) = b1p2(i,j,k,3) + b0(i,j,k,3) 
               do 10 m=1,3
!     uf2(i,j,k,m) = uf(i,j,k,m)
                  b12(i,j,k,m) = b1(i,j,k,m)
                  b1(i,j,k,m) = b1p2(i,j,k,m)
 10            continue
               
!               call obstacle_boundary_B(b0,bt)
!               call obstacle_boundary_B(b0,b12)
!               call obstacle_boundary_B(b0,b1)

      
      return
      end SUBROUTINE f_update_tlev
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      SUBROUTINE crossf(aa,bbmf,cc)
! The cross product is formed at the main cell center.  a is assumed
! be main cell contravarient (cell face) and b is assumed to be
! main cell covarient (cell edge).  The result is main cell
! contravarient (cell face).

! Can only center vectors on main cell for loops 3 to n...and can
! only extrapolate back for loops 2 to n-1.  Must handle other cases
! separately.

! The magnetic field does not require extrapolation to cell centers
! on boundaries since dB/dx = 0 is the boundary condition.  That is
! just copy the interior values to the boundary.
!----------------------------------------------------------------------
CVD$R VECTOR

!      include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real bbmf(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)

      real ax,ay,az,bx,by,bz    !dummy vars
      real temp                 !used to vectorize loop
      real zfrc(nz)             !0.5*dz_grid(k)/dz_cell(k)
!      real ct(nx,ny,nz,3)       !temp main cell center cross product
      real aac(3),bbc(3)


! extrapolate(/interpolate) to main cell center and do cross product


      call periodic(aa)
      call periodic(bbmf)


      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue



      do 10 k=2,nz-1      
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               im = i-1         !assume daa/dxyz = 0 at boundary
               jm = j-1         !bbmf is given on boundary
               km = k-1

               ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
               bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

               ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
               by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

               az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
               bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

               ct(i,j,k,1) = ay*bz - az*by
               ct(i,j,k,2) = az*bx - ax*bz
               ct(i,j,k,3) = ax*by - ay*bx

 10            continue

       call periodic(ct)

! extrapolate back to main cell contravarient positions.
! ...just average across cells since cell edges are centered
! about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      call periodic(cc)


      return
      end SUBROUTINE crossf
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE crossf2(aa,btc,cc)
! The cross product is formed at the main cell center.  aa and btc must
! be given already extrapolated to the main cell center.
!----------------------------------------------------------------------
      !include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real btc(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)

      real ax,ay,az,bx,by,bz    !dummy vars


      call periodic(aa)
      call periodic(btc)

      do 10 k=2,nz-1      
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               im = i-1    
               jm = j-1    
               km = k-1

               ax = aa(i,j,k,1) 
               bx = btc(i,j,k,1)

               ay = aa(i,j,k,2)
               by = btc(i,j,k,2)

               az = aa(i,j,k,3)
               bz = btc(i,j,k,3)

               ct(i,j,k,1) = ay*bz - az*by
               ct(i,j,k,2) = az*bx - ax*bz
               ct(i,j,k,3) = ax*by - ay*bx

 10            continue

       call periodic(ct)

! extrapolate back to main cell contravarient positions.
! ...just average across cells since cell edges are centered
! about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      call periodic(cc)


      return
      end SUBROUTINE crossf2
!----------------------------------------------------------------------



!----------------------------------------------------------------------
      SUBROUTINE cov_to_contra(bt,btmf)
! Converts total magnetic field from main cell covarient positions
! to main cell contravarient positions.  This is then used in the
! fluid velocity update routines.  This routine assumes that cell 
! edges and cell centers are "topologically centered".  So the grid
! points do not reside at the cell centers...rather they are offset
! a little so that the cell edges are equidistant from the k and k-1
! grid points.  In extrapolating the coventient vectors to the 
! contravarient vector positions, this assymetry is accounted for
! using a linear interpolation of the k and k-1 values to the grid
! point location.
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real bt(nx,ny,nz,3),   !main cell covarient
     x     btmf(nx,ny,nz,3)  !main cell contravarient

      real bx1, bx2, by1, by2, bz1, bz2  !main cell center fields
      real zrat           !ratio for doing linear interpolation
                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b_j, b_jm, b_i, b_im !intermediate step in average process

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1
               im = i-1
               jm = j-1
               km = k-1

! The x component of B resides at the k and k-1 edges, so this
! requires the non-uniform grid interpolation

               zplus = (qz(k+1) + qz(k))/2.0
               zminus = (qz(k) + qz(k-1))/2.0
               zrat = (qz(k) - zminus)/(zplus - zminus)
   
               b_j = bt(i,j,km,1) 
     x               + zrat*(bt(i,j,k,1) - bt(i,j,km,1)) 
               b_jm = bt(i,jm,km,1)
     x                + zrat*(bt(i,jm,k,1) - bt(i,jm,km,1))
               bx1 = (b_j + b_jm)/2.0

               b_j = bt(ip,j,km,1) 
     x               + zrat*(bt(ip,j,k,1) - bt(ip,j,km,1)) 
               b_jm = bt(ip,jm,km,1)
     x                + zrat*(bt(ip,jm,k,1) - bt(ip,jm,km,1))
               bx2 = (b_j + b_jm)/2.0

               
               b_i = bt(i,j,km,2) 
     x               + zrat*(bt(i,j,k,2) - bt(i,j,km,2)) 
               b_im = bt(im,j,km,2)
     x                + zrat*(bt(im,j,k,2) - bt(im,j,km,2))           
               by1 = (b_i + b_im)/2.0

               b_i = bt(i,jp,km,2) 
     x               + zrat*(bt(i,jp,k,2) - bt(i,jp,km,2)) 
               b_im = bt(im,jp,km,2)
     x                + zrat*(bt(im,jp,k,2) - bt(im,jp,km,2))
               by2 = (b_i + b_im)/2.0


               bz1 = 0.25*(bt(i,j,k,3) + bt(i,jm,k,3) +
     x                     bt(im,jm,k,3) + bt(im,j,k,3))
               bz2 = 0.25*(bt(i,j,kp,3) + bt(i,jm,kp,3) +
     x                     bt(im,jm,kp,3) + bt(im,j,kp,3))

               btmf(i,j,k,1) = 0.5*(bx1+bx2)
               btmf(i,j,k,2) = 0.5*(by1+by2)
               btmf(i,j,k,3) = 0.5*(bz1+bz2)

 10            continue

!      call boundaries(btmf)
      call periodic(btmf)

      return
      end SUBROUTINE cov_to_contra
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE curlB2(b1,np,aj)
! Calculates curl B / n*alpha.  The resulting "current" is called aj
! which is used in several other places in the code.  This curl is 
! performed on the main cell where B is covarient.  The resulting
! current is main cell contravarient.  Note that dz_cell is used for
! the cell dimensions since dz_grid is not equal to dz_cell on non-
! uniform grid.
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real b1(nx,ny,nz,3),
!     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)

      real curl_B(3)      !dummy for holding curl vector
      real ntot(3)        !total density, np + nf

!      call periodic_scalar(np)
!      call periodic_scalar(nf)
      call periodic(b1)
!c     call fix_normal_b(b1)

      do 10 k=2,nz-1   
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

!               if (ip .gt. nx) then ip = nx
!               if (jp .gt. ny) then jp = ny
!               if (kp .gt. nz) then kp = nz

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
!     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
!     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
!     x                 + 0.5*(np(i,j,k)+np(i,j,kp))


               ntot(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(np(i,j,k)+np(i,j,kp))

               curl_B(1) = (b1(i,j,k,3)/dy) - (b1(i,j-1,k,3)/dy) 
     x                    + (b1(i,j,k-1,2)/dz_cell(k))  
     x                    - (b1(i,j,k,2)/dz_cell(k))
               curl_B(2) = (b1(i,j,k,1)/dz_cell(k)) 
     x                     - (b1(i,j,k-1,1)/dz_cell(k))
     x                     - (b1(i,j,k,3)/dx) + (b1(i-1,j,k,3)/dx)
               curl_B(3) = (b1(i,j,k,2)/dx) - (b1(i-1,j,k,2)/dx) 
     x                     + (b1(i,j-1,k,1)/dy) - (b1(i,j,k,1)/dy)

               do 10 m=1,3
                  aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)

 10            continue

!      call periodic(aj)

      return
      end SUBROUTINE curlB2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE curlB(b1,np,aj)
! Calculates curl B / n*alpha.  The resulting "current" is called aj
! which is used in several other places in the code.  This curl is 
! performed on the main cell where B is covarient.  The resulting
! current is main cell contravarient.  Note that dz_cell is used for
! the cell dimensions since dz_grid is not equal to dz_cell on non-
! uniform grid.
!----------------------------------------------------------------------
CVD$R VECTOR
      !include 'incurv.h'

      real b1(nx,ny,nz,3),
!     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)

      real curl_B(3)      !dummy for holding curl vector
      real ntot(3)        !total density, np + nf

!      call periodic_scalar(np)
!      call periodic_scalar(nf)
      call periodic(b1)
!c     call fix_normal_b(b1)

      do 10 k=2,nz-1   
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

!               if (ip .gt. nx) then ip = nx
!               if (jp .gt. ny) then jp = ny
!               if (kp .gt. nz) then kp = nz

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
!     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
!     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
!     x                 + 0.5*(np(i,j,k)+np(i,j,kp))


               ntot(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(np(i,j,k)+np(i,j,kp))

               curl_B(1) = (b1(i,j,k,3) - 
     x              b1(i,j-1,k,3))/dy_cell(j) +
     x              (b1(i,j,k-1,2) - b1(i,j,k,2))/dz_cell(k)
               curl_B(2) = (b1(i,j,k,1) - 
     x              b1(i,j,k-1,1))/dz_cell(k) +
     x              (b1(i-1,j,k,3) - b1(i,j,k,3))/dx_cell(i)
               curl_B(3) = (b1(i,j,k,2) - 
     x              b1(i-1,j,k,2))/dx_cell(i) + 
     x              (b1(i,j-1,k,1) - b1(i,j,k,1))/dy_cell(j)


               do 10 m=1,3
                  aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)
 10            continue

!      call periodic(aj)

      return
      end SUBROUTINE curlB
!----------------------------------------------------------------------




!----------------------------------------------------------------------
      SUBROUTINE curlE2(E,curl_E)
! E is dual cell covarient, and curl_E will be returned as main
! cell covarient...as all magnetic fields are.  All i,j,k exclude
! boundaries.  Boundaries are taken care of in main fluid code.
!----------------------------------------------------------------------
!      include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges

!      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               lx = qx(i+1) - qx(i)
               ly = qy(j+1) - qy(j)
               lz = qz(k+1) - qz(k)

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/ly) - (E(i,j,k,3)/ly)
     x                       + (E(i,j,k,2)/lz) - (E(i,j,k+1,2)/lz)
               curl_E(i,j,k,2) =  (E(i,j,k,3)/lx) - (E(i+1,j,k,3)/lx)
     x                       + (E(i,j,k+1,1)/lz) - (E(i,j,k,1)/lz)
               curl_E(i,j,k,3) =  (E(i,j,k,1)/ly) - (E(i,j+1,k,1)/ly)
     x                       + (E(i+1,j,k,2)/lx) - (E(i,j,k,2)/lx)

 10          continue

!      call periodic(curl_E)

      return
      end SUBROUTINE curlE2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE curlE(E,curl_E)
! E is dual cell covarient, and curl_E will be returned as main
! cell covarient...as all magnetic fields are.  All i,j,k exclude
! boundaries.  Boundaries are taken care of in main fluid code.
!----------------------------------------------------------------------
      !include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges

!      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

!               lx = qx(i+1) - qx(i)
!               ly = qy(j+1) - qy(j)
!               lz = qz(k+1) - qz(k)

!               lx = dx_grid(i)
!               ly = dy_grid(j)
!               lz = dz_grid(k)

!               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/dy_grid(j)) - 
!     x              (E(i,j,k,3)/dy_grid(j))
!     x              + (E(i,j,k,2)/dz_grid(k)) - 
!     x              (E(i,j,k+1,2)/dz_grid(k))
!               curl_E(i,j,k,2) =  (E(i,j,k,3)/dx_grid(i)) - 
!     x              (E(i+1,j,k,3)/dx_grid(i))
!     x              + (E(i,j,k+1,1)/dz_grid(k)) - 
!     x              (E(i,j,k,1)/dz_grid(k))
!               curl_E(i,j,k,3) =  (E(i,j,k,1)/dy_grid(j)) - 
!     x              (E(i,j+1,k,1)/dy_grid(j))
!     x              + (E(i+1,j,k,2)/dx_grid(i)) - 
!     x              (E(i,j,k,2)/dx_grid(i))


               curl_E(i,j,k,1) =  (E(i,j+1,k,3)- 
     x              E(i,j,k,3))/dy_grid(j)
     x              + (E(i,j,k,2)- 
     x              E(i,j,k+1,2))/dz_grid(k)
               curl_E(i,j,k,2) =  (E(i,j,k,3) - 
     x              E(i+1,j,k,3))/dx_grid(i)
     x              + (E(i,j,k+1,1) - 
     x              E(i,j,k,1))/dz_grid(k)
               curl_E(i,j,k,3) =  (E(i,j,k,1) - 
     x              E(i,j+1,k,1))/dy_grid(j)
     x              + (E(i+1,j,k,2) - 
     x              E(i,j,k,2))/dx_grid(i)



 10          continue

!      call periodic(curl_E)

      return
      end SUBROUTINE curlE
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_E(E,b0,bt,aj,up,np,nu)
! E must be at time level m. We have uf at levels m-1/2 and m+1/2, so
! the average value is used for uf in the calculation of ui.
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
!     x     btmf(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
!     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
!     x     gradP(nx,ny,nz,3)

      real btc(nx,ny,nz,3)

      real ntot(3)         !total density np + nf
      real fnp(3),fnf(3)   !fraction np and nf of n
      real npave(3)

!      real a(nx,ny,nz,3), 
!     x     c(nx,ny,nz,3)  !dummy vars for doing cross product

      real aa(nx,ny,nz,3)

!      call periodic_scalar(np)
!      call periodic_scalar(nf)

      do 10 k=2,nz-1    
         do 10 j=2,ny-1
            do 10 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!               if (ip .gt. nx) then ip = nx
!               if (jp .gt. ny) then jp = ny
!               if (kp .gt. nz) then kp = nz

!               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
!               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
!               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

!               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
!               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
!               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
               
!               fnp(1) = npave(1)/ntot(1)
!               fnp(2) = npave(2)/ntot(2)
!               fnp(3) = npave(3)/ntot(3)

!               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
!               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
!               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

!               ntot = np(i,j,k) + nf(i,j,k)
!               fnp = np(i,j,k)/ntot
!               fnf = nf(i,j,k)/ntot

               do 10 m=1,3
                  a(i,j,k,m) = aj(i,j,k,m) - up(i,j,k,m)
!                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
!     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
!                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
!     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
 10               continue




!      call crossf(a,btmf,c)
      call face_to_center(a,aa)
      call edge_to_center(bt,btc)
      call crossf2(aa,btc,c)


      do 20 k=2,nz-1      
         do 20 j=2,ny-1   
            do 20 i=2,nx-1
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
!     x                         + nuei*aj(i,j,k,m) !- gradP(i,j,k,m)
!     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue


!      call fix_tangential_E(E)
      call periodic(E)

!      call obstacle_boundary(E)

!      call fix_tangential_E(E)

!      E(nx-1:nx,:,:,3) = -vsw*b0(nx-1:nx,:,:,2)
!      E(nx-1:nx,:,:,2) = -vsw*b0(nx-1:nx,:,:,3)
!      E(:,:,nz,2) = 0.0
!      E(:,:,nz,1) = 0.0

!      E(:,:,1,3) = -vsw*q*b0(:,:,1,2)/mO
!      E(:,:,1,2) = -vsw*q*b0(:,:,1,3)/mO
!      E(:,:,1,2) = 0.0
!      E(:,:,1,1) = 0.0


!      E(nx-1:nx,:,:,3) = vsw*q*b0_init/mO

!      E(nx-1:nx,:,:,2) = 0.0
!c      E(nx-1:nx,:,:,1) = 0.0

      return
      end SUBROUTINE get_E
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu)
! Predictor step in magnetic field update.
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
!     x     btmf(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
!     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
!     x     gradP(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)   !curl of E

!      call cov_to_contra(bt,btmf) 
!      call edge_to_center(bt,btc)
!      call get_E(E,b0,bt,btmf,aj,up,np,nu)  !E at time level m 
      call get_E(E,b0,bt,aj,up,np,nu)  !E at time level m 

      call curlE(E,curl_E)
!      call fix_tangential_E(E)
!      call periodic(E)
!      call fix_tangential_E(E)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

!                  b1p2(i,j,k,m)=lww1*(b12(i+1,j,k,m)+
!     x                 b12(i-1,j,k,m)+
!     x                 b12(i,j+1,k,m)+b12(i,j-1,k,m)+
!     x                 b12(i,j,k+1,m)+b12(i,j,k-1,m))+
!     x                 lww2*b12(i,j,k,m) -
!     x                 2.0*dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b12(i,j,k,m) - 
     x                            2.0*dtsub*curl_E(i,j,k,m)
 10               continue

!      call boundaries(b1p2)
!      call damp(b1p2)
      call periodic(b1p2)
!      call obstacle_boundary_B(b0,b1p2)
!      call fix_normal_b(b1p2)

      return
      end SUBROUTINE predict_B
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)
! The main feature here is that E must be calculated at time level
! m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
! calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
! m + 1/2, so they are used as is. 
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
!     x     gradP(nx,ny,nz,3),
!     x     bdp(nx,ny,nz,3)

      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real btc(nx,ny,nz,3) 
      real aa(nx,ny,nz,3) 
    

      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)

!      real a(nx,ny,nz,3),
!     x     c(nx,ny,nz,3)    !dummy vars for doing cross product


      do 5 k=1,nz
         do 5 j=1,ny
            do 5 i=1,nx
               btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1)) +
     x              b0(i,j,k,1)
               b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2)) + 
     x              b0(i,j,k,2)
               b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               btp1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3)) + 
     x              b0(i,j,k,3)
               b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
 5             continue

      call curlB(b1p1,np,aj)
!      call obstacle_boundary_B(b0,b1p1)

!      call periodic_scalar(np)
!      call periodic_scalar(nf)

      do 10 k=2,nz-1       
         do 10 j=2,ny-1
            do 10 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!               if (ip .gt. nx) then ip = nx
!               if (jp .gt. ny) then jp = ny
!               if (kp .gt. nz) then kp = nz

!               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
!               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
!               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

!               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
!               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
!               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
               
!               fnp(1) = npave(1)/ntot(1)
!               fnp(2) = npave(2)/ntot(2)
!               fnp(3) = npave(3)/ntot(3)

!               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
!               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
!               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

!               ntot = np(i,j,k) + nf(i,j,k)
!               fnp = np(i,j,k)/ntot
!               fnf = nf(i,j,k)/ntot

               do 10 m=1,3
!                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
!     x                                       fnf(m)*uf(i,j,k,m)
                  a(i,j,k,m) = aj(i,j,k,m) - up(i,j,k,m)
!                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
!     x                           fnf(m)*uf(i,j,k,m)
 10               continue

!      call cov_to_contra(btp1,btp1mf)
      call face_to_center(a,aa)
!      call edge_to_center(btp1,btc)
      call edge_to_face(btp1,btp1mf)   !add shift to cell face for smoothing
      call face_to_center(btp1mf,btc)
!      call face_to_center(btp1mf,btc)

!       call crossf(a,btp1mf,c)
      call crossf2(aa,btc,c)
       
      do 20 k=2,nz-1       
         do 20 j=2,ny-1     
            do 20 i=2,nx-1  
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
!     x                         + nuei*aj(i,j,k,m) !- gradP(i,j,k,m)
!     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue

!      call fix_tangential_E(E)
      call periodic(E)
!      call obstacle_boundary(E)
!      call fix_tangential_E(E)

!      E(nx-1:nx,:,:,3) = -vsw*b0(nx-1:nx,:,:,2)
!      E(nx-1:nx,:,:,2) = -vsw*b0(nx-1:nx,:,:,3)

!      E(:,:,nz,1) = 0.0

!      E(:,:,1,3) = -vsw*q*b0(:,:,1,2)/mO
!      E(:,:,1,2) = 0.0
!      E(:,:,1,1) = 0.0

!      E(nx-1:nx,:,:,3) = vsw*q*b0_init/mO
!      E(nx-1:nx,:,:,2) = 0.0
!c      E(nx-1:nx,:,:,1) = 0.0

      return
      end SUBROUTINE get_Ep1
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,np,nu)
! Corrector step in magnetic field update.
!----------------------------------------------------------------------
CVD$R VECTOR
!      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
!     x     gradP(nx,ny,nz,3),
!     x     bdp(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)            !curl of E

      call get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)  
                                                   !E at time level m 

      call curlE(E,curl_E)
!      call fix_tangential_E(E)
!      call periodic(E)
!      call fix_tangential_E(E)

!      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

!                  b1p2(i,j,k,m)=lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+
!     x                 b1(i,j+1,k,m)+b1(i,j-1,k,m)+
!     x                 b1(i,j,k+1,m)+b1(i,j,k-1,m))+
!     x                 lww2*b1(i,j,k,m) -
!     x                 dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b1(i,j,k,m) - 
     x                            dtsub*curl_E(i,j,k,m)
 10               continue

!      call boundaries(b1p2)
!      call damp(b1p2)
      call periodic(b1p2)
!      call obstacle_boundary_B(b0,b1p2)
!      call fix_normal_b(b1p2)

      return
      end SUBROUTINE correct_B
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      SUBROUTINE Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     E(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     pup(3),
     x     puf(3),
     x     peb(3),
     x     input_p(3)

      real vol
      real mom_flux
      real exb(nx,ny,nz,3)
      real npave(3),nfave(3)


      call crossf2(E,b1,exb)

      do 5 m=1,3
         pup(m) = 0
         puf(m) = 0
         peb(m) = 0
 5              continue

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               ip = i+1
               jp = j+1
               kp = k+1
               if (ip .eq. nx) then ip = nx-1
               if (jp .eq. ny) then jp = ny-1
               if (kp .eq. nz) then kp = nz-1
               vol = dx*dy*dz_cell(k)
               npave(1) = 0.5*(np(i,j,k) + np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k) + np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k) + np(i,j,kp))
               nfave(1) = 0.5*(nf(i,j,k) + nf(ip,j,k))
               nfave(2) = 0.5*(nf(i,j,k) + nf(i,jp,k))
               nfave(3) = 0.5*(nf(i,j,k) + nf(i,j,kp))
               do 10 m=1,3
!                  pup(m) = pup(m) + npave(m)*vol*mBa*up(i,j,k,m)
                  pup(m) = pup(m) + np(i,j,k)*vol*mBa*up(i,j,k,m)
!                  puf(m) = puf(m) + nfave(m)*vol*mO*uf(i,j,k,m)
                  puf(m) = puf(m) + nf(i,j,k)*vol*mO*uf(i,j,k,m)
                  peb(m) = peb(m) + epsilon*1e3*exb(i,j,k,m)*vol*(mO/q)
 10               continue


      return
      end SUBROUTINE Momentum_diag
!----------------------------------------------------------------------


!----------------------------------------------------------------      
      SUBROUTINE check_time_step(bt,np)
!----------------------------------------------------------------      
      
      real bt(nx,ny,nz,3)
      real np(nx,ny,nz)

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               ak = 2./dx
               btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
     x              bt(i,j,k,3)**2)
               a1 = ak**2*Btot/(alpha*(np(i,j,k)))
               a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
               womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
               phi = womega/ak
               deltat = dx/phi
               if(deltat .le. 2.0*dtsub) then 
                  write(*,*) 'time stepping error...',i,j,k
                  dtsub = dtsub/2.0
                  ntf = ntf*2.0
!     if (mindt .gt. deltat) then
!     deltat = mindt
!     write(*,*) 'mindt...',mindt
!     endif
!     do while (2.0*dtsub .gt. deltat)
!     dtsub = dtsub/2.0
!     ntf = ntf*2.0
!     write(*,*) 'Changing subcycle time step...',dtsub,deltat,ntf
!     enddo
               endif
               
!     if(deltat .gt. 2.0*dtsub) then 
!     write(*,*) 'time stepping error...'
!     dtsub = 2.0*dtsub
!     ntf = ntf/2.0
!     write(*,*) 'subcycle time steps...',ntf,dtsub/dt
!     endif
            enddo
         enddo
      enddo

      if (ntf .gt. 100) then 
         ntf = 100
      endif
      
!      write(*,*) 'subcycle time steps...',ntf,dtsub
      

      return
      end subroutine check_time_step
!----------------------------------------------------------------      
      end MODULE gutsf
