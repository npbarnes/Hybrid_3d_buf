
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




!c----------------------------------------------------------------------
!      SUBROUTINE get_um_dot_BB(u,b,cc)
!c uf and btmf are gathered at main cell center and uf.B*B 
!c calculated.  Result returned to main cell contravarient
!c postion.
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
!     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
!     x     cc(nx,ny,nz,3)   !(uf.B)*B

!      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
!      real temp                !used to vectorize loop
!c      real ct(nx,ny,nz,3)      !result are main cell center
!      real udotb               !u dot b
!      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

!! first gather everything at center

!      call periodic(u)
!      call periodic(b)
!c      call fix_normal_b(b)

!      do 5 k=1,nz
!         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
! 5       continue

!      do 10 k=2,nz-1      
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1

!               im = i-1
!               jm = j-1     
!               km = k-1

!               ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
!               bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))

!               uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
!               by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))

!               uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
!               bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)

!               udotb = ux*bx + uy*by + uz*bz

!               ct(i,j,k,1) = udotb*bx
!               ct(i,j,k,2) = udotb*by
!               ct(i,j,k,3) = udotb*bz

! 10            continue

!      call periodic(ct)


!c extrapolate back to main cell contravarient positions.
!c ...just average across cells.

!      do 60 k=2,nz-1
!         do 60 j=2,ny-1
!            do 60 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
!               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
!               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

! 60            continue

!      call periodic(cc)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_um_dot_BB_old(u,b,cc)
!c uf and btmf are gathered at main cell center and uf.B*B 
!c calculated.  Result returned to main cell contravarient
!c postion.
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
!     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
!     x     cc(nx,ny,nz,3)   !(uf.B)*B

!      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
!      real temp                !used to vectorize loop
!      real ct(nx,ny,nz,3)      !result are main cell center
!      real udotb               !u dot b
!      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

!! first gather everything at center

!      call periodic(u)
!      call periodic(b)
!c      call fix_normal_b(b)

!      do 5 k=1,nz
!         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
! 5       continue

!      do 10 i=2,nx-1      
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1

!               im = i-1
!               jm = j-1     
!               km = k-1

!               if (i .eq. 2) then 
!                  ux = 2.0*u(2,j,k,1) - 
!     x                 0.5*(u(3,j,k,1) + u(2,j,k,1))
!                  bx = 2.0*b(2,j,k,1) - 
!     x                 0.5*(b(3,j,k,1) + b(2,j,k,1))
!               else
!                  ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
!                  bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))
!               endif

!               if (j .eq. 2) then 
!                  uy = 2.0*u(i,2,k,2) - 
!     x                 0.5*(u(i,3,k,2) + u(i,2,k,2))
!                  by = 2.0*b(i,2,k,2) - 
!     x                 0.5*(b(i,3,k,2) + b(i,2,k,2))
!               else
!                  uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
!                  by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))
!               endif

!               if (k .eq. 2) then
!                  uz = 2.0*u(i,j,2,3) - 
!     x                 zfrc(k)*(u(i,j,3,3) - u(i,j,2,3)) + u(i,j,2,3) 
!c                  uz = 2.0*u(i,j,2,3) - 
!c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
!                  bz = b(i,j,2,3)
!c                  bz = 2.0*b(i,j,2,3) - 
!c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
!               else
!                  uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
!                  bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)
!c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
!c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))            
!               endif

!               udotb = ux*bx + uy*by + uz*bz

!               ct(i,j,k,1) = udotb*bx
!               ct(i,j,k,2) = udotb*by
!               ct(i,j,k,3) = udotb*bz

! 10            continue

!      call periodic(ct)

!c extrapolate back to main cell contravarient positions.
!c ...just average across cells.

!      do 60 i=2,nx
!         do 60 j=2,ny
!            do 60 k=2,nz

!               ip = i+1
!               jp = j+1
!               kp = k+1

!               if (i .eq. nx-1) then 
!                  cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
!     x                           0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
!               else
!                  cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
!               endif

!               if (j .eq. ny-1) then 
!                  cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
!     x                           0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
!               else
!                  cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
!               endif
                  
!               if (k .eq. nz-1) then
!c                  temp = 2.0*ct(i,j,nz,3) - 
!c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
!                   temp = 0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) +
!     x                   (2.0*dz_cell(nz)/dz_grid(nz))*(ct(i,j,nz,3) -
!     x                    0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)))
!               else
!                  temp = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
!               endif             
  
!               cc(i,j,k,3) = temp

! 60            continue

!      call periodic(cc)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_ugradu_Lax(uf,ugradu,delta_t)
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real uf(nx,ny,nz,3),
!     x     ugradu(nx,ny,nz,3)

!      real ufc(nx,ny,nz,3)
!      real ax1,ax2,ay1,ay2,az1,az2       !misc const
!      real u1,u2,u3                      !temp vars

!      parameter(ad = 0.001)                 !coefficient to add extra
                                         !diffusion
!      call periodic(uf)

!      call face_to_center(uf,ufc)

!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               ip=i+1
!               jp=j+1
!               kp=k+1
!               im=i-1
!               jm=j-1
!               km=k-1

!c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!               ax1 = 0.5*ufc(i,j,k,1)/dx
!               ax2 = ad*abs(ufc(ip,j,k,1) - ufc(im,j,k,1))
!               ay1 = 0.5*ufc(i,j,k,2)/dy
!               ay2 = ad*abs(ufc(ip,j,k,2) - ufc(im,j,k,2))
!               u1 = ax1*(ufc(ip,j,k,1) - ufc(im,j,k,1)) - 
!     x              ax2*(ufc(im,j,k,1) - 2.0*ufc(i,j,k,1) +
!     x                      ufc(ip,j,k,1))
!               u2 = ay1*(ufc(i,jp,k,1) - ufc(i,jm,k,1)) - 
!     x              ay2*(ufc(i,jm,k,1) - 2.0*ufc(i,j,k,1) +
!     x                      ufc(i,jp,k,1)) 
!               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
!               az2 = ad*abs(ufc(ip,j,k,3) - ufc(im,j,k,3))
!               u3 = az1*(ufc(i,j,kp,1)-ufc(i,j,km,1)) -
!     x              az2*(ufc(i,j,km,1) - 2.0*ufc(i,j,k,1) +
!     x                   ufc(i,j,kp,1))
!               ct(i,j,k,1) = u1 + u2 + u3

!c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

!               ax1 = 0.5*ufc(i,j,k,1)/dx
!               ax2 = ad*abs(ufc(i,jp,k,1) - ufc(i,jm,k,1))
!               ay1 = 0.5*ufc(i,j,k,2)/dy
!               ay2 = ad*abs(ufc(i,jp,k,2) - ufc(i,jm,k,2))
!               u1 = ax1*(ufc(ip,j,k,2) - ufc(im,j,k,2)) - 
!     x              ax2*(ufc(im,j,k,2) - 2.0*ufc(i,j,k,2) +
!     x                      ufc(ip,j,k,2))
!               u2 = ay1*(ufc(i,jp,k,2) - ufc(i,jm,k,2)) - 
!     x              ay2*(ufc(i,jm,k,2) - 2.0*ufc(i,j,k,2) +
!     x                      ufc(i,jp,k,2)) 
!               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
!               az2 = ad*abs(ufc(i,jp,k,3) - ufc(i,jm,k,3))
!               u3 = az1*(ufc(i,j,kp,2)-ufc(i,j,km,2)) -
!     x              az2*(ufc(i,j,km,2) - 2.0*ufc(i,j,k,2) +
!     x                   ufc(i,j,kp,2))
!               ct(i,j,k,2) = u1 + u2 + u3

!c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
!               ax1 = 0.5*ufc(i,j,k,1)/dx
!               ax2 = ad*abs(ufc(i,j,kp,1) - ufc(i,j,km,1))
!               ay1 = 0.5*ufc(i,j,k,2)/dy
!               ay2 = ad*abs(ufc(i,j,kp,2) - ufc(i,j,km,2))
!               u1 = ax1*(ufc(ip,j,k,3) - ufc(im,j,k,3)) - 
!     x              ax2*(ufc(im,j,k,3) - 2.0*ufc(i,j,k,3) +
!     x                      ufc(ip,j,k,3))
!               u2 = ay1*(ufc(i,jp,k,3) - ufc(i,jm,k,3)) - 
!     x              ay2*(ufc(i,jm,k,3) - 2.0*ufc(i,j,k,3) +
!     x                      ufc(i,jp,k,3)) 
!               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
!               az2 = ad*abs(ufc(i,j,kp,3) - ufc(i,j,km,3))
!               u3 = az1*(ufc(i,j,kp,3)-ufc(i,j,km,3)) -
!     x              az2*(ufc(i,j,km,3) - 2.0*ufc(i,j,k,3) +
!     x                   ufc(i,j,kp,3))
!               ct(i,j,k,3) = u1 + u2 + u3

! 10            continue

!      call periodic(ct)

!c interpolate back to contravarient positions.

!      do 20 k=2,nz-1
!         do 20 j=2,ny-1
!            do 20 i=2,nx-1
!               ip = i+1
!               jp = j+1
!               kp = k+1
!               ugradu(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(ip,j,k,1))
!               ugradu(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,jp,k,2))
!               ugradu(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,kp,3))
! 20         continue

!       call periodic(ugradu)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE get_Ef(Ef,aj,np,nf,up,uf,btmf,nu,ugradu,delta_t,
!     x                  gradPf)
!c Need to treat boundaries separately!!
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real Ef(nx,ny,nz,3),
!     x     aj(nx,ny,nz,3),
!     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
!     x     up(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
!     x     btmf(nx,ny,nz,3),
!     x     nu(nx,ny,nz),
!     x     ugradu(nx,ny,nz,3),
!     x     gradPf(nx,ny,nz,3)
 
!      real ntot(3)                !total plasma density
!      real fnp(3)                 !fraction, np/n
!      real aac(3),bbc(3),ccc(3)
!      real cc(nx,ny,nz,3)

!      call periodic_scalar(np)
!      call periodic_scalar(nf)

!      do 10 k=2,nz-1 
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!c               if (ip .gt. nx) then ip = nx
!c               if (jp .gt. ny) then jp = ny
!c               if (kp .gt. nz) then kp = nz

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
!     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
!     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
!     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

!               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
!               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
!               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
               
!               do 10 m=1,3
!                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m)
! 10               continue

!      call crossf(a,btmf,c)

!      call get_ugradu_Lax(uf,ugradu,delta_t)

!      do 20 k=2,nz-1 
!         do 20 j=2,ny-1
!            do 20 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!c               if (ip .gt. nx) then ip = nx
!c               if (jp .gt. ny) then jp = ny
!c               if (kp .gt. nz) then kp = nz

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
!     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
!     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
!     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

!               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
!               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
!               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

!               do 20 m=1,3
!                  Ef(i,j,k,m) = c(i,j,k,m) - ugradu(i,j,k,m) 
!     x                          + nu(i,j,k)*fnp(m)*up(i,j,k,m)
!     x                          + nuei*aj(i,j,k,m) - gradPf(i,j,k,m)
!c     x                          + etar(i,j,k,m)*aj(i,j,k,m)
!c                  Ef(i,j,k,m) = c(i,j,k,m) + 
!c     x                          nu(i,j,k)*fnp(m)*up(i,j,k,m)
! 20            continue

!      call periodic(Ef)
!c      call fix_tangential_E(Ef)

!      Ef(nx-1:nx,:,:,3) = 0.0
!      Ef(nx-1:nx,:,:,2) = 0.0
!      Ef(nx-1:nx,:,:,1) = 0.0

!      return
!      end
!c----------------------------------------------------------------------



!c----------------------------------------------------------------------
!      SUBROUTINE get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
!     x                            delta_t)
!c This is the heart of the fluid velocity update.  It solves eqn. 18
!c (Dan's paper) for uf+
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real Ef(nx,ny,nz,3),
!     x     btmf(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3), !using particle velocities at t level n-1/2
!     x     nu(nx,ny,nz),
!     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
!     x     uplus(nx,ny,nz,3),
!     x     uminus(nx,ny,nz,3)

!      real a1,b1,b2,b3      !4 coefficients (see notebook solution)
!      real PP,QQ            !intermediate variables for a1,b1,b2,b3
!      real eta              !intermediate variable for PP,QQ
!c      real B(3)             !B for cross product call
!c      real Bsqrd                     !B*B
!      real um_x_B(nx,ny,nz,3)        !uf- X B
!      real um_dot_BB(nx,ny,nz,3)     !uf- . B
!      real ntot(3)                   !total density np + nf
!      real npave(3)
!      real btc(nx,ny,nz,3)
!      real bsqrd(nx,ny,nz),bsq(3)

!      do 10 i=2,nx-1    
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1
!               do 10 m=1,3
!                  uminus(i,j,k,m) = uf(i,j,k,m) + 
!     x                              0.5*delta_t*Ef(i,j,k,m)
! 10            continue

!      call crossf(uminus, btmf, um_x_B)
!      call get_um_dot_BB(uminus , btmf, um_dot_BB)

!      call face_to_center(btmf,btc)

!      do 15 k=2,nz-1 
!         do 15 j=2,ny-1
!            do 15 i=2,nx-1
!               bsqrd(i,j,k) =  btc(i,j,k,1)**2 + btc(i,j,k,2)**2 + 
!     x               btc(i,j,k,3)**2
! 15            continue

!      call periodic_scalar(bsqrd)
              
!      do 20 k=2,nz-1
!         do 20 j=2,ny-1
!            do 20 i=2,nx-1

!               ip = i+1
!               jp = j+1
!               kp = k+1

!c               if (ip .gt. nx) then ip = nx
!c               if (jp .gt. ny) then jp = ny
!c               if (kp .gt. nz) then kp = nz

!               bsq(1) = 0.5*(bsqrd(i,j,k) + bsqrd(ip,j,k))
!               bsq(2) = 0.5*(bsqrd(i,j,k) + bsqrd(i,jp,k))
!               bsq(3) = 0.5*(bsqrd(i,j,k) + bsqrd(i,j,kp))

!               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
!               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
!               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

!               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + npave(1)
!               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + npave(2)
!               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + npave(3)

!               do 20 m=1,3

!                  eta = npave(m)*delta_t/(2.0*ntot(m))
!                  QQ = eta/(1.0+nu(i,j,k)*eta)
!                  PP = (1.0-nu(i,j,k)*eta)/(1.0+nu(i,j,k)*eta)
!                  a1 = 1.0/(1.0 + QQ*QQ*bsq(m))
!                  b1 = PP - (QQ*QQ*bsq(m))
!                  b2 = (QQ*PP) + QQ
!                  b3 = (QQ*QQ*PP) + (QQ*QQ)

!                  uplus(i,j,k,m) = a1*(b1*uminus(i,j,k,m) + 
!     x                             b2*um_x_B(i,j,k,m) + 
!     x                             b3*um_dot_BB(i,j,k,m))

! 20            continue

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus, 
!     x                      uminus,ugradu,up,gradP,nuin,bdp,pf1)
!c Calculate the fluid velocity, uf,  at the new time step and replace
!c uf1 with the new value, uf, in preparation for the next time step.
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real Ef(nx,ny,nz,3),
!     x     b0(ny),
!     x     b1(nx,ny,nz,3),
!     x     b12(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
!     x     uf2(nx,ny,nz,3),
!     x     ufp2(nx,ny,nz,3),
!     x     nu(nx,ny,nz),
!     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
!     x     uplus(nx,ny,nz,3), 
!     x     uminus(nx,ny,nz,3),
!     x     ugradu(nx,ny,nz,3),
!     x     up(nx,ny,nz,3),
!c     x     gradP(nx,ny,nz,3),
!     x     nuin(nx,ny,nz),
!     x     bdp(nx,ny,nz,3),
!     x     pf1(nx,ny,nz)

!      real b1h(nx,ny,nz,3)
!      real bth(nx,ny,nz,3)
!      real btmfh(nx,ny,nz,3)
!      real ajh(nx,ny,nz,3)
!      real gradPf(nx,ny,nz,3)

!      real m_den
!      real delta_t

!      delta_t = 2.0*dtsub

!      do k=2,nz-1
!         do j=2,ny-1
!            do i=2,nx-1
!               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
!               gradPf(i,j,k,1) = (pf1(i+1,j,k)-pf1(i,j,k))/(dx*m_den)
!               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
!               gradPf(i,j,k,2) = (pf1(i,j+1,k)-pf1(i,j,k))/(dy*m_den)
!               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
!               gradPf(i,j,k,3) = (pf1(i,j,k+1)-
!     x                             pf1(i,j,k))/(dz_grid(k)*m_den)
!            enddo
!         enddo
!      enddo


!      do 10 k=1,nz
!         do 10 j=1,ny
!            do 10 i=1,nx
!               b1h(i,j,k,1) = 0.5*(b1(i,j,k,1) + b12(i,j,k,1))
!               b1h(i,j,k,2) = 0.5*(b1(i,j,k,2) + b12(i,j,k,2))
!               b1h(i,j,k,3) = 0.5*(b1(i,j,k,3) + b12(i,j,k,3))
!               bth(i,j,k,1) = b1h(i,j,k,1) + bdp(i,j,k,1)
!               bth(i,j,k,2) = b1h(i,j,k,2) + b0(j) + bdp(i,j,k,2)
!               bth(i,j,k,3) = b1h(i,j,k,3) + bdp(i,j,k,3)
! 10            continue

!      call cov_to_contra(bth,btmfh)
!      call curlB(b1h,nf,np,ajh)

!      call get_Ef(Ef,ajh,np,nf,up,uf,btmfh,nu,ugradu,delta_t,gradPf)
!      call get_uplus_uminus(Ef,btmfh,uf2,nu,np,nf,uplus,uminus,
!     x                      delta_t)

!      do 20 k=2,nz-1
!         do 20 j=2,ny-1
!            do 20 i=2,nx-1
!               do 20 m=1,3
!                  ufp2(i,j,k,m) = uplus(i,j,k,m) + 
!     x                            0.5*delta_t*Ef(i,j,k,m) !-
!c     x                        0.5*delta_t*nuin(i,j,k)*uplus(i,j,k,m)
! 20            continue

!c      ufp2(nx-1:nx,:,:,1) = -vsw

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
!     x                      ugradu,aj,up,ufp1,gradP,nuin,pf)
!c Calculate the fluid velocity, uf,  at the new time step and replace
!c uf1 with the new value, uf, in preparation for the next time step.
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real Ef(nx,ny,nz,3),
!     x     btmf(nx,ny,nz,3),
!     x     uf(nx,ny,nz,3),
!     x     uf2(nx,ny,nz,3),
!     x     ufp2(nx,ny,nz,3),
!     x     nu(nx,ny,nz),
!     x     np(nx,ny,nz),
!     x     nf(nx,ny,nz),
!     x     uplus(nx,ny,nz,3), 
!     x     uminus(nx,ny,nz,3),
!     x     ugradu(nx,ny,nz,3),
!     x     aj(nx,ny,nz,3),
!     x     up(nx,ny,nz,3),
!     x     ufp1(nx,ny,nz,3),
!c     x     gradP(nx,ny,nz,3),
!     x     nuin(nx,ny,nz),
!     x     pf(nx,ny,nz)
 
!      real m_den
!      real delta_t
!      real gradPf(nx,ny,nz,3)

!      delta_t = dtsub

!      do k=2,nz-1
!         do j=2,ny-1
!            do i=2,nx-1
!               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
!               gradPf(i,j,k,1) = (pf(i+1,j,k)-pf(i,j,k))/(dx*m_den)
!               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
!               gradPf(i,j,k,2) = (pf(i,j+1,k)-pf(i,j,k))/(dy*m_den)
!               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
!               gradPf(i,j,k,3) = (pf(i,j,k+1)-
!     x                             pf(i,j,k))/(dz_grid(k)*m_den)
!            enddo
!         enddo
!      enddo

      
!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               do 10 m=1,3
!                  ufp1(i,j,k,m) = 0.5*(uf(i,j,k,m) + ufp2(i,j,k,m))
! 10              continue

!      call get_Ef(Ef,aj,np,nf,up,ufp1,btmf,nu,ugradu,delta_t,gradPf)
!      call get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
!     x                      delta_t)

!      do 20 k=2,nz-1
!         do 20 j=2,ny-1
!            do 20 i=2,nx-1
!               do 20 m=1,3
!                  uf2(i,j,k,m) = uf(i,j,k,m)
!                  uf(i,j,k,m) = uplus(i,j,k,m) + 0.5*dtsub*Ef(i,j,k,m)
!c     x                          - 0.5*dtsub*nuin(i,j,k)*uplus(i,j,k,m)
! 20            continue

!c      uf(nx-1:nx,:,:,1) = -vsw

!      return
!      end
!c----------------------------------------------------------------------


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


!c----------------------------------------------------------------------
!      SUBROUTINE predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real nf(nx,ny,nz),
!     x     nf1(nx,ny,nz),
!     x     nf3(nx,ny,nz),
!     x     nfp1(nx,ny,nz),
!     x     uf(nx,ny,nz,3),
!     x     divu(nx,ny,nz),
!     x     b1(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real flx(nx,ny,nz,3)
!      real ufc(nx,ny,nz,3)
!      real b1c(nx,ny,nz,3)

!c      call face_to_center(uf,ufc)
!c      call face_to_center(b1,b1c)
 
!      minnf = 10.0e20
!      maxnf = 10.0e20

!      call periodic_scalar(nf1)

!      do 5 i=2,nx-1
!         do 5 j=2,ny-1
!            do 5 k=2,nz-1
!               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
!               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
!               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
! 5             continue

!      call periodic(flx)

!      do 10 i=2,nx-1
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1
!               nfp1(i,j,k) = nf3(i,j,k) 
!     x         - (2.0*dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
!     x         - (2.0*dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
!     x         - (2.0*dtsub/(dz_cell(k)))*
!     x           (flx(i,j,k,3)-flx(i,j,k-1,3))
!               divu(i,j,k) = (uf(i,j,k,1) - uf(i-1,j,k,1))/dx +
!     x                       (uf(i,j,k,2) - uf(i,j-1,k,2))/dy +
!     x                       (uf(i,j,k,3) - uf(i,j,k-1,3))/dz_cell(k)
! 10            continue

!      do 50 i=2,nx-1
!         do 50 j=2,ny-1
!            do 50 k=2,nz-1
!               nf(i,j,k) = 0.5*(nfp1(i,j,k) + nf1(i,j,k))
!               nf3(i,j,k) = nf1(i,j,k)
! 50            continue

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE correct_nf(nf,nf1,ufp1)
!c----------------------------------------------------------------------

!      include 'incurv.h'

!      real nf(nx,ny,nz),
!     x     nf1(nx,ny,nz),
!     x     ufp1(nx,ny,nz,3) 

!      real flx(nx,ny,nz,3)
!      real minnf,maxnf
!      real ufp1c(nx,ny,nz,3)

!      minnf = 10.0000000000000000e20
!      maxnf = 10.0000000000000000e20

!c      call face_to_center(ufp1,ufp1c)

!      do 5 i=2,nx-1
!         do 5 j=2,ny-1
!            do 5 k=2,nz-1
!             flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
!             flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
!             flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
! 5           continue

!      call periodic(flx)

!      do 10 i=2,nx-1
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1
!               nf1(i,j,k) = nf1(i,j,k) 
!     x            - (dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
!     x            - (dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
!     x            - (dtsub/(dz_cell(k)))*
!     x              (flx(i,j,k,3)-flx(i,j,k-1,3))
!               if (nf(i,j,k) .lt. minnf) then minnf = nf(i,j,k)
!               if (nf(i,j,k) .ge. maxnf) then maxnf = nf(i,j,k)
! 10            continue


!      call periodic_scalar(nf1)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE trans_nf_Lax(nf,nf1,nfp1,uf)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real nf(nx,ny,nz),
!     x     nf1(nx,ny,nz),
!     x     nfp1(nx,ny,nz),
!     x     uf(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real flx(nx,ny,nz,3)
!      real ufc(nx,ny,nz,3)

!c      call face_to_center(uf,ufc)
 
!      minnf = 10.0e20
!      maxnf = 10.0e20

!      do 3 i=2,nx-1
!         do 3 j=2,ny-1
!            do 3 k=2,nz-1
!               nf1(i,j,k) = nfp1(i,j,k)
! 3             continue      

!      call periodic_scalar(nf1)

!      do 5 i=2,nx-1
!         do 5 j=2,ny-1
!            do 5 k=2,nz-1
!               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
!               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
!               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
! 5             continue

!      call periodic(flx)

!      do 10 i=2,nx-1
!         do 10 j=2,ny-1
!            do 10 k=2,nz-1
!               nfp1(i,j,k) = (1.0/6.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
!     x                               nf1(i,j+1,k)+nf1(i,j-1,k)+
!     x                               nf1(i,j,k+1)+nf1(i,j,k-1))  
!     x         - 0.5*(dtsub/dx)*(flx(i+1,j,k,1)-flx(i-1,j,k,1))
!     x         - 0.5*(dtsub/dy)*(flx(i,j+1,k,2)-flx(i,j-1,k,2))
!     x         - (dtsub/(dz_grid(k)+dz_grid(k+1)))*
!     x           (flx(i,j,k+1,3)-flx(i,j,k-1,3))
! 10            continue

!      do 50 i=2,nx-1
!         do 50 j=2,ny-1
!            do 50 k=2,nz-1
!               nf(i,j,k) = 0.5*(nf1(i,j,k) + nfp1(i,j,k))
! 50            continue

!      call periodic_scalar(nf)

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE trans_nf_LaxWend1(nf,nf1,nfp1,uf)
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real nf(nx,ny,nz),
!     x     nf1(nx,ny,nz),
!     x     nfp1(nx,ny,nz),
!     x     uf(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real flx(nx,ny,nz,3)

!      do 3 k=2,nz-1
!         do 3 j=2,ny-1
!            do 3 i=2,nx-1
!               nf1(i,j,k) = nfp1(i,j,k)
! 3             continue   
   
!      call periodic_scalar(nf1)

!      do 5 k=2,nz-1
!         do 5 j=2,ny-1
!            do 5 i=2,nx-1
!               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
!               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
!               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
! 5             continue

!      call periodic(flx)
!      call periodic_scalar(nf1)

!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               nf(i,j,k) = (1.0/12.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
!     x                     nf1(i,j+1,k)+nf1(i,j-1,k)+
!     x                     nf1(i,j,k+1)+nf1(i,j,k-1)+6.0*nf1(i,j,k))  
!     x         - 0.5*(dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
!     x         - 0.5*(dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
!     x         - 0.5*(dtsub/dz_grid(k))*
!     x           (flx(i,j,k,3)-flx(i,j,k-1,3))
! 10            continue

!      call periodic_scalar(nf)
!      nf(nx-1:nx,:,:) = nf_init*0.01

!      return
!      end
!c----------------------------------------------------------------------



!c----------------------------------------------------------------------
!      SUBROUTINE trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real nf(nx,ny,nz),
!     x     nf1(nx,ny,nz),
!     x     nfp1(nx,ny,nz),
!     x     ufp1(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real flx(nx,ny,nz,3)

!      call periodic_scalar(nf)

!      do 5 k=2,nz-1
!         do 5 j=2,ny-1
!            do 5 i=2,nx-1
!               flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
!               flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
!               flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
! 5             continue

!      call periodic(flx)
!      call periodic_scalar(nf1)

!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               nfp1(i,j,k) = nf1(i,j,k) +
!     x                     0.002*(nf1(i+1,j,k)+nf1(i-1,j,k)+
!     x                     nf1(i,j+1,k)+nf1(i,j-1,k)+
!     x                     nf1(i,j,k+1)+nf1(i,j,k-1)-6.0*nf1(i,j,k))
!     x         - (dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
!     x         - (dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
!     x         - (dtsub/dz_grid(k))*(flx(i,j,k,3)-flx(i,j,k-1,3))
! 10            continue

!      call periodic_scalar(nfp1)
!      nfp1(nx-1:nx,:,:) = nf_init*0.01

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE trans_pf_LaxWend1(pf,pf1,uf)
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real pf(nx,ny,nz),
!     x     pf1(nx,ny,nz),
!     x     uf(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real t1(nx,ny,nz,3)
!      real t1c(nx,ny,nz)
!      real t2(nx,ny,nz)

!      parameter(gamma = 5./3.)

!      call periodic_scalar(pf1)

!      do 5 k=2,nz-1
!         do 5 j=2,ny-1
!            do 5 i=2,nx-1
!               t1(i,j,k,1) = uf(i,j,k,1)*(pf1(i+1,j,k)-pf1(i,j,k))/dx 
!               t1(i,j,k,2) = uf(i,j,k,2)*(pf1(i,j+1,k)-pf1(i,j,k))/dy
!               t1(i,j,k,3) = 
!     x               uf(i,j,k,3)*(pf1(i,j,k+1)-pf1(i,j,k))/dz_grid(k)
!               t2(i,j,k) = 
!     x            gamma*pf1(i,j,k)*((uf(i,j,k,1)-uf(i-1,j,k,1))/dx +
!     x                   (uf(i,j,k,2)-uf(i,j-1,k,2))/dy +
!     x                   (uf(i,j,k,3)-uf(i,j,k-1,3))/dz_cell(k))
! 5             continue

!      call periodic(t1)

!      do k = 2,nz-1
!         do j = 2,ny-1
!            do i = 2,nx-1
!               t1c(i,j,k) = (1./2.)*(t1(i,j,k,1)+t1(i-1,j,k,1)) +
!     x                      (1./2.)*(t1(i,j,k,2)+t1(i,j-1,k,2)) +
!     x                      (1./2.)*(t1(i,j,k,3)+t1(i,j,k-1,3))
!            enddo
!         enddo
!      enddo

!c      call periodic_scalar(t1c)
!c      call periodic_scalar(t2)

!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               pf(i,j,k) = (1.0/12.0)*(pf1(i+1,j,k)+pf1(i-1,j,k)+
!     x                     pf1(i,j+1,k)+pf1(i,j-1,k)+
!     x                     pf1(i,j,k+1)+pf1(i,j,k-1)+6.0*pf1(i,j,k))  
!     x         - 0.5*dtsub*(t1c(i,j,k)+t2(i,j,k))
!               if (pf(i,j,k) .lt. 0.0) then
!                 ! write(*,*) 'error...',t1c(i,j,k),t2(i,j,k)
!               endif

! 10            continue

!c      pf = abs(pf)

!      call periodic_scalar(pf)
!      pf(nx-1:nx,:,:) = nf_init*0.01*kboltz*tempf0

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE trans_pf_LaxWend2(pf,pf1,ufp1)
!c----------------------------------------------------------------------
!CVD$R VECTOR
!      include 'incurv.h'

!      real pf(nx,ny,nz),
!     x     pf1(nx,ny,nz),
!     x     ufp1(nx,ny,nz,3)
 
!      real minnf,maxnf
!      real t1(nx,ny,nz,3)
!      real t1c(nx,ny,nz)
!      real t2(nx,ny,nz)
!      real pfp1(nx,ny,nz)

!      parameter(gamma = 5./3.)

!c      call periodic_scalar(pf1)

!      do 5 k=2,nz-1
!         do 5 j=2,ny-1
!            do 5 i=2,nx-1
!               t1(i,j,k,1) = ufp1(i,j,k,1)*(pf(i+1,j,k)-pf(i,j,k))/dx 
!               t1(i,j,k,2) = ufp1(i,j,k,2)*(pf(i,j+1,k)-pf(i,j,k))/dy
!               t1(i,j,k,3) = 
!     x               ufp1(i,j,k,3)*(pf(i,j,k+1)-pf(i,j,k))/dz_grid(k)
!               t2(i,j,k) = 
!     x           gamma*pf(i,j,k)*((ufp1(i,j,k,1)-ufp1(i-1,j,k,1))/dx +
!     x                   (ufp1(i,j,k,2)-ufp1(i,j-1,k,2))/dy +
!     x                   (ufp1(i,j,k,3)-ufp1(i,j,k-1,3))/dz_cell(k))
! 5             continue

!      call periodic(t1)

!      do k = 2,nz-1
!         do j = 2,ny-1
!            do i = 2,nx-1
!               t1c(i,j,k) = (1./2.)*(t1(i,j,k,1)+t1(i-1,j,k,1)) +
!     x                      (1./2.)*(t1(i,j,k,2)+t1(i,j-1,k,2)) +
!     x                      (1./2.)*(t1(i,j,k,3)+t1(i,j,k-1,3))
!            enddo
!         enddo
!      enddo

!c      call periodic_scalar(t1c)
!c      call periodic_scalar(t2)

!      do 10 k=2,nz-1
!         do 10 j=2,ny-1
!            do 10 i=2,nx-1
!               pfp1(i,j,k) =  
!     x                     0.002*(pf1(i+1,j,k)+pf1(i-1,j,k)+
!     x                     pf1(i,j+1,k)+pf1(i,j-1,k)+
!     x                     pf1(i,j,k+1)+pf1(i,j,k-1)-6.0*pf1(i,j,k))
!     x         +   pf1(i,j,k)
!     x         - dtsub*(t1c(i,j,k)+t2(i,j,k))
! 10            continue

!      pf1 = pfp1
!c      pf1 = abs(pf1)

!      call periodic_scalar(pf1)
!      pf1(nx-1:nx,:,:) = nf_init*0.01*kboltz*tempf0

!      return
!      end
!c----------------------------------------------------------------------


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

!      write(*,*) 'Momentum conservation...'
!      write(*,*) '  Particles.............',pup(1),pup(2),pup(3)
!      write(*,*) '  Fluid.................',puf(1),puf(2),puf(3)
!      write(*,*) '  ExB...................',peb(1),peb(2),peb(3)
!      write(*,*) '  Normalized............',
!     x                     (pup(1)+puf(1)+peb(1))/input_p(1),
!     x                     (pup(2)+puf(2)+peb(2))/input_p(2),
!     x                     (pup(3)+puf(3)+peb(3))/input_p(3)

! Momentum flux through boundary faces

! i = 2 face

!      do 20 j=2,ny
!         do 20 k=2,nz
!            m=1
!            i=2
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
! 20         continue

! i = nx face

!      do 30 j=2,ny
!         do 30 k=2,nz
!            m=1
!            i=nx
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
! 30         continue
!**********************
! j = 2 face

!      do 40 i=2,nx
!         do 40 k=2,nz
!            m=2
!            j=2
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
! 40         continue

! j = ny face

!      do 50 i=2,nx
!         do 50 k=2,nz
!            m=2
!            j=ny
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
! 50         continue
!****************
! k = 2 face

!      do 60 i=2,nx
!         do 60 j=2,ny
!            m=3
!            k=2
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
! 60         continue

! k = nz face

!      do 70 i=2,nx
!         do 70 j=2,ny
!            m=3
!            k=nz
!            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
!            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
!            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
! 70         continue

!      write(*,*) 'Normalized x momentum...',(pup(1)+puf(1))/input_p(1)
!      write(*,*) 'Normalized y momentum...',(pup(2)+puf(2))/input_p(2)
!      write(*,*) 'Normalized z momentum...',(pup(3)+puf(3))/input_p(3)

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













