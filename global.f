      MODULE GLOBAL


      USE inputs
      USE dimensions
      save

! contains the common data structures to be used by the various
! simulation subroutines 

!      include 'para.h'
!      include 'mpif.h'

! raw grid coordinate data
  
      real gz(nz)
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction 
      real dx_grid(nx),dx_cell(nx)
      real dy_grid(ny),dy_cell(ny)
      real xrat(nx), yrat(ny), zrat(nz)
      
      common /coords/ qx,qy,qz,gz,lambda,ri,rj,rk,dz_grid,dz_cell


! Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot, Ni_tot_sys, Ni_tot_sw, Ni_tot_buf,
     x          Ni_tot_out_buf       

! Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical in_bounds(Ni_max)
      logical in_bounds_buf(Ni_max_buf)
!      logical ionized(Ni_max)
!      real mix_ind(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, Ni_tot_sys, Ni_tot_sw,
     x       in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf

!      real seed
      integer*4 seed
      common /rndseed/ seed

! Indices for boundary cells
      integer ip,im,jp,jm,kp,km

! Total input energy (ions) 
      real input_E,input_Eb,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux, 
     x                input_chex, input_bill,input_Eb

! Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot

! Dummy variables for doing vector operations and particle flux
! calculation....etc.

      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct

! Variable time step for fluid velocity update used in both
! particle update and fluid update

      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub

! Weight variables for trilinear interpolation

      real wght(Ni_max,8)!, wquad(Ni_max,3)
      common /weights/ wght!, wquad

! variable for anomlous resistivity

!      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei

! variable for particle scale

      real beta, beta_p(Ni_max), dNi, dNi_sw
      common /scale/ beta, beta_p, dNi, dNi_sw


! mass array for multi-ion species
      real mrat,mrat_buf,beta_p_buf 
      common /mass/ mrat(Ni_max),
     x     mrat_buf(Ni_max_buf),beta_p_buf(Ni_max_buf)

! mixing array
      real mixed(nx,ny,nz)
      common /mix/ mixed


! parallel processor info
!      integer procnum,my_rank
!      common /procinfo/ procnum, my_rank

! parallel processor info
      integer procnum,my_rank,nbrs(4),cart_rank,cartcomm
      common /procinfo/ procnum, my_rank,nbrs,cart_rank,cartcomm

      end MODULE GLOBAL
