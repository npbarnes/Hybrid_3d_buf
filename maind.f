      PROGRAM MAIND
     
!----------------------------------------------------------------------
! maind.f
! Parallel version with no ion fluid, Nov 24, 2004
!----------------------------------------------------------------------

      
      USE global
      USE dimensions
      USE inputs
      USE mpi
      USE initial
      USE misc
      USE gutsp_dd
      USE gutsp_buf
      USE gutsf
      USE part_init
      USE grid_interp
      USE chem_rates

!      include 'incurv.h'

!----------------------------------------------------------------------
! Listing of all declared variables
!
! Note that the format for specifying variables at different time
! levels is based upon the number of 1/2 time steps that the varible
! is behind the current time level.  For example, 
! uf2 is 1 time step behind uf, ufp2 is 1 time step ahead of uf,
! and b12 is 1 full time step behind b1 (not 12 1/2 steps behind b 
! since b does not exist...right). b1p2 is an exception...this is a
! temporary holder for b1 at m+1 in the predictor/corrector update
! of the magnetic field.
!----------------------------------------------------------------------
!      integer time, t1, t2    !keep track of run time
!      external time
      save

      real b0(nx,ny,nz,3),            !ambient magnetic field
     x     b1(nx,ny,nz,3),    !1st order magnetic field
     x     b12(nx,ny,nz,3),   !b1 at previous time step
     x     b1p2(nx,ny,nz,3),  !temporary b1 at time level m+1
     x     bt(nx,ny,nz,3),    !total magnetic field..mc covarient
!     x     btmf(nx,ny,nz,3),  !main cell contravarient bt field
     x     btc(nx,ny,nz,3),   !btmf at cell center for particle move
!     x     bdp(nx,ny,nz,3),   !dipole magnetic field
!     x     nf(nx,ny,nz),      !ambient fixed fluid density
!     x     nf1(nx,ny,nz),     !nf at n-1/2
!     x     nf3(nx,ny,nz),     !nf at n-3/2
!     x     nfp1(nx,ny,nz),    !nf at n+1/2
!     x     nn(nx,ny,nz),      !neutral cloud density
!     x     nnd(nx,ny,nz),     !neutral cloud density decrement
     x     np(nx,ny,nz),      !particle ion den at time level n, n+1/2
!     x     np_1(nx,ny,nz),
     x     np_2(nx,ny,nz),
!     x     np3(nx,ny,nz,3),
     x     vp(Ni_max,3),      !particle velocity at t level n+1/2
     x     vp1(Ni_max,3),     !particle velocity at t level n
     x     vplus(Ni_max,3),   !v+ used in velocity update
     x     vminus(Ni_max,3),  !v- used in velocity update
     x     up(nx,ny,nz,3),    !particle flow at time level n, n+1/2
     x     xp(Ni_max,3),      !coordinates of ion particles
!     x     xp1(Ni_max,3),     !coordinates of ion particles at previous time step
!     x     uf(nx,ny,nz,3),    !fluid velocity
!     x     uf2(nx,ny,nz,3),   !fluid velcity at time level n-1
!     x     ufp1(nx,ny,nz,3),  !fluid velocity at time level n+1/2
!     x     ufp2(nx,ny,nz,3),  !fluid velocity at time level n+1
!     x     ui(nx,ny,nz,3),    !total ion flow velocity
     x     aj(nx,ny,nz,3),    !curlB/(alpha*n) 
     x     nu(nx,ny,nz),      !collision frequency
!     x     nuin(nx,ny,nz),    !ion-neutral collision frequency
     x     Ep(Ni_max,3),      !Ion particle electric field
!     x     Ef(nx,ny,nz,3),    !fluid electric field
     x     E(nx,ny,nz,3)     !electric field from electron mom eqn
!     x     uplus(nx,ny,nz,3), !u plus used in velocity update
!     x     uminus(nx,ny,nz,3),!u minus used in velocity update
!     x     pf(nx,ny,nz),      !fluid pressure at n
!     x     pf1(nx,ny,nz)      !fluid pressure at n-1/2


      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)
      real Ep_buf(Ni_max_buf,3)
      real vplus_buf(Ni_max_buf,3)
      real vminus_buf(Ni_max_buf,3)

      real xp_out_buf(Ni_max_buf,3)
      real vp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)
!      real m_arr_out_buf(Ni_max_buf)

      real temp_p(nx,ny,nz)
!     x     temp_p_1(nx,ny,nz),
!     x     temp_p_2(nx,ny,nz)

      real Evp,       !total particle kinetic energy
     x     Euf,       !total fluid kinetic energy
     x     EB1,       !total magnetic field energy
     x     EB1x,      !total b1x energy
     x     EB1y,      !total b1y energy
     x     EB1z,      !total b1z energy
     x     EE,        !total electric field energy
     x     EeP        !total electron pressure energy

      real pup(3),      !total particle momentum
     x     puf(3),      !total fluid momentum
     x     peb(3),      !total momentum carried by E and B fields
     x     input_p(3)   !input momentum

!      integer np_t_flg(Ni_max)
!      integer np_b_flg(Ni_max)
!      real np_t(nx,ny,nz)
!      real np_b(nx,ny,nz)

      real mr

      real chex_rate
      real bill_rate
      real satnp
!      real gradP(nx,ny,nz,3)
!      real etemp(nx,ny,nz)
!      real ugradu(nx,ny,nz,3)
!      real minnf,maxnf
!      real divu(nx,ny,nz)
      real mindt
      integer*4 t1,t2,cnt_rt
      real time
      integer ierr

      real ndot(nx,ny,nz)
      
      integer seedsize
      integer, dimension(:), allocatable :: seeder

      real recvbuf
      integer count
      
!      character filenum
      character(len=:), allocatable::filenum
      character(len=10) :: arg



!----------------------------------------------------------------------

      call readInputs()
      call initparameters()

!      stop

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

      filenum = int_to_str(my_rank+1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!      write(filenum,"(I2)") my_rank + 1
!      filenum = adjustl(filenum)
!      filenum = adjustr(adjustl(filenum))

      write(*,*) 'filenum...',filenum,my_rank+1

      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)

! create virtual topology (set dimensions in para.h)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periods, 
     x     reorder,cartcomm, ierr)

      call MPI_COMM_RANK(cartcomm, cart_rank, ierr)
      call MPI_CART_COORDS(cartcomm, cart_rank, cart_dims, cart_coords, 
     x                     ierr)
      call MPI_CART_SHIFT(cartcomm,0,1,nbrs(n_up),nbrs(n_down),ierr)
      call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(n_left), nbrs(n_right), 
     &     ierr)

      call system_clock(t1,cnt_rt)
!      seed = float(t1)

!----------------------------------------------------------------------
! Initialize all variables
!----------------------------------------------------------------------
      write(*,*) 'initializing variables...'

      Ni_tot = Ni_tot_0
      Ni_tot_sw = Ni_tot
!      Ni_tot_sys = Ni_tot*procnum
      Ni_tot_sys = Ni_tot
      print *,'Ni_tot_sys, Ni_tot..',Ni_tot_sys,Ni_tot,Ni_tot_sw

      if (my_rank .eq. 0) then
         call check_inputs(my_rank)
         write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*ny*nz)
         write(*,*) ' '
      endif

!      stop

!      Ni_tot = 6
      mstart = 0
      ndiag = 0
      prev_Etot = 1.0
      nuei = 0.0

! initialize seed for each processor

!      call random_seed
!      call random_seed(size = seedsize)
!      allocate(seeder(seedsize))
!      do n = 0,procnum-1 
!         if (my_rank .eq. n) then 
!            call random_seed(get=seeder)
!            call random_seed(put=seeder)
!         endif
!      enddo

      seed = t1 +my_rank*100

      ! Make sure all ranks are initialized and sychronized before attempting
      ! to make a system call
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
      
      call get_command_argument(number=1,value=arg,status=ierr)
      if(ierr .eq. 0 .and. trim(arg) .eq. "debug") then
            call debug_random_initialize()
      else 
            call random_initialize()
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 


      if (.not.(restart)) then
         do 66 i=1,nx
            do 66 j=1,ny
               do 66 k=1,nz
!                  pf(i,j,k) = nf_init*0.05*kboltz*tempf0
!                  pf1(i,j,k) = nf_init*0.05*kboltz*tempf0
!                  nf(i,j,k) = nf_init*0.0
!                  nf1(i,j,k) = nf_init*0.05  
!                  nf3(i,j,k) = nf_init*0.05 
!                  nfp1(i,j,k) = nf_init*0.05  
                  input_E = 0.0
                  input_p = 0.0
                  input_chex = 0.0
                  input_bill = 0.0
 66               continue
               endif

!      do 68 i = 1,nx
!         do 68 j = 1,ny
!            do 68 k = 1,nz
!               uf(i,j,k,1) = -vsw
!               uf2(i,j,k,1) = -vsw
!               ufp1(i,j,k,1) = -vsw
!               ufp2(i,j,k,1) = -vsw
! 68            continue

!      Ni_tot = 4000000


      if (.not.(restart)) then
!         m_arr(1:Ni_tot) = mproton
!         m_arr(Ni_tot+1:) = m_pu*mproton !mass N_2+ = 28.0
         mrat(1:Ni_tot) = 1.0
         mrat(Ni_tot+1:) = 1.0/m_pu !mass N_2+ = 28.0
         beta_p(1:Ni_tot) = 1.0
         beta_p(Ni_tot+1:) = beta_pu
      endif

      call grd8()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu)

!      call obstacle_boundary_nu(nu)

      call get_beta()

!         input_E = 0.0
!      do i = 1,nx
!         do j = 1,ny
!            do k = 1,nz
!               input_E = input_E + 
!     x          0.5*dx*dy*dz_grid(k)*nf_init*0.01*mO*(vsw*km_to_m)**2
!            enddo
!         enddo
!      enddo


      if (.not.(restart)) then
!      call sw_part_setup_temp(np,vp,vp1,xp,input_p,up)
!      call sw_part_setup_maxwl(np,vp,vp1,xp,xp1,input_p,up,np_t_flg,
!     x                         np_b_flg)
         call sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up)

         call part_setup_buf(xp_buf,vp_buf)
         
         call part_setup_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf,b0)
                  
!         call get_ndot(ndot)
!         call predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu) 
!         call correct_B(b0,b1,b1p2,E,aj,up,np,nu)

         call f_update_tlev(b1,b12,b1p2,bt,b0)
      endif

!----------------------------------------------------------------------




!----------------------------------------------------------------------
! check for restart flag
!----------------------------------------------------------------------
      write(*,*) 'restart status....',restart
      if (restart) then 
         write(*,*) 'opening restart.vars......'
         open(210,file='restart.vars'//filenum,
     x        status='unknown',
     x        form='unformatted')
         
         read(210)  b0,b1,b12,b1p2,bt,btc,np,np3,
     x        up,aj,nu,E,input_E,input_p,mstart,input_EeP,prev_Etot
     x        Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP

!         open(210,file='restart.vars',status='unknown',
!     x            form='unformatted')
         write(*,*) 'reading restart.vars......'
!         read(210) b0,b1,b12,b1p2,bt,btmf,btc,np,np3,vp,vp1,vplus,
!     x        vminus,up,xp,xp1,aj,nu,Ep,E,Evp,Euf,
!     x        EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,Ni_tot_sys,
!     x        ijkp,mstart,input_p,input_EeP,prev_Etot,
!     x        mrat,m_arr,xp_buf,vp_buf,Ep_buf,vplus_buf,
!     x        vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
!     x        B_out_buf,mrat_out_buf,m_arr_out_buf,
!     x        in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
!     x        m_arr_buf,mrat_buf
!         write(*,*) 'restarting hybrid.....'

!         if (my_rank .ge. 0) then 
          open(211,file='restart.part'//filenum,
     x           status='unknown',form='unformatted')
          read(211) vp,vp1,vplus,vminus,
     x         xp,Ep,Ni_tot,
     x         Ni_tot_sys,ijkp,
     x         mrat,
     x         xp_buf,vp_buf,Ep_buf,vplus_buf,
     x         vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
     x         B_out_buf,mrat_out_buf,
     x         in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
     x         mrat_buf



!          open(211,file='restart.part'//filenum,
!     x            status='unknown',form='unformatted')
!          read(211)  b0,b1,b12,b1p2,bt,btmf,btc,np,np3,
!     x         vp,vp1,vplus,vminus,
!     x         up,xp,xp1,aj,nu,Ep,E,input_E,Ni_tot,
!     x         Ni_tot_sys,ijkp,mstart,input_p,input_EeP,prev_Etot,
!     x         mrat,m_arr,
!     x         xp_buf,vp_buf,Ep_buf,vplus_buf,
!     x         vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
!     x         B_out_buf,mrat_out_buf,m_arr_out_buf,
!     x         in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
!     x         m_arr_buf,mrat_buf
!c         endif
      endif

      close(210)
      close(211)



!      write(*,*) 'mstart...',mstart
!      stop

!      write(*,*) 'restart status....',restart
!      if (restart) then 
!         write(*,*) 'opening restart.vars......'
!         open(210,file='restart.vars',status='unknown',
!     x            form='unformatted')
!         write(*,*) 'reading restart.vars......'
!         read(210) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,vplus,vminus,
!     x            up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
!     x            EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
!     x            ijkp,mstart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1,
!     x            input_chex,input_bill,pf,pf1,mrat,m_arr
!         write(*,*) 'restarting hybrid.....'

!         if (my_rank .gt. 0) then 
!          open(211,file='restart.part'//filenum,
!     x            status='unknown',form='unformatted')
!          read(211) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,
!     x              ijkp,input_p,mrat,m_arr
!         endif
!      endif
      
!      close(211)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! write para.h file

      if (my_rank .eq. 0) then


         open(109, file=trim(out_dir)//'para.dat',
     x        status='unknown',form='unformatted')
         
         write(109) nx,ny,nz,dx,dy,delz
         write(109) nt,dtsub_init,ntsub,dt,nout
         write(109) out_dir
!         write(109) model_choice
!         write(109) nf_init,b0_init
!         write(109) nu_init,lww2,lww1
!         write(109) Mdot,Mdot_part
         write(109) vtop,vbottom
         write(109) Ni_max
         write(109) mproton,m_pu,m_heavy
         write(109) np_top,np_bottom
         write(109) b0_top,b0_bottom
         write(109) vth_top,vth_bottom
!         write(109) RIo
         write(109) alpha,beta
!         write(109) comm_sz
         write(109) RIo
         close(109)

      endif
 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
! Initialize diagnostic output files
!----------------------------------------------------------------------


      open(110,file=trim(out_dir)//
     x     'c.np_'//filenum//'.dat',
     x     status='unknown',form='unformatted')
      open(111,file=trim(out_dir)//
     x     'c.np_3d_'//filenum//'.dat',
     x     status='unknown',form='unformatted')

      open(115,file=trim(out_dir)//
     x     'c.np_3d_1_'//filenum//'.dat',
     x     status='unknown',form='unformatted')
      open(116,file=trim(out_dir)//
     x     'c.np_3d_2_'//filenum//'.dat',
     x     status='unknown',form='unformatted')


!      open(115,file=trim(out_dir)//
!     x     'c.nf_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')
!      open(116,file=trim(out_dir)//
!     x     'c.nf_3d_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!      open(120,file=trim(out_dir)//
!     x     'c.uf_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')
!      open(121,file=trim(out_dir)//
!     x     'c.uf_3d_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

      open(130,file=trim(out_dir)//
     x     'c.b1_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')
      open(131,file=trim(out_dir)//
     x     'c.b1_3d_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(140,file=trim(out_dir)//
     x     'c.aj_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(150,file=trim(out_dir)//
     x     'c.E_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

!      open(160,file=trim(out_dir)//
!     x     'c.energy_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!      open(170,file=trim(out_dir)//
!     x     'c.chex_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!      open(172,file=trim(out_dir)//
!     x     'c.bill_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!      open(175,file=trim(out_dir)//
!     x     'c.satnp_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

      open(180,file=trim(out_dir)//
     x     'c.up_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')
      open(181,file=trim(out_dir)//
     x     'c.up_3d_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

!      open(190,file=trim(out_dir)//
!     x     'c.momentum_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!      open(192,file=trim(out_dir)//
!     x     'c.p_conserve_'//filenum//'.dat',
!     x     status='unknown',               
!     x     form='unformatted')                 

      open(300,file=trim(out_dir)//
     x     'c.temp_p_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')
      open(301,file=trim(out_dir)//
     x     'c.temp_p_3d_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(305,file=trim(out_dir)//
     x     'c.xp_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(310,file=trim(out_dir)//
     x     'c.vp_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(315,file=trim(out_dir)//
     x     'c.beta_p_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(320,file=trim(out_dir)//
     x     'c.mrat_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(320,file=trim(out_dir)//
     x     'c.mrat_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')


      open(330,file=trim(out_dir)//
     x     'c.temp_p_3d_1_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      open(331,file=trim(out_dir)//
     x     'c.temp_p_3d_2_'//filenum//'.dat',
     x     status='unknown',
     x     form='unformatted')

      

!      open(320,file='c.uf2.dat',status='unknown',
!     x         form='unformatted')

!      open(330,file='c.ufp2.dat',status='unknown',
!     x         form='unformatted')

!      open(340,file='c.eta.dat',status='unknown',
!     x         form='unformatted')

!      open(350,file=trim(out_dir)//
!     x     'c.pf_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')
!      open(351,file=trim(out_dir)//
!     x     'c.pf_3d_'//filenum//'.dat',
!     x     status='unknown',
!     x     form='unformatted')

!----------------------------------------------------------------------


!======================================================================
!  MAIN LOOP!
!======================================================================

      do 1 m = mstart+1, nt

         if (my_rank .eq. 0) then
            write(*,*) 'time...', m, m*dt
         endif

         !Calculate neutral density


         !Ionize cloud and calculate ion density
         write(*,*) 'Ni_tot...',Ni_tot,Ni_max,my_rank
         call separate_np(np_2,mr)
         if (Ni_tot .lt. 0.9*Ni_max) then
!          call Ionize_Io(np,vp,vp1,xp,xp1,up,ndot)
            mr = 1.0/m_pu
            call Ionize_pluto_mp(np,np_2,vp,vp1,xp,m,input_p,up)
         endif

         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
         call update_np_boundary(np)

         !energy diagnostics
         
         call get_bndry_Eflux(b1,E)
         call Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
     x                    EeP,nu,up,np)

         call curlB(b1,np,aj)
!         call obstacle_boundary_B(b0,b1)

!         call cov_to_contra(bt,btmf)
!         call face_to_center(btmf,btc)       !interp bt to cell center
         
         call edge_to_center(bt,btc)

         call extrapol_up(up,vp,vp1,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call improve_up(vp1,vplus,vminus,up,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call get_vp_final(Ep,vp,vp1,vplus)
         
         call move_ion_half(xp,vp,vp1,input_p,Ep)

                  !1/2 step ion move to n+1/2
!         call check_min_den_boundary(np,xp,vp,up)

         call get_Ep_buf(Ep_buf,b0,xp_buf,up)
         call get_vplus_vminus_buf(Ep_buf,vp_buf,vplus_buf,
     x        vminus_buf,b0)
         call get_vp_buf_final(Ep_buf,vp_buf,vplus_buf)
         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
c         call part_setup_buf(xp_buf,vp_buf)
c         call move_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
c     x        B_out_buf,mrat_out_buf)
         call exchange_ion_in(xp,vp,vp1,input_p,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,input_p,xp_buf,vp_buf,
     x        E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)
c         call exchange_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
c     x        B_out_buf,mrat_out_buf,xp,vp,vp1)

         call part_setup_buf(xp_buf,vp_buf)


         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
         ndiag = ndiag + 1
         if (ndiag .eq. nout) then         
            call get_temperature(xp,vp,np,temp_p)
!            mr = 1.0
!            call separate_np(np_1,mr)
!            mr = 1.0/m_pu
!            call separate_np(np_2,mr)
!            mr = 1.0
!            call separate_temp(vp,temp_p_1,mr)
!            mr = 1.0/m_pu
!            call separate_temp(vp,temp_p_2,mr)
         endif
         call update_np_boundary(np)

         
!**********************************************************************
! SUBCYCLING LOOP!
!**********************************************************************

         dtsub = dtsub_init
         ntf = ntsub
!all MPI_Barrier(MPI_COMM_WORLD,ierr)

         
         call check_time_step(bt,np)

         count = 1

         call MPI_ALLREDUCE(ntf,recvbuf,count,
     x        MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)

!         write(*,*) 'nft max...',recvbuf
         ntf = recvbuf

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 2 n = 1, int(ntf)

!         write(*,*) 'subcycle step...',n,ntf

         !convert main cell covarient bt to main cell contravarient
!         call cov_to_contra(bt,btmf) 
!         call edge_to_center(bt,btc)
         call curlB(b1,np,aj)     
!         call obstacle_boundary_B(b0,b1)

         !update fluid velocity, uf 

! only need predict_uf when calculating ugradu

!c         call trans_nf_Lax(nf,nf1,nfp1,uf) 
!         call trans_nf_LaxWend1(nf,nf1,nfp1,uf)
!         call trans_pf_LaxWend1(pf,pf1,uf)

!         call get_nuin(nuin,nn,uf)
!         call predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf1,uplus, 
!     x                   uminus,ugradu,up,gradP,nuin,bdp,pf1)

!         call predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)  

!         call get_nuin(nuin,nn,uf)
!         call correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
!     x                   ugradu,aj,up,ufp1,gradP,nuin,pf)

!         call trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
!         call trans_pf_LaxWend2(pf,pf1,ufp1)

         !update magnetic field, b1
!         call predict_B(b1,b12,b1p2,bt,btmf,E,aj,up,uf,uf2,np,nf,nu,
!     x                  gradP) 


         call predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu) 
!         call predict_B(b0,b1,b12,b1p2,bt,btmf,E,aj,up,np,nu) 


!         call correct_nf(nf,nf1,ufp1)

!         call correct_B(b0,b1,b1p2,E,aj,up,uf,np,nfp1,nu,gradP,bdp)
         call correct_B(b0,b1,b1p2,E,aj,up,np,nu)


!         call f_update_tlev(uf,uf2,b1,b12,b1p2,bt,b0,bdp)
         call f_update_tlev(b1,b12,b1p2,bt,b0)

!         call Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
!         call check_momentum(uf,nf,bt,ugradu)
!         write(192) m
!         write(192) n
!         write(192) surf_tot, graduu_tot, ugradu_tot



 2     continue
!**********************************************************************


         call move_ion_half(xp,vp,vp1,input_p,Ep)

         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
!         call part_setup_buf(xp_buf,vp_buf)

!         call move_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
!     x        B_out_buf,mrat_out_buf)
         
         call exchange_ion_in(xp,vp,vp1,input_p,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,input_p,xp_buf,vp_buf,
     x        E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)
!         call exchange_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
!     x        B_out_buf,mrat_out_buf,xp,vp,vp1)

         call part_setup_buf(xp_buf,vp_buf)
!         call check_min_den_boundary(np,xp,vp,up)

         call check_min_den(np,xp,vp,vp1,up,bt)

!         if (Ni_tot .lt. 0.9*Ni_max) then
!            call res_chex(xp,vp,vp1)
!         endif

!         endif

!         write(*,*) 'Momentum conservation...'
!         write(*,*) '  Particles.............',pup(1),pup(2),pup(3)
!         write(*,*) '  Fluid.................',puf(1),puf(2),puf(3)
!         write(*,*) '  ExB...................',peb(1),peb(2),peb(3)
!         write(*,*) '  Normalized............',
!     x                        (pup(1)+puf(1)+peb(1))/input_p(1),
!     x                        (pup(2)+puf(2)+peb(2))/input_p(2),
!     x                        (pup(3)+puf(3)+peb(3))/input_p(3)

!         call get_np3(np,np3)

!         call update_mixed


!----------------------------------------------------------------------
! diagnostic output
!----------------------------------------------------------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!         if (my_rank .eq. 0) then
!            write(160) m
!            write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
!     x           EeP,input_chex,input_bill
!            write(190) m
!            write(190) pup, puf, peb, input_p
!            write(320) np(ri-20,rj,rk),np(ri-40,rj,rk),
!     x                 np(ri-40,rj,rk+50),np(ri+5,rj,rk)
!         endif


         if (ndiag .eq. nout) then

!            call separate_np(np_1,1.0)
!            call separate_np(np_2,1/m_pu)
!            call get_temperature(xp,vp,np,temp_p)
!            call separate_temp(temp_p_1,1.0)
!            call separate_temp(temp_p_2,1/m_pu)

            nproc_2rio = nint(100*rio/(delz*nz))
!            write(*,*) 'nproc_2rio....',nproc_2rio,
!     x           (comm_sz/2)-nproc_2rio

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! save 3d arrays------------------------
            if ((my_rank .ge. comm_sz/2-nproc_2rio) .and. 
     x           (my_rank .lt. comm_sz/2+nproc_2rio)) then
               write(111) m
               write(111) np
!               write(115) m
!               write(115) np_1
!               write(116) m
!               write(116) np_2
               write(131) m
               write(131) bt
!               write(140) m
!               write(140) aj
!               write(150) m
!               write(150) E
               write(181) m
               write(181) up
               write(301) m
               write(301) temp_p/1.6e-19
               write(305) m
               write(305) xp
               write(310) m
               write(310) vp
               write(315) m
               write(315) beta_p
               write(320) m
               write(320) mrat
!               write(330) m 
!               write(330) temp_p_1/1.6e-19
!               write(331) m 
!               write(331) temp_p_2/1.6e-19

           endif


! save 2d arrays----------------------
               write(110) m
               write(110) np(:,ny/2,:),np(:,:,2)
               write(130) m
               write(130) bt(:,ny/2,:,:),bt(:,:,2,:)
!               write(140) m
!               write(140) aj
!               write(150) m
!               write(150) E
               write(180) m
               write(180) up(:,ny/2,:,:),up(:,:,2,:)
               write(300) m
               write(300) temp_p(:,ny/2,:)/1.6e-19,
     x                    temp_p(:,:,2)/1.6e-19
               ndiag = 0

         endif

!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Write restart file
!----------------------------------------------------------------------

!         if (my_rank .eq. 0) then
!            if (m .eq. mrestart) then
!               open(220,file='restart.vars.new',status='unknown',
!     x              form='unformatted')
!               write(220) b0,b1,b12,b1p2,bt,btmf,btc,np,np3,vp,vp1,
!     x              vplus,vminus,
!     x              up,xp,xp1,aj,nu,Ep,E,Evp,Euf,
!     x              EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
!     x              Ni_tot_sys,
!     x              ijkp,mrestart,input_p,input_EeP,prev_Etot,
!     x              mrat,m_arr,xp_buf,vp_buf,Ep_buf,vplus_buf,
!     x              vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
!     x              B_out_buf,mrat_out_buf,m_arr_out_buf,
!     x              in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
!     x              m_arr_buf,mrat_buf
!c                    write(220) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,
!c     x           vplus,vminus,
!c     x          up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
!c     x          EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
!c     x          ijkp,mrestart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1,
!c     x          input_chex,input_bill,pf,pf1,mrat,m_arr  
!            endif
!         endif
         
         if (my_rank .ge. 0) then
         if (m .eq. mrestart) then
            write(*,*) 'writing restart file....',
     x  'restart.part'//filenum//'.new'

          open(220,file='restart.vars'//filenum//'.new',
     x           status='unknown',
     x           form='unformatted')

          write(220)  b0,b1,b12,b1p2,bt,btc,np,np3,
     x         up,aj,nu,E,input_E,input_p,mrestart,input_EeP,prev_Etot,
     x         Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP

          open(221,file='restart.part'//filenum//'.new',
     x           status='unknown',form='unformatted')
          write(221) vp,vp1,vplus,vminus,
     x         xp,Ep,Ni_tot,
     x         Ni_tot_sys,ijkp,
     x         mrat,
     x         xp_buf,vp_buf,Ep_buf,vplus_buf,
     x         vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
     x         B_out_buf,mrat_out_buf,
     x         in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
     x         mrat_buf
            
!       write(211) b0,b1,b12,b1p2,bt,btmf,btc,np,np3
!     x         vp,vp1,vplus,vminus,xp,xp1,Ep,input_E,Ni_tot,
!     x         ijkp,input_p,mrat,m_arr,
!     x         xp_buf,vp_buf,Ep_buf,vplus_buf,
!     x         vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
!     x         B_out_buf,mrat_out_buf,m_arr_out_buf,
!     x         in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
!     x         m_arr_buf,mrat_buf
            
         endif
         endif

!----------------------------------------------------------------------


         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


 1     continue
!======================================================================

!       if(my_rank .eq. 0) then

          close(110)
          close(115)
          close(116)
          close(120)
          close(130)
          close(140)
          close(150)
          close(160)
          close(170)
          close(172)
          close(175)
          close(180)
          close(190)
          close(192)
          close(210)
          close(211)
          close(220)
          close(221)
          close(300)
          close(305)
          close(310)
          close(315)
          close(320)
          close(330)
          close(331)
!     close(340)
          close(350)

!       endif

          deallocate(filenum)

       call system_clock(t2,cnt_rt)
       time = (real(t2) - real(t1))/real(cnt_rt)
       if (my_rank .eq. 0) then
          write(*,*) 
          write(*,*) 'Elapsed time....',time,' sec'
          write(*,*)
       endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


       call MPI_FINALIZE(ierr)

       stop 
       end



















