      MODULE part_init


      USE global
      USE dimensions
      USE misc
      USE gutsp_dd
      USE mpi

      contains


!----------------------------------------------------------------------
      SUBROUTINE Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,
     x                       EE,EeP,nu,up,np)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real vp(Ni_max,3),
!     x     uf(nx,ny,nz,3),
!     x     nf(nx,ny,nz),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
!     x     etemp(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      real mO_q


      real Evp                  !kinetic energy of particles
      real Euf                  !kinetic energy of fluid flow
      real EB1,EB1x,EB1y,EB1z   !Magnetic field energy 
      real EE                   !Electric field energy 
      real EeP                  !Electron pressure energy
      real total_E              !total energy
      real aveEvp               !average particle energy
      real norm_E               !normalized energy
      real vol                  !volume of cell
      real denf                 !fluid density

      real recvbuf
      integer count
      count = 1

      mO_q = mion/q

      Euf = 0.0
      EB1 = 0.0
      EB1x = 0.0
      EB1y = 0.0
      EB1z = 0.0
      EE = 0.0
      EeP = 0.0


      do 10 i=1,nx-1
!         j = 2
         do 10 j=1,ny-1
            do 10 k=1,nz-1
               vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
               EB1x = EB1x + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,1))**2 
               EB1y = EB1y + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,2))**2 
               EB1z = EB1z + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,3))**2 
!               EeP = EeP + kboltz*etemp(i,j,k)
               do 10 m=1,3
                  denf = np(i,j,k)/(km_to_m**3)
                  Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
!                  EB1 = EB1 + 
!     x              (vol/(2.0*mu0))*(mO_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                  EB1 = EB1 + 
     x              (vol/(2.0*mu0))*(mO_q*b1(i,j,k,m))**2
                  EE = EE + (epsilon*vol/2.0)*
     x                      (mO_q*E(i,j,k,m)*km_to_m)**2
 10               continue

!      input_EeP = input_EeP + EeP

!      write(*,*) 'Energy diag...',Ni_tot,m_arr(2000000)
 
      Evp = 0.0
      do 15 l=1,Ni_tot
         do 15 m=1,3
            Evp = Evp + 0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 15   continue

!      write(*,*) 'Energy diag 2...',Ni_tot,m_arr(2000000)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(Evp,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_Evp = recvbuf

!      write(*,*) 'recvbuf...',recvbuf,Evp

      call MPI_ALLREDUCE(input_E,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_input_E = recvbuf

      call MPI_ALLREDUCE(bndry_Eflux,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_bndry_Eflux = recvbuf

!      total_E = S_Evp+EE+EB1
      total_E = S_Evp+EB1
      aveEvp = S_Evp/S_input_E
      
      if (my_rank .eq. 0) then

!      write(*,*) 'Input energy (J).............',S_input_E
!c      write(*,*) 'Input EeP energy (J).........',input_EeP
!      write(*,*) 'Total vp energy (J)..........',S_Evp
!      write(*,*) 'Total up energy (J)..........',Euf
!      write(*,*) 'Total B energy (J)...........',EB1/S_input_E
!      write(*,*) 'Total E energy (J)...........',EE/S_input_E
!c      write(*,*) 'Total EeP energy (J).........',EeP
!      write(*,*) 'Total energy (J).............',total_E
!c      write(*,*) 'Total energy w/ eP (J).......',total_E+EeP
!      write(*,*) 'Energy thru boundaries.......',bndry_Eflux/S_input_E
      write(*,*) 'Normalized particle energy...',aveEvp
      write(*,*) 'Normalized energy............',total_E/S_input_E,
     x   my_rank
      write(*,*) 'Normalized energy (bndry)....',
!     x                S_bndry_Eflux/total_E
     x                (total_E)/(S_input_E+S_bndry_Eflux)
!      write(*,*) 'Normalized energy (no b1z)...',(S_Evp+Euf+EE+EB1x+
!     x                                            EB1y)/S_input_E
!c      write(*,*) 'Normalized energy (w/ eP)....',
!c     x                             (total_E+EeP)/(input_E + input_EeP)
!      write(*,*) ' '

      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      norm_E = total_E/S_input_E

!      if (prev_Etot .eq. 0.0) then prev_Etot = norm_E
!      do 20 i=1,nx 
!         do 20 j=1,ny
!            do 20 k=1,nz
!               nu(i,j,k) = nu(i,j,k) + 
!     x                 nu(i,j,k)*2.0*((norm_E - prev_Etot)/norm_E)
! 20            continue
      prev_Etot = norm_E

      return
      end SUBROUTINE Energy_diag
!----------------------------------------------------------------------



!c----------------------------------------------------------------------
!      SUBROUTINE sw_part_setup(np,vp,vp1,xp,xp1,input_p,up)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real vp(Ni_max,3)
!      real vp1(Ni_max,3)
!      real xp(Ni_max,3)
!      real xp1(Ni_max,3)
!      real input_p(3)
!      real up(nx,ny,nz,3)

!c      Ni_tot = 120000

!      do 10 l = 1,Ni_tot

!         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
!         xp(l,2) = qy((ny/2)-8)+(1.0-pad_ranf())*(qy((ny/2)+8)-
!     x                                        qy((ny/2)-8))
!         xp(l,3) = qz((nz/2)-8)+(1.0-pad_ranf())*(qz((nz/2)+8)-
!     x                                        qz((nz/2)-8))

!c         i = nint(nx*pad_ranf())
!c         j = nint((ny/2) + 16.*(0.5-pad_ranf()))
!c         k = nint((nz/2) + 16.*(0.5-pad_ranf()))
!cc         write(*,*) 'l...',l,i,j,k

!c         xp(l,1) = qx(i)+dx*(0.5-pad_ranf())
!c         xp(l,2) = qy(j)+dy*(0.5-pad_ranf())
!c         xp(l,3) = qz(k)+dz_grid(k)*(0.5-pad_ranf())

!         vp(l,1) = -vsw
!         vp(l,2) = 0.0
!         vp(l,3) = 0.0


!         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!         ijkp(l,2) = nint(xp(l,2)/dy)
         
!         k=1
!         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!            ijkp(l,3) = k       !grid
!            k=k+1
! 50      continue
!         k=ijkp(l,3)
!         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!            ijkp(l,3) = k+1
!         endif
!         do 45 m=1,3
!            vp1(l,m) = vp(l,m)
!            input_E = input_E + 
!     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
!            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
! 45      continue
! 10      continue

        
 
!      write(*,*) 'get interp weights...'
!      call get_interp_weights(xp,xp1)
!      write(*,*) 'update_np...'
!      call update_np(np)
!      write(*,*) 'update_up...'
!      call update_up(vp,np,up)
!      write(*,*) 'update_up complete...'

!      return
!      end
!c----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE sw_part_setup_temp(np,vp,vp1,xp,xp1,input_p,up)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real vp(Ni_max,3)
!      real vp1(Ni_max,3)
!      real xp(Ni_max,3)
!      real xp1(Ni_max,3)
!      real input_p(3)
!      real up(nx,ny,nz,3)
!      real phi,theta,rnd,f,v
!      real rand

!c      Ni_tot = 120000

!c      do n = 0,procnum-1
!c      if (my_rank .eq. n) then
!      do 10 l = 1,Ni_tot
!c         write(*,*) 'procnum, random number...',n,pad_ranf()

!         phi = 2.0*pi*pad_ranf()
!         flg = 0
!         do 30 while (flg .eq. 0)
!            theta = pi*pad_ranf()
!            f = sin(theta)
!            rnd = pad_ranf()
!            if (f .ge. rnd) flg = 1
! 30      continue


!         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
!         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
!         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))


!         flg = 0
!         do 40 while (flg .eq. 0)
!            v = (100*pad_ranf())
!c            f = (vth**2/exp(1.0))*v**2*exp(-(v)**2 / vth**2)
!            f = exp(-(v)**2 / vth**2)
!            rnd = pad_ranf()
!            if (f .ge. rnd) then 
!               flg = 1
!               vp(l,1) = vsw + v*cos(phi)*sin(theta)
!               vp(l,2) = v*sin(phi)*sin(theta)
!               vp(l,3) = v*cos(theta)
!            endif

!c         vp(l,1) = -vsw
!c         vp(l,2) = 0.0
!c         vp(l,3) = 0.0

! 40      continue


!         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!         ijkp(l,2) = nint(xp(l,2)/dy)
         
!         k=1
!         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!            ijkp(l,3) = k       !grid
!            k=k+1
! 50      continue
!         k=ijkp(l,3)
!         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!            ijkp(l,3) = k+1
!         endif
!         do 45 m=1,3
!            vp1(l,m) = vp(l,m)
!            input_E = input_E + 
!     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
!            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
! 45      continue
! 10      continue
!c      endif
!c      enddo

 
!      write(*,*) 'get interp weights...'
!      call get_interp_weights(xp,xp1)
!      write(*,*) 'update_np...'
!      call update_np(np)
!      write(*,*) 'update_up...'
!      call update_up(vp,np,up)
!      write(*,*) 'update_up complete...'

   

!      return
!      end
!c----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE maxwl_init(vthm,vx,vy,vz)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real vthm
      integer flg
      real v,vx,vy,vz
      real rnd,f

      vx = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())
      vy = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())
      vz = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())


!      flg = 0
!      do 40 while (flg .eq. 0)
!         v = (2*vth_max*pad_ranf())-vth_max
!         f = exp(-(v)**2 / vthm**2)
!         rnd = pad_ranf()
!         if (f .ge. rnd) then 
!            flg = 1
!            vx = v
!         endif
! 40   continue
!      flg = 0
!      do 42 while (flg .eq. 0)
!         v = (2*vth_max*pad_ranf())-vth_max
!         f = exp(-(v)**2 / vthm**2)
!         rnd = pad_ranf()
!         if (f .ge. rnd) then 
!            flg = 1
!            vy = v
!         endif
! 42   continue
!      flg = 0
!      do 44 while (flg .eq. 0)
!         v = (2*vth_max*pad_ranf())-vth_max
!         f = exp(-(v)**2 / vthm**2)
!         rnd = pad_ranf() 
!         if (f .ge. rnd) then 
!            flg = 1
!            vz = v
!         endif
! 44   continue

      return
      end SUBROUTINE maxwl_init
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up)
!----------------------------------------------------------------------
!      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
!      integer np_t_flg(Ni_max)
!      integer np_b_flg(Ni_max)

      integer flg
      real nprat
      real bwght
      real mpart

      v1 = 1.0

!      np_t_flg(:) = 0
!      np_b_flg(:) = 0


      nprat = np_bottom/np_top

!      Ni_tot_O = Ni_tot*(2./3.)
!      Ni_tot_S = Ni_tot*(1./3.)

      Ni_tot = Ni_tot + Ni_tot*nint((1/nprat)-1)/2

! initialize protons


      bwght = 1.0
      vth = vth_top
      mpart = mproton

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         

         i=0
 11      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 11 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


!         ijkp(l,2) = floor(xp(l,2)/dy) 

         j=0
 13      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 13 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j

         
         k=2
 12      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 12 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k


!         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!         ijkp(l,2) = nint(xp(l,2)/dy)
         
!         k=1
!         do 12 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!            ijkp(l,3) = k       !grid
!            k=k+1
! 12      continue
!         k=ijkp(l,3)
!         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!            ijkp(l,3) = k+1
!         endif

         call maxwl_init(vth,vx,vy,vz)

         ii = ijkp(l,1)
         kk = ijkp(l,3)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz 

!         m_arr(l) = mpart
         mrat(l) = mproton/mpart
         beta_p(l) = bwght

         do 20 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*bwght)
            input_p(m) = input_p(m) + mpart*vp(l,m) / 
     x           (beta*bwght)
 20      continue
                 
 10      continue

! initialize He++

      Ni_tot_1 = Ni_tot + 1
      
      Ni_tot = Ni_tot + f_mq_2*Ni_tot

      do 30 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


!         ijkp(l,2) = floor(xp(l,2)/dy) 

         j=0
 33      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 33 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j


         k=2
 32      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 32 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k

!         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!         ijkp(l,2) = nint(xp(l,2)/dy)
         
!         k=1
!         do 32 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!            ijkp(l,3) = k       !grid
!            k=k+1
! 32      continue
!         k=ijkp(l,3)
!         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!            ijkp(l,3) = k+1
!         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         call maxwl_init(vth,vx,vy,vz)

!         ii = ijkp(l,1)
!         kk = ijkp(l,3)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz 

!         m_arr(l) = 2*mproton
         mrat(l) = 1.0/2.0
         beta_p(l) = b_mq_2

         do 40 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
            input_p(m) = input_p(m) + (mion/mrat(l))*vp(l,m) / 
     x           (beta*beta_p(l))
 40      continue
 30   continue


! add shell distribution

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + f_shl*Ni_tot

         do 69 l = Ni_tot_1,Ni_tot


            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
            xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
            
!            m_arr(l) = mproton
            mrat(l) = 1.0
            beta_p(l) = b_shl


            i=0
 71         continue
            i = i + 1
            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
            i = i-1
            ijkp(l,1)= i
            
            
!            ijkp(l,2) = floor(xp(l,2)/dy) 


            j=0
 73         continue
            j = j + 1
            if (xp(l,2) .gt. qy(j)) go to 73 !find i on non-uniform 
            j = j-1
            ijkp(l,2)= j

            
            k=2
 72         continue
            k = k + 1
            if (xp(l,3) .gt. qz(k)) go to 72 !find k on non-uniform 
            k = k-1
            ijkp(l,3)= k

!            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!            ijkp(l,2) = nint(xp(l,2)/dy)
            
!            k=1
!            do 70 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!               ijkp(l,3) = k    !grid
!               k=k+1
! 70         continue
!            k=ijkp(l,3)
!            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!               ijkp(l,3) = k+1
!            endif
            
!            ii = ijkp(l,1)
!            kk = ijkp(l,3)

            theta = pad_ranf()*PI
            phi = pad_ranf()*2*PI
!            vp(l,1) = 1.0*vsw*cos(theta) + vx !+dvx
!            vp(l,2) = vy 
!            vp(l,3) = 1.0*vsw*sin(theta) + vz        !+dvz 

            vp(l,1) = -vsw+vsw*cos(phi)*sin(theta) !+dvx
            vp(l,2) = vsw*sin(phi)*sin(theta) !+dvz 
            vp(l,3) = vsw*cos(theta)

!            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
!            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              (beta*b_shl)
               input_p(m) = input_p(m)+(mion/mrat(l))*vp(l,m)/
     x              (beta*b_shl)
            enddo
            

 69      enddo



      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl
!----------------------------------------------------------------------


!c----------------------------------------------------------------------
!      SUBROUTINE sw_part_setup_maxwl_1(np,vp,vp1,xp,xp1,input_p,up,
!     x     np_t_flg,np_b_flg)
!c----------------------------------------------------------------------
!      include 'incurv.h'

!      real np(nx,ny,nz)
!      real vp(Ni_max,3)
!      real vp1(Ni_max,3)
!      real xp(Ni_max,3)
!      real xp1(Ni_max,3)
!      real input_p(3)
!      real up(nx,ny,nz,3)
!      real phi,theta,rnd,f,v
!      real rand
!      real vx,vy,vz
!      real dvx,dvz,v1
!      integer np_t_flg(Ni_max)
!      integer np_b_flg(Ni_max)

!      integer flg
!      real nprat

!      v1 = 1.0

!      np_t_flg(:) = 0
!      np_b_flg(:) = 0


!      nprat = np_bottom/np_top

!c      Ni_tot = 120000

!c      do n = 0,procnum-1
!c      if (my_rank .eq. n) then
!      Ni_tot = Ni_tot + Ni_tot*nint((1/nprat)-1)/2

!      do 10 l = 1,Ni_tot
!c         write(*,*) 'procnum, random number...',n,pad_ranf()

!c         phi = 2.0*pi*pad_ranf()
!c         flg = 0
!c         do 30 while (flg .eq. 0)
!c            theta = pi*pad_ranf()
!c            f = sin(theta)
!c            rnd = pad_ranf()
!c            if (f .ge. rnd) flg = 1
!c 30      continue


!         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
!         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
!         flg = 0
!         do 20 while (flg .eq. 0)
!            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
!             rnd = (1-nprat)* 
!     x             ((1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) +
!     x             nprat
!c            write(*,*) 'rnd...',rnd
!             if (pad_ranf() .le. rnd) flg = 1
!             if (xp(l,3) .ge. qz(nz/2)) np_t_flg(l) = 1
!             if (xp(l,3) .lt. qz(nz/2)) np_b_flg(l) = 1

! 20      continue

!         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
!         ijkp(l,2) = nint(xp(l,2)/dy)
         
!         k=1
!         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
!            ijkp(l,3) = k       !grid
!            k=k+1
! 50      continue
!         k=ijkp(l,3)
!         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
!            ijkp(l,3) = k+1
!         endif

!         vth = 0.5*(vth_top + vth_bottom) + 
!     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

!         flg = 0
!         do 40 while (flg .eq. 0)

!            vy = (400*pad_ranf())-200
!            vz = (400*pad_ranf())-200
!            vx = (400*pad_ranf())-200
            
!            v = sqrt(vx**2 + vy**2 + vz**2)
!            f = exp(-(v)**2 / vth**2)
!            rnd = pad_ranf() 
!            if (f .gt. rnd) then 
!               flg = 1
!            endif
! 40      continue

!         ii = ijkp(l,1)
!         kk = ijkp(l,3)
!         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
!     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
!         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
!     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
!         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
!         vp(l,2) = vy 
!         vp(l,3) = vz !+dvz 

!         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
!         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
!         do 45 m=1,3
!            vp1(l,m) = vp(l,m)
!            input_E = input_E + 
!     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
!            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
! 45      continue
! 10      continue
!c      endif
!c      enddo



 
!      write(*,*) 'get interp weights...'
!      call get_interp_weights(xp,xp1)
!      write(*,*) 'update_np...'
!      call update_np(np)
!      write(*,*) 'update_up...'
!      call update_up(vp,np,up)
!      write(*,*) 'update_up complete...'

   

!      return
!      end
!c----------------------------------------------------------------------


!----------------------------------------------------------------------
      SUBROUTINE load_Maxwellian(np,vp,vp1,xp,input_p,up,vth,Ni_tot_1)
!----------------------------------------------------------------------

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
!      real phi,theta,rnd,f,v
!      real rand
!      real vx,vy,vz
!      real dvx,dvz,v1
!      integer np_t_flg(Ni_max)
!      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real vth

!      integer flg
!      real nprat

      v1 = 1.0

!      np_t_flg(:) = 0
!      np_b_flg(:) = 0

      do 10 l = 1,Ni_tot_1

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
!         m_arr(l) = mion
         mrat(l) = 1.0

!         ijkp(l,1) = floor(xp(l,1)/dx) 

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


!         ijkp(l,2) = floor(xp(l,2)/dy) 


         j=0
 33      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 33 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j
        
         k=0
 30      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 30 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
         
!         vth = vth_bottom

         vx = vsw+vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
         vy = vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
         vz = vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

!         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
!         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + (mion/mrat(l))*vp(l,m) / beta
 45      continue
 10      continue

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)

      return
      end SUBROUTINE load_Maxwellian
!----------------------------------------------------------------------





      end MODULE part_init









