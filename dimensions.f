      MODULE DIMENSIONS

! simulation domain dimensions

      PARAMETER (nx = 79, ny = 39, nz = 7)
! particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 2000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
