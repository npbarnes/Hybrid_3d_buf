      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 100, ny = 75, nz = 10)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 1200000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
