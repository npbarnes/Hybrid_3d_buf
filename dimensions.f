      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 300, ny = 200, nz = 20)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 20000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
