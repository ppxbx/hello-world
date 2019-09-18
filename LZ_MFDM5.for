      module dish_type
        type dish
          integer         row
          double complex  value
        end type
      end module dish_type

      use dish_type

      integer, parameter :: n1 = 801

      real*8,  parameter :: delta = 0.0d0,    pi    = 3.1415926d0,
     &                      k0    = 2.0d-1,   natom = 1.0d3

      double complex c1, c2, c3

      double complex za(n1), za1(n1), zb1(n1), za_temp(n1),
     &               zb(n1), za2(n1), zb2(n1), zb_temp(n1),
     &               zc(n1), za3(n1), zb3(n1),
     &                       za4(n1), zb4(n1)

      real*8 b1, dt, momentum,  coherence, fragmentation, time,
     &       b2, dk, momentum1, momentum3, omega, g11, g22, g12,
     &       b3,     momentum2

      real*8 trap1(n1), zd1(n1), zd(n1),
     &       trap2(n1), zd2(n1)

      integer kvp, kvp1, int1, step, kvp4, lzr, position,
     &        int, kvp2, int2,
     &        num, kvp3

      character(len=80) filename, filepath

      character(len=4) value

      type(dish), allocatable, dimension(:) :: zaa, zbb


      g11   = 2.0d-5
      g22   = 2.0d-5
      g12   = 1.0d-5
      dk    = 2.0d-2

      position = 1   !@@@@@ position refers to how many arguments
                     !are taken from command line
      call getarg(position, value)

      read(value, '(i4)') step

      if      (step.eq.1) then
         omega = 2.0d-1
      else if (step.eq.2) then
         omega = 4.0d-1
      else if (step.eq.3) then
         omega = 6.0d-1
      else if (step.eq.4) then
         omega = 8.0d-1
      else if (step.eq.5) then
         omega = 1.0d0
      else if (step.eq.6) then
         omega = 1.2d0
      else if (step.eq.7) then
         omega = 1.4d0
      else if (step.eq.8) then
         omega = 1.6d0
      else
      end if


      filepath ='/state/partition1/ppxbx/data'//trim(value)//'/'

      trap1 = 0.0d0
      trap2 = 0.0d0
      do i = 1, n1
         momentum = dble(i-1)*dk - 8.0d0
         trap1(i) = (momentum + 1.0d0)**2 - delta/2.0d0
         trap2(i) = (momentum - 1.0d0)**2 + delta/2.0d0
      end do

      open(1, file = trim(filepath)//'trap.dat',
     &    position = 'append')
      do i = 1, n1
         write(1, "(i5, 2f16.9)") i, trap1(i), trap2(i)
      end do
      close(1)

      za = (0.0d0, 0.0d0)
      zb = (0.0d0, 0.0d0)   !@@@@@ the initial state is a pure spin state
      do i = 1, n1
         momentum = dble(i-1)*dk - 8.0d0
         zb(i)    = exp(-(momentum + 1.0d0)**2/(2.0d0*k0**2))
     &            * (1.0, 0.0)
      end do
      call normalization(zb, n1, natom)

      time = 0.0d0
      kvp1 = 0
      kvp2 = 0
      kvp3 = 0
      do kvp = 0, 40000

         dt   = 1.0d-3
         time = time + dt


!!!!!    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ part one @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         za1 = (0.0d0, 0.0d0)
         do i = 1, n1
            za1(i) = za1(i) - (0.0, 1.0)*dt*omega*zb(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(za, zb, dk, dt, g11, g12, n1, zc)
         za1 = za1 + zc

         do i = 1, n1
            if (i.eq.1) then
               za1(i) = za1(i) - (0.0,1.0)*dt*trap1(i)*za(i) + (0.0,1.0)
     &                         * dt*k0**4*(za(i+1) - 2.0d0*za(i))/dk**2
            else if (i.eq.n1) then
               za1(i) = za1(i) - (0.0,1.0)*dt*trap1(i)*za(i) + (0.0,1.0)
     &                         * dt*k0**4*(za(i-1) - 2.0d0*za(i))/dk**2
            else
               za1(i) = za1(i) - (0.0,1.0)*dt*trap1(i)*za(i) + (0.0,1.0)
     &                         * dt*k0**4*(za(i-1) - 2.0d0*za(i)
     &                         + za(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step one : k1 = g[f(t_{i}), t_{i}] \delta t !!!!!!!!!!!!!

         za_temp = (0.0d0, 0.0d0)
         za_temp = za + 5.0d-1*za1

         za2 = (0.0d0, 0.0d0)
         do i = 1, n1
            za2(i) = za2(i) - (0.0, 1.0)*dt*omega*zb(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(za_temp, zb, dk, dt, g11, g12, n1, zc)
         za2 = za2 + zc

         do i = 1, n1
            if (i.eq.1) then
               za2(i) = za2(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i+1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else if (i.eq.n1) then
               za2(i) = za2(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else
               za2(i) = za2(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i)+za_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step two: k2 = g[f(t_{i})+k1/2, t_{i}/2] \delta t !!!!!!!!!!!!!

         za_temp = (0.0d0, 0.0d0)
         za_temp = za + 5.0d-1*za2

         za3 = (0.0d0, 0.0d0)
         do i = 1, n1
            za3(i) = za3(i) - (0.0, 1.0)*dt*omega*zb(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(za_temp, zb, dk, dt, g11, g12, n1, zc)
         za3 = za3 + zc

         do i = 1, n1
            if (i.eq.1) then
               za3(i) = za3(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i+1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else if (i.eq.n1) then
               za3(i) = za3(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else
               za3(i) = za3(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i)+za_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step three: k3 = g[f(t_{i})+k2/2, t_{i}/2] \delta t !!!!!!!!!!!!!

         za_temp = (0.0d0, 0.0d0)
         za_temp = za + za3

         za4 = (0.0d0, 0.0d0)
         do i = 1, n1
            za4(i) = za4(i) - (0.0, 1.0)*dt*omega*zb(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(za_temp, zb, dk, dt, g11, g12, n1, zc)
         za4 = za4 + zc

         do i = 1, n1
            if (i.eq.1) then
               za4(i) = za4(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i+1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else if (i.eq.n1) then
               za4(i) = za4(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i))/dk**2
            else
               za4(i) = za4(i) - (0.0, 1.0)*dt*trap1(i)*za_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(za_temp(i-1)
     &                         - 2.0d0*za_temp(i)+za_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step four: k4 = g[f(t_{i})+k3, t_{i}/2] \delta t !!!!!!!!!!!!!

         za_temp = (0.0d0, 0.0d0)
         za_temp = za + (za1 + 2.0d0*(za2 + za3) + za4)/6.0d0
         za      = (0.0d0, 0.0d0)
         za      = za_temp


!!!!!    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ part two @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         zb1 = (0.0d0, 0.0d0)
         do i = 1, n1
            zb1(i) = zb1(i) - (0.0, 1.0)*dt*omega*za(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(zb, za, dk, dt, g22, g12, n1, zc)
         zb1 = zb1 + zc

         do i = 1, n1
            if (i.eq.1) then
               zb1(i) = zb1(i) - (0.0,1.0)*dt*trap2(i)*zb(i) + (0.0,1.0)
     &                         * dt*k0**4*(zb(i+1) - 2.0d0*zb(i))/dk**2
            else if (i.eq.n1) then
               zb1(i) = zb1(i) - (0.0,1.0)*dt*trap2(i)*zb(i) + (0.0,1.0)
     &                         * dt*k0**4*(zb(i-1) - 2.0d0*zb(i))/dk**2
            else
               zb1(i) = zb1(i) - (0.0,1.0)*dt*trap2(i)*zb(i) + (0.0,1.0)
     &                         * dt*k0**4*(zb(i-1) - 2.0d0*zb(i)
     &                         + zb(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step one : k1 = g[f(t_{i}), t_{i}] \delta t !!!!!!!!!!!!!

         zb_temp = (0.0d0, 0.0d0)
         zb_temp = zb + 5.0d-1*zb1

         zb2 = (0.0d0, 0.0d0)
         do i = 1, n1
            zb2(i) = zb2(i) - (0.0, 1.0)*dt*omega*za(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(zb_temp, za, dk, dt, g22, g12, n1, zc)
         zb2 = zb2 + zc

         do i = 1, n1
            if (i.eq.1) then
               zb2(i) = zb2(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i+1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else if (i.eq.n1) then
               zb2(i) = zb2(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else
               zb2(i) = zb2(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i)+zb_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step two: k2 = g[f(t_{i})+k1/2, t_{i}/2] \delta t !!!!!!!!!!!!!

         zb_temp = (0.0d0, 0.0d0)
         zb_temp = zb + 5.0d-1*zb2

         zb3 = (0.0d0, 0.0d0)
         do i = 1, n1
            zb3(i) = zb3(i) - (0.0, 1.0)*dt*omega*za(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(zb_temp, za, dk, dt, g22, g12, n1, zc)
         zb3 = zb3 + zc

         do i = 1, n1
            if (i.eq.1) then
               zb3(i) = zb3(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i+1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else if (i.eq.n1) then
               zb3(i) = zb3(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else
               zb3(i) = zb3(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i)+zb_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step three: k3 = g[f(t_{i})+k2/2, t_{i}/2] \delta t !!!!!!!!!!!!!

         zb_temp = (0.0d0, 0.0d0)
         zb_temp = zb + zb3

         zb4 = (0.0d0, 0.0d0)
         do i = 1, n1
            zb4(i) = zb4(i) - (0.0, 1.0)*dt*omega*za(i)/2.0d0
         end do

         zc = (0.0d0, 0.0d0)
         call interaction(zb_temp, za, dk, dt, g22, g12, n1, zc)
         zb4 = zb4 + zc

         do i = 1, n1
            if (i.eq.1) then
               zb4(i) = zb4(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i+1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else if (i.eq.n1) then
               zb4(i) = zb4(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i))/dk**2
            else
               zb4(i) = zb4(i) - (0.0, 1.0)*dt*trap2(i)*zb_temp(i)
     &                         + (0.0, 1.0)*dt*k0**4*(zb_temp(i-1)
     &                         - 2.0d0*zb_temp(i)+zb_temp(i+1))/dk**2
            end if
         end do
!!!!!!!!!!!!!!! step four: k4 = g[f(t_{i})+k3, t_{i}/2] \delta t !!!!!!!!!!!!!

         zb_temp = (0.0d0, 0.0d0)
         zb_temp = zb + (zb1 + 2.0d0*(zb2 + zb3) + zb4)/6.0d0
         zb      = (0.0d0, 0.0d0)
         zb      = zb_temp

         b1 = 0.0d0
         do i = 1, n1
            b1 = b1 + conjg(za(i))*za(i)
         end do

         b2 = 0.0d0
         do i = 1, n1
            b2 = b2 + conjg(zb(i))*zb(i)
         end do

C         b3 = 0.0d0
C        b4 = 0.0d0
C        do i = 1, n1
C            momentum = dble(i-1)*dk - 8.0d0
C           if (b2.lt.1.0d-4) then
C              b4 = 0.0d0
C           else
C               b4 = b4 + momentum*conjg(zb(i))*zb(i)/b2
C         end if
C           if (b1.lt.1.0d-4) then
C          b3 = 0.0d0
C           else
C          b3 = b3 + momentum*conjg(za(i))*za(i)/b1
C            end if
C        end do

C         if (abs(b4).lt.1.0d-5) then
C           open(3, file = trim(filepath)//'LZR.dat',
C     &          position = 'append')
C         zd  = 0.0d0
C           lzr = 0.0d0
C           call phase(zb, n1, zd)
C           do i = 3, n1
C              lzr = lzr + (zd(i) - zd(i-2))/(2.0d0*dk)
C            end do
C           b1  = lzr
C           lzr = exp(-pi*omega**2/(4*abs(lzr)))
C           write(3, "(3f16.9)") b4, b1, lzr
C           close(3)
C        else
C        end if

         kvp2 = kvp1*400
         if (kvp.eq.kvp2) then
            kvp1 = kvp1 + 1

            zd1 = 0.0d0
            zd2 = 0.0d0
            zd1 = conjg(za)*za
            zd2 = conjg(zb)*zb

          if (kvp1.lt.10) then
               write(filename, '(i1)') kvp1
            else if (kvp1.ge.10.and.kvp1.lt.100) then
               write(filename, '(i2)') kvp1
            else
               write(filename, '(i3)') kvp1
            end if

           open(1, file = trim(filepath)//'den'//trim(filename)//'.dat',
     &         position = 'append')
           do i = 1, n1
              write(1, "(1x, i5, 2f16.9)") i, zd1(i), zd2(i)
           end do
           close(1)

            open(2, file = trim(filepath)//'quality.dat',
     &          position = 'append')
            b1 = 0.0d0
            do i = 1, n1
               b1 = b1 + conjg(za(i))*za(i)
            end do

            b2 = 0.0d0
            do i = 1, n1
               b2 = b2 + conjg(zb(i))*zb(i)
            end do

            b3 = 0.0d0
            b4 = 0.0d0
            do i = 1, n1
               momentum = dble(i-1)*dk - 8.0d0
               if (b2.lt.1.0d-4) then
                  b4 = 0.0d0
               else
                  b4 = b4 + momentum*conjg(zb(i))*zb(i)/b2
             end if
           if (b1.lt.1.0d-4) then
              b3 = 0.0d0
               else
              b3 = b3 + momentum*conjg(za(i))*za(i)/b1
               end if
            end do

            write(2, "(5f16.9)") time, b1, b2, b3, b4
            close(2)


         end if


      end do

      stop

      end


      subroutine interaction(za, zb, dk, dt, g11, g12, n1, zc)

      use            dish_type
      double complex za(n1), zb(n1), zc(n1)
      integer        n1,     kvp1,   kvp2
      double complex c1
      real*8         momentum, momentum1, momentum2, momentum3, g11,
     &               g12, dk,  dt

      type(dish), allocatable, dimension(:) :: zaa, zbb

      kvp1 = 0
      kvp2 = 0
      do i = 1, n1
         b1 = abs(za(i))
         b2 = abs(zb(i))
         if (b1.gt.3.0d-1) then
            kvp1 = kvp1 + 1
         end if
         if (b2.gt.3.0d-1) then
            kvp2 = kvp2 + 1
         end if
      end do
      allocate(zaa(kvp1))
      allocate(zbb(kvp2))
      zaa  = dish(0, (0.0, 0.0))
      zbb  = dish(0, (0.0, 0.0))
      kvp1 = 0
      kvp2 = 0
      do i = 1, n1
         b1 = abs(za(i))
         b2 = abs(zb(i))
         if (b1.gt.3.0d-1) then
            kvp1      = kvp1 + 1
            zaa(kvp1) = dish(i, za(i))
         end if
         if (b2.gt.3.0d-1) then
            kvp2      = kvp2 + 1
            zbb(kvp2) = dish(i, zb(i))
         end if
      end do

      momentum  = 0.0d0
      momentum1 = 0.0d0
      momentum2 = 0.0d0
      momentum3 = 0.0d0
      zc        = (0.0d0, 0.0d0)
      do i = 1, kvp1
         momentum = dble(zaa(i)%row-1)*dk - 8.0d0 !@@@ this term is from physics
         c1       = (0.0d0, 0.0d0)
         do j = 1, kvp1
            momentum1 = dble(zaa(j)%row-1)*dk - 8.0d0
            do k = 1, kvp1
               momentum2 = dble(zaa(k)%row-1)*dk - 8.0d0
               do l = 1, kvp1
                  momentum3 = dble(zaa(l)%row-1)*dk - 8.0d0
            b1        = momentum  + momentum1
            b2        = momentum2 + momentum3
                  if (abs(b1-b2).lt.1.0d-4) then
               c1     = c1 - (0.0, 1.0)*dt*g11
     &                           * conjg(zaa(j)%value)
     &                           * zaa(k)%value*zaa(l)%value
                  end if
             end do
          end do
         end do
         zc(zaa(i)%row) = zc(zaa(i)%row) + c1
      end do

      momentum  = 0.0d0
      momentum1 = 0.0d0
      momentum2 = 0.0d0
      momentum3 = 0.0d0
      do i = 1, kvp1
         momentum = dble(zaa(i)%row-1)*dk - 8.0d0
         c1       = (0.0d0, 0.0d0)
         do j = 1, kvp2
            momentum1 = dble(zbb(j)%row-1)*dk - 8.0d0
            do k = 1, kvp2
               momentum2 = dble(zbb(k)%row-1)*dk - 8.0d0
               do l = 1, kvp1
                  momentum3 = dble(zaa(l)%row-1)*dk - 8.0d0
            b1        = momentum  + momentum1
            b2        = momentum2 + momentum3
                  if (abs(b1-b2).lt.1.0d-4) then
               c1     = c1 - (0.0, 1.0)*dt*g12
     &                           * conjg(zbb(j)%value)
     &                           * zbb(k)%value*zaa(l)%value
                  end if
             end do
          end do
         end do
         zc(zaa(i)%row) = zc(zaa(i)%row) + c1
      end do

      deallocate(zaa)
      deallocate(zbb)

      return

      end



      subroutine normalization(za, n1, N)
      integer         n1
      real*8          N
      double complex  za(n1)
      b1 = 0.0d0
      do i = 1, n1
         b1 = b1 + conjg(za(i))*za(i)
      end do
      if (abs(b1).lt.1.0d-6) then
         za = (0.0d0, 0.0d0)
      else
         b1 = sqrt(dble(N)/b1)
         za = za*b1
      end if

      return

      end


      subroutine phase(za, n1, zd)
      real*8,  parameter ::  pi = 3.1415926d0
      double complex za(n1)
      real*8  zd(n1)
      integer n1
      real*8 b1, b2, b3

      b1    = 0.0d0
      b2    = 0.0d0
      b3    = 0.0d0
      theta = 0.0d0
      zd    = 0.0d0
      do i = 1, n1
         b1    = real(za(i))
         b2    = aimag(za(i))
         b3    = b2/b1
         if ((b1>0.0d0).and.(b2>0.0d0)) then
            zd(i) = atan(b3)
         else if((b1<0.0d0).and.(b2/=0.0d0)) then
            zd(i) = pi+atan(b3)
         else if((b1>0.0d0).and.(b2<0.0d0))  then
            zd(i) = 2*pi+atan(b3)
         else if((b1==0.0d0).and.(b2>0.0d0)) then
            zd(i) = 0.5*pi
         else if((b1==0.0d0).and.(b2<0.0d0)) then
            zd(i) = 1.5*pi
         else if((b1>0.0d0).and.(b2==0.0d0)) then
            zd(i) = 0.0d0
         else if((b1<0.0d0).and.(b2==0.0d0)) then
            zd(i) = pi
         else if((b1==0.0d0).and.(b2==0.0d0)) then
            zd(i) = 0.0d0
         else
         end if
      end do

      return

      end


