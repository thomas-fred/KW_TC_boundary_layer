module cal_slp_mod
   contains

   ! Ming Ge Feb. 2 2021 calculate surface pressure from prx, pry,,txp,typ,pct
   ! Modified by James Done Mar 5 2021

   subroutine cal_slp(slp, IL, JL, ds, prx, pry, txc, tyc, pc)
      implicit none
      integer ii, jj, IL, JL, ic, jc
      real txc, tyc, pc, ds
      real slp(IL, JL), prx(IL, JL), pry(IL, JL)

      ic = int(txc)
      jc = int(tyc)

      slp(ic, jc) = pc
   !      Standard Atmosphere states the density of air is 1.225kg/m3
   !      for hurrican, (P=980hPa, T= 300, R= 287) rho = 1.13
   !       rho = 1.13
   !       Print*, "ds = ", ds
   !       Print*, "ds_rho = ", ds_rho

   !      we divide by 100 to convert from pascals to hPa
   !      we actually divide by 113 because we divide by 100 and by rho

      do ii = 1, ic - 1
         slp(ic-ii, jc) =  slp(ic-ii+1, jc) - prx(ic-ii, jc)*ds/100.
      enddo

      do ii = ic + 1, IL
         slp(ii, jc) =  slp(ii-1, jc) + prx(ii, jc)*ds/100.
      enddo

      do ii = 1, IL
         do jj = 1, jc - 1
            slp(ii, jc-jj) =  slp(ii, jc-jj+1) - pry(ii, jc-jj)*ds/100.
         enddo
         do jj = jc + 1, JL
            slp(ii, jj) =  slp(ii, jj-1) + pry(ii, jj)*ds/100.
         enddo
      enddo

      return
   end
end module cal_slp_mod
