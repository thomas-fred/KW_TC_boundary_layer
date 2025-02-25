A*******************************************************************
c  The model is designed to simulate the boundary layer response to *
c  a presure field of a tropical cyclones in gradient gradient wind *
c  balance at the top of the boundary layer with a multi-level      *
c  hydrostatic, elastic model on an f-plane or on a beta-plan       *
c  The model was designed and coded by Yuqing Wang                  *
c  The code was modified by James Done and Ming Ge in 2019-2021     *
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      program PBL_MODEL
c&&&&&&&&&&&&&&&&&&s&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c Compile for Intel:  ifort -O3 -cpp -qopenmp -ftz -zero -o footprint_model.exe footprint_model.f
c--------------------------------------------------------------------
      implicit none
      integer lq, lp, km, nbd, nx, ny, nx_c, ny_c, n_bdy, n_lu
      parameter(km=18, nbd=10, n_lu=21)
c 
      integer nt,jt0,it0,ji,ij,jtop,itop,jj,ii
      integer i,j,k,kk,inhr,indx,kfi,hr_o
      integer nhr,npl,khr,nst,ntt,nmh
      real y,y0,yct,beta,pi1,ap,xd,yd
      real txp,typ,ds,dt,txpc,typc
      real um1,vm1,rm1,ux,vy,ux0,vy0, ux1, vy1
      real spd_uv, spd_sfc_sym, wind_sfc_sym
      real dir_gl, dir_gl_0, dir_sfc, dir_sfc_0, spd_gl, spd_gl_0
      real bnd(nbd)

      real, allocatable, dimension(:):: vm_1d, pc_1d, rm_1d, r34_1d, 
     +                           pe_1d, typ_1d, txp_1d, lat_1d, lon_1d
  
      real, allocatable, dimension(:,:)::f_2d, prx, pry, prr, bw, slp,
     +                                top, cd, rr_tc, topo, lat_2d, top4
      real, allocatable, dimension(:,:,:):: u, v, w, u4, v4, akv

      integer, allocatable, dimension(:,:):: land, landuse

c  Model variables as usually used
      real zz(km+1),zu(km),dzz(km),dzu(km), sfz0(n_lu)
      real vmt,rmt,r34t,pct,pet
      real tvmt,trmt,tr34t,tpct,tpet
 
      character fnm*20,fname*10, nhr_s*3, nt_s*3, pro_wind*10

      NAMELIST /INIT/lp,lq,ds,dt,kfi,khr,indx,nst,fnm,pro_wind
c   jt0, it0, x- and y-grid indice of TC center at initiall time
c   ds - grid spacing
c   dt - time step
c   kfi - total integration time in hours
c   khr - Output units, 1 means every hour
c   index - 0-f-plane, 1-beta-plane
c   nst - whether the periodic zonal boundary condition is used
c          1 - cyclic zonal boundary condition is used
c          0 - no cyclic zonal boundary condition is used
c   fnm - output fine name prefix
c   pro_wind - wind profile, H10 and WILLOUBY for now

      open(7,file='namelist',form='formatted',status='old')
      read(7,INIT)
      close(7)

      nx=lq
      ny=lp 

      nx_c  = (nx+1)/2
      ny_c  = (ny+1)/2

      ! 10-min
      n_bdy = kfi*6 + 1

      allocate(land(lq,lp), landuse(lq,lp), slp(lq, lp))
      allocate(f_2d(lq,lp),prx(lq,lp),pry(lq,lp),lat_2d(lq,lp),
     +      prr(lq,lp), bw(lq,lp), top(lq,lp), cd(lq,lp), topo(lq,lp))
      allocate(u(lq,lp,km),v(lq,lp,km), akv(lq,lp,km))
      allocate(w(lq,lp,km+1))
      allocate(u4(nx,ny,km),v4(nx,ny,km),rr_tc(nx,ny),top4(nx,ny))
      allocate(typ_1d(n_bdy),txp_1d(n_bdy))
      allocate(vm_1d(n_bdy),pc_1d(n_bdy), rm_1d(n_bdy), r34_1d(n_bdy))
      allocate(pe_1d(n_bdy),lat_1d(n_bdy),lon_1d(n_bdy))

c     read time series of TC info. from ebtrk every 10 minutes
      open(17,file='bdy_10min.txt', form='formatted', status='old')
      read(17, 27)
      do k = 1, n_bdy
        read(17, 27) txp_1d(k), typ_1d(k), vm_1d(k), pc_1d(k)
     +             , rm_1d(k), r34_1d(k), pe_1d(k), lat_1d(k), lon_1d(k)
      enddo
      close(17)
  27  format(2f10.2, 5f10.1, 2f10.2 ) 
c JAMES check read in is ok   
      Print*, txp_1d
      Print*, typ_1d
      Print*, vm_1d
      Print*, pc_1d
      Print*, rm_1d
      Print*, r34_1d
      Print*, pe_1d
      Print*, lat_1d
      Print*, lon_1d
c JAMES check read in is ok   



c JAMES test for N. Hemisphere storms
c      lat_1d=lat_1d*-1.
c      Print*, lat_1d
c JAMES test for N. Hemisphere storms    


c-----------------------------------
c to read lat-lon, landuse and terrain data 
c-----------------------------------
      open (16, file="lat_2d.dat", form="unformatted")
      read (16) lat_2d
       
c JAMES test for N. Hemisphere storms
c      lat_2d = lat_2d*-1.
c JAMES test for N. Hemisphere storms

      open (18, file="topo.dat", form="unformatted")
      read (18) topo

      open (19, file="landuse.dat", form="unformatted")
      read (19) landuse

      open(20, file='landuse.tbl', form='formatted', status='old')
      read(20, 29)
      read(20, 29)
      read(20, 29)
      do k = 1, n_lu
        read(20, 29) sfz0(k)
      enddo
      close(20)
29    format(27x, f4.0)


      do i=1,lp
      do j=1,lq
        top(j,i)= real(topo(j,i))
        land(j,i)=0
        if(top(j,i).le.0.12) top(j,i)=0.0
        if(top(j,i).gt.0.12) land(j,i)=1
        cd(j,i)=0.4*0.4*log(10.0/(sfz0(landuse(j,i))*0.01))**(-2) 
      enddo
      enddo

c if rmax is missing,calculates Rmax using Eq. 7a from Willoughby et al. (2006):
c fix spd_gl 
      do k = 1, n_bdy
        if(rm_1d(k).ge.800.0e3.or.rm_1d(k).le.9.0e3) then
          !spd_gl_0 =spd_gl(vm_1d(k),land(int(txp_1d(k)),int(txp_1d(k))))
          spd_gl_0 =vm_1d(k)/0.9
          rm_1d(k) = max(rm_1d(k), 
     +         46.4 *exp(-0.0155*spd_gl_0+ 0.0169*abs(lat_1d(k)))*1.0e3)
        end if
      end do
      
      zz(1)=0.0
      do k=2,km+1
!       zz(k)=zz(k-1)+20.0*(real(k-1))**0.5
!  Yuqing version
!       zz(k)=zz(k-1)+20.0*(real(k-1))**0.67026
! Ming version   if topo is too hight, change following parameter     
! 0.9->2550  0.95 -> 2875  0.98->3088
      zz(k)=zz(k-1)+20.0*(real(k-1))**0.98
      enddo
      do k=1,km
        zu(k)=0.5*(zz(k)+zz(k+1))
        dzz(k)=zz(k+1)-zz(k)
        write(6,*) zu(k)
      enddo
      do k=1,km-1
        dzu(k)=zu(k+1)-zu(k)
      enddo
      dzu(km)=2.0*zz(km+1)-zu(km)

c     inhr: time steps in an hour
c     hr_o: time steps for 10-minute output hr_o 
      inhr = int(3600.0/dt + 0.1)
      hr_o = int(60*10/dt + 0.1)

c--------------------------------
c to calculate Coriolis parameter
c--------------------------------
c f = 2*omega*sin(phi)
      pi1=atan(1.d0)/45.0d0
      do i=1,lp
      do j=1,lq
        f_2d(j,i) = 1.4584d-4*sin(lat_2d(j,i)*pi1)  
      end do
      end do

c----------------------------------------------------------
c to specify the lateral boundary layer buffer zone width
c----------------------------------------------------------
      do kk=1,nbd
        bnd(kk)=0.125d0*real(nbd-kk)/real(nbd-1)
      enddo
c-----------------------------------------------------------
c lateral boundary damping coefficient in the buffer zone
c-----------------------------------------------------------
      do 80 kk=1,nbd
        ap=(real(nbd-kk)/real(nbd-1))*sqrt(real(nbd-kk+1))
      do 80 i=kk,lp-kk+1
      if(nst.eq.1) then
        do j=1,lq
          bw(j,i)=1.d0+ap
        enddo
      else
        do j=kk,lq-kk+1
          bw(j,i)=1.d0+ap
        enddo
      endif
  80  continue

c-----------------------------------------
c  to calculate the initial moving speed of the TC
c-----------------------------------------
      nhr   = 0
      ux    = (txp_1d(2) - txp_1d(1))*ds/600.0
      vy    = (typ_1d(2) - typ_1d(1))*ds/600.0
      tvmt  = (vm_1d(2)  - vm_1d(1))/600.0
      trmt  = (rm_1d(2)  - rm_1d(1))/600.0
      tr34t = (r34_1d(2) - r34_1d(1))/600.0
      tpct  = (pc_1d(2)  - pc_1d(1))/600.0
      tpet  = (pe_1d(2)  - pe_1d(1))/600.0

      ux0=ux
      vy0=vy

c-----------------------------------------------
c  to initialise the cyclone centre at each mesh
c-----------------------------------------------
      txp =txp_1d(1)
      typ =typ_1d(1)
      vmt =vm_1d(1)
      rmt =rm_1d(1)
      r34t=r34_1d(1)
      pct =pc_1d(1)
      pet =pe_1d(1)
     
c------------------------------------------------------------
c to initialize the model wind & pressure gradient fields
c-------------------------------------------------------------
      akv=0.0
      w=0.0

      if(pro_wind.eq."H10") then
        call H10(u,v,f_2d,lq,lp,km,ds,vmt,rmt,txp,typ,
     +        r34t, pct, pet, ux, vy)
      else if (pro_wind.eq.'WILLOUBY')then
        call willouby(u,v,f_2d,lq,lp,km,ds,vmt,rmt,txp,typ,ux,vy,
     +                     lat_1d(1),land)
      else
        print*, "not a valid wind profile"
        stop
      end if

c-----------------------------------------------------------------------
      open(10,file=trim(fnm)//'.dat',form='formatted',status='unknown')
c
c  Start time integration cycles
c
      nmh = 1
      ntt = 1
      nt  = nhr*inhr
c
  100 nt = nt+1 
c
      if(mod(nt-1,inhr).eq.0) then
        nhr=(nt-1)/inhr
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  to determine the cyclone parameters such as center, maximum wind.
c-------------------------------------------------------------------
        txpc=txp_1d(nmh)
        typc=typ_1d(nmh)
        
        call centre(u,v,lq,lp,km,vm1,um1,rm1,txpc,typc,ds)
      
        write(6, 121) nhr,txpc,typc,vm1,rm1
        write(10,121) nhr,txpc,typc,vm1,rm1
 121    format(1x,'hr=',i3,', txp=',f8.3,', typ=',f8.3,', vm=',f5.1,
     +         ', rm=',f6.1)

      end if

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  to write the results during the time integration such as K.E
c-------------------------------------------------------------------
        ! update every 10-minute
        if(mod(nt-1, hr_o).eq.0.and.nmh.lt.n_bdy.and.nt.ne.1) then
          nmh  = nmh+1
          ux   = (txp_1d(nmh+1)-txp_1d(nmh-1))*ds/1200.0
          vy   = (typ_1d(nmh+1)-typ_1d(nmh-1))*ds/1200.0
          tvmt = (vm_1d(nmh+1)-vm_1d(nmh-1))/1200.0
          trmt = ( rm_1d(nmh+1)-rm_1d(nmh-1))/1200.0
          tr34t= (r34_1d(nmh+1)-r34_1d(nmh-1))/1200.0
          tpct = (pc_1d(nmh+1)-pc_1d(nmh-1))/1200.0
          tpet = (pe_1d(nmh+1)-pe_1d(nmh-1))/1200.0
c         
          txp  = txp_1d(nmh) + ux*dt/ds
          typ  = typ_1d(nmh) + vy*dt/ds
          vmt  = vm_1d (nmh) + tvmt*dt
          rmt  = rm_1d (nmh) + trmt*dt
          r34t = r34_1d(nmh) + tr34t*dt
          pct  = pc_1d (nmh) + tpct*dt
          pet  = pe_1d (nmh) + tpet*dt
        else
          txp  = txp + ux*dt/ds
          typ  = typ + vy*dt/ds
          vmt  = vmt + tvmt*dt
          rmt  = rmt + trmt*dt
          r34t = r34t + tr34t*dt
          pct  = pct + tpct*dt
          pet  = pet + tpet*dt
        end if

c      write out after 10 minute      
       if(mod(nt-1, hr_o).eq.0) then

c        call cal_slp(slp, lq, lp, ds, prx, pry, txp, typ, pct)

         if(pro_wind.eq."H10") then
           call H10_slp(slp,lq,lp,ds,vmt,rmt,txp,typ,pct,pet)
         end if


        do i=1,ny
        do j=1,nx
          jj = j
            ii = i
 
          top4(j,i)=top(jj,ii)
          do k=1,km
            u4(j,i,k)=u(jj,ii,k)
            v4(j,i,k)=v(jj,ii,k)
          enddo

        end do
        end do
c     
        write (nt_s, '(I3.3)') (nt-1)/hr_o
        open(21, file=trim(fnm)//'_'//trim(pro_wind)//"_"//trim(nt_s)
     +     //'.d', form='unformatted')

        do k=1,1
          write(21) ((u4(j,i,k),j=1,nx),i=1,ny)
        enddo
        do k=1,1
          write(21) ((v4(j,i,k),j=1,nx),i=1,ny)
        enddo

        write(21) ((top4(j,i),j=1,nx),i=1,ny)
 
        write(21) ((cd(j,i),j=1,nx),i=1,ny)

        write(21) ((slp(j,i),j=1,nx),i=1,ny)

c        write(21) ((prx(j,i),j=1,nx),i=1,ny)

c        write(21) ((pry(j,i),j=1,nx),i=1,ny)

        close(21)
      endif

      if(nhr.lt.kfi) then
c       wind update 10min at 4,9,14,19
c       if(mod(nt, hr_o).eq.hr_o/2) then
        if(mod(nt, hr_o).eq.hr_o/2) then
          if(pro_wind.eq."H10") then
            call H10(u,v,f_2d,lq,lp,km,ds,vmt,rmt,txp,typ,
     +        r34t, pct, pet, ux, vy)
          else if (pro_wind.eq.'WILLOUBY')then
            call willouby(u,v,f_2d,lq,lp,km,ds,vmt,rmt,txp,typ,ux,vy,
     +                    lat_1d(1),land)
          end if
        end if

        if(pro_wind.eq."H10") then
          call prxy_H10(u,v,prx,pry,f_2d,lq,lp,km, ds,vmt,rmt,
     +    r34t, txp, typ, pct, pet, ux, vy)
        else
          call prxy_willouby(prx,pry,f_2d,lq,lp,km,ds,vmt,rmt,txp,typ
     +                     ,ux,vy,lat_1d(nmh),land)
        end if

       
        call inner(u,v,w,f_2d,prx,pry,akv,lq,lp,ds,dt,top,
     +     land,cd,bw,bnd,km,zz,dzz,zu,dzu,nbd,ntt,nst)

        goto 100
      end if

      deallocate(f_2d, prx, pry, lat_2d, prr, bw, top, cd, topo)
      deallocate(u, v, akv, w)
      deallocate(u4,v4,rr_tc,top4,slp)
      deallocate(typ_1d,txp_1d)
      deallocate(vm_1d,pc_1d, rm_1d, r34_1d,pe_1d,lat_1d,lon_1d)

      close(10)
      stop 
      end 
c Ming Ge Feb. 2 2021 calculate surface pressure from prx, pry,,txp,typ,pct
c Modified by James Done Mar 5 2021

       subroutine cal_slp(slp, IL, JL, ds, prx, pry, txc, tyc, pc)
       implicit none
       integer ii, jj, IL, JL, ic, jc
       real txc, tyc, pc, ds
       real slp(IL, JL), prx(IL, JL), pry(IL, JL)

       ic = int(txc)
       jc = int(tyc)

       slp(ic, jc) = pc
c      Standard Atmosphere states the density of air is 1.225kg/m3 
c      for hurrican, (P=980hPa, T= 300, R= 287) rho = 1.13
c       rho = 1.13 
c       Print*, "ds = ", ds
c       Print*, "ds_rho = ", ds_rho

c      we divide by 100 to convert from pascals to hPa
c      we actually divide by 113 because we divide by 100 and by rho
      
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

c





c
        subroutine inner(u,v,w,f_2d,prx,pry,akv,IL,JL,ds,dt,top,land,
     +             cd, bw,bnd,km,zz,dzz,zu,dzu,nbn,ntt,nst)
        implicit none
        integer i,j,k,IL,JL,km,nbn,ntt,nst
        real uu,vv,ds,dt,vadu,vadv,coru,corv
        real ww,wwa,wwb
        real u(IL,JL,km),v(IL,JL,km),prx(IL,JL),pry(IL,JL),f_2d(IL,JL) 
        real up(IL,JL,km),vp(IL,JL,km),akh(IL,JL,km),akv(IL,JL,km)
        real ut(IL,JL,km),vt(IL,JL,km),bw(IL,JL),bnd(nbn)
        real w(IL,JL,km+1),zz(km+1),dzz(km),zu(km),dzu(km)
        real hadu(IL,JL,km),hadv(IL,JL,km),top(IL,JL),cd(IL,JL)
        integer land(IL,JL)

c       inner performs explicit time integration
c       computes time tendencies for all variables:
c       scheme employed for time-differencing is due to
c       m. miller (qjrms,1978) or matsuno (via option)
c
       ut=0.0
       vt=0.0





       call hdiff(ut,vt,u,v,akh,bw,IL,JL,km,ds,dt,nst)


       call vdiff(ut,vt,u,v,akv,IL,JL,km,dt,zz,dzz,dzu,top,land,cd)


       call hadva(hadu,u,v,u,IL,JL,km,ds,nst)


       call hadva(hadv,u,v,v,IL,JL,km,ds,nst)



c$omp parallel do default(shared)
c$omp& private(i,j,k,uu,vv,vadu,vadv,coru,corv,ww,wwa,wwb)
        do i=1,IL
         do k=1,km
          do j=1,JL
           uu=u(i,j,k)
           vv=v(i,j,k)
c    vertical advection first-order upwind scheme
          ww=0.5*(w(i,j,k)+w(i,j,k+1))
          wwa=0.5*(ww-abs(ww))/dzu(k)
          wwb=0.5*(ww+abs(ww))/dzu(k-1)
          if(k.eq.1) then
           vadu=-wwa*(u(i,j,k+1)-uu)
           vadv=-wwa*(v(i,j,k+1)-vv)
          else if(k.eq.km) then
           vadu=-wwb*(uu-u(i,j,k-1))
           vadv=-wwb*(vv-v(i,j,k-1))
          else
           vadu=-wwa*(u(i,j,k+1)-uu)-wwb*(uu-u(i,j,k-1))
           vadv=-wwa*(v(i,j,k+1)-vv)-wwb*(vv-v(i,j,k-1))
          endif
c   coriolis term
           coru= f_2d(i,j)*vv
           corv=-f_2d(i,j)*uu
c   update variables
           up(i,j,k)=uu+dt*(hadu(i,j,k)+vadu-prx(i,j)+coru+ut(i,j,k))
           vp(i,j,k)=vv+dt*(hadv(i,j,k)+vadv-pry(i,j)+corv+vt(i,j,k))
         enddo
        enddo
      enddo
c$omp end parallel do
      if (ntt.eq.1) then
c    update u,v
c$omp parallel do default(shared)
c$omp& private(i,j,k)
       do i=1,IL
        do j=1,JL
         do k=1,km
           u(i,j,k)=up(i,j,k)
           v(i,j,k)=vp(i,j,k)
          enddo
         enddo
       enddo
c$omp end parallel do
       ntt=0
       go to 8
      endif
c
c    for even time-steps complete a matsuno time integration

       call www(w,up,vp,top,IL,JL,km,ds,zz,zu,dzz,nst)

       call hadva(hadu,up,vp,up,IL,JL,km,ds,nst)

       call hadva(hadv,up,vp,vp,IL,JL,km,ds,nst)

c$omp parallel do default(shared)
c$omp& private(i,j,k,uu,vv,vadu,vadv,coru,corv,ww,wwa,wwb)
      do i=1,IL
       do k=1,km
        do j=1,JL
         uu=up(i,j,k)
         vv=vp(i,j,k)
c   vertical advection first-order upwind sckeme
         ww=0.5*(w(i,j,k)+w(i,j,k+1))
         wwa=0.5*(ww-abs(ww))/dzu(k)
         wwb=0.5*(ww+abs(ww))/dzu(k-1)
         if(k.eq.1) then
          vadu=-wwa*(up(i,j,k+1)-uu)
          vadv=-wwa*(vp(i,j,k+1)-vv)
         else if(k.eq.km) then
          vadu=-wwb*(uu-up(i,j,k-1))
          vadv=-wwb*(vv-vp(i,j,k-1))
         else
          vadu=-wwa*(up(i,j,k+1)-uu)-wwb*(uu-up(i,j,k-1))
          vadv=-wwa*(vp(i,j,k+1)-vv)-wwb*(vv-vp(i,j,k-1))
         endif
c   coriolis terms
        coru= f_2d(i,j)*vv
        corv=-f_2d(i,j)*uu
c   update variables
        u(i,j,k)=u(i,j,k)+dt*(hadu(i,j,k)+vadu-prx(i,j)+coru+ut(i,j,k))
        v(i,j,k)=v(i,j,k)+dt*(hadv(i,j,k)+vadv-pry(i,j)+corv+vt(i,j,k))
        enddo
       enddo
      enddo
c$omp end parallel do

       ntt=1

  8   continue
c
       call bound(u,IL,JL,bnd,nbn,km,nst)
       call bound(v,IL,JL,bnd,nbn,km,nst)
       call www(w,u,v,top,IL,JL,km,ds,zz,zu,dzz,nst)
       return                                                      
       end    
c
c----------------------------------------------------------
c  To calculate horizontal advection tendency of u, v winds 
c----------------------------------------------------------

       subroutine hadva(had,u,v,a,IL,JL,km,ds,nnest)

c----------------------------------------------------------
       implicit none
       real ds,dsi,C1,C2
       parameter(C1=1.0/6.0,C2=1.0/60.0)
       integer i,j,k,IL,JL,km,nnest,im1,im2,im3,ip1,ip2,ip3
       real advx(IL),u(IL,JL,km),a(IL,JL,km),had(IL,JL,km)
       real advy(JL),v(IL,JL,km)

       dsi=1.0/ds

       if(nnest.eq.1) then

c$omp parallel do default(shared)
c$omp& private(i,j,k,advx,im1,im2,im3,ip1,ip2,ip3)
        do j=1,JL
         do k=1,km
           do i=1,IL
           ip1=i+1
           ip2=i+2
           ip3=i+3
           im1=i-1
           im2=i-2
           im3=i-3
           if(i.eq.1) then
           im1=IL
           im2=IL-1
           im3=IL-2
           endif
           if(i.eq.2) then
           im2=IL
           im3=IL-1
           endif
           if(i.eq.3) then
           im3=IL
           endif
           if(i.eq.IL) then
           ip1=1
           ip2=2
           ip3=3
           endif
           if(i.eq.IL-1) then
           ip2=1
           ip3=2
           endif
           if(i.eq.IL-2) then
           ip3=1
           endif
           if(u(i,j,k).ge.0.0) then
            advx(i)=(-dsi*C2*u(i,j,k))*(-3.0*a(ip2,j,k)+30.0*
     +          a(ip1,j,k)+20.0*a(i,j,k)-60.0*a(im1,j,k)+15.0
     +          *a(im2,j,k)-2.0*a(im3,j,k))  
           else
            advx(i)=(-dsi*C2*u(i,j,k))*(2.0*a(ip3,j,k)-15.0
     +          *a(ip2,j,k)+60.0*a(ip1,j,k)-20.0*a(i,j,k)-30.0*
     +          a(im1,j,k)+3.0*a(im2,j,k))
           endif
         enddo
         do i=1,IL
           had(i,j,k)=advx(i)
         enddo
c
         enddo
         enddo
c$omp end parallel do

      else

c$omp parallel do default(shared)
c$omp& private(i,j,k,advx)
        do j=1,JL
         do k=1,km
          if (u(1,j,k).ge.0.0) then
          advx(1)=0.0
          else
          advx(1)=-u(1,j,k)*(a(2,j,k)-a(1,j,k))*dsi
          endif
          if (u(2,j,k).ge.0.0) then
          advx(2)=-u(2,j,k)*(a(2,j,k)-a(1,j,k))*dsi
          else
          advx(2)=-u(2,j,k)*(a(3,j,k)-a(2,j,k))*dsi
          endif
          if (u(3,j,k).ge.0.0) then
          advx(3)=-u(3,j,k)*(2.0*a(4,j,k)+3.0*a(3,j,k)-6.0*a(2,j,k)
     +         +a(1,j,k))*C1*dsi
          else
          advx(3)=-u(3,j,k)*(-2.0*a(2,j,k)-3.0*a(3,j,k)+6.0*a(4,j,k)
     +         -a(5,j,k))*C1*dsi
          endif
          if (u(IL,j,k).le.0.0) then
          advx(IL)=0.0
          else
          advx(IL)=-u(IL,j,k)*(a(IL,j,k)-a(IL-1,j,k))*dsi
          endif
          if (u(IL-1,j,k).ge.0.0) then
          advx(IL-1)=-u(IL-1,j,k)*(a(IL-1,j,k)-a(IL-2,j,k))*dsi
          else
          advx(IL-1)=-u(IL-1,j,k)*(a(IL,j,k)-a(IL-1,j,k))*dsi
          endif
          if (u(IL-2,j,k).ge.0.0) then
          advx(IL-2)=-u(IL-2,j,k)*(2.0*a(IL-1,j,k)+3.0*a(IL-2,j,k)-
     +         6.0*a(IL-3,j,k)+a(IL-4,j,k))*C1*dsi
          else
          advx(IL-2)=-u(IL-2,j,k)*(-2.0*a(IL-3,j,k)-3.0*a(IL-2,j,k)+
     +         6.0*a(IL-1,j,k)-a(IL,j,k))*C1*dsi
          endif
           do i=4,IL-3
           if(u(i,j,k).ge.0.0) then
            advx(i)=(-dsi*C2*u(i,j,k))*(-3.0*a(i+2,j,k)+30.0*
     +          a(i+1,j,k)+20.0*a(i,j,k)-60.0*a(i-1,j,k)+15.0
     +          *a(i-2,j,k)-2.0*a(i-3,j,k))  
           else
            advx(i)=(-dsi*C2*u(i,j,k))*(2.0*a(i+3,j,k)-15.0
     +          *a(i+2,j,k)+60.0*a(i+1,j,k)-20.0*a(i,j,k)-30.0*
     +          a(i-1,j,k)+3.0*a(i-2,j,k))
           endif
         enddo
         do i=1,IL
           had(i,j,k)=advx(i)
         enddo
c
         enddo
         enddo
c$omp end parallel do
 
      endif

c$omp parallel do default(shared)
c$omp& private(i,j,k,advy)
        do i=1,IL
         do k=1,km
          if (v(i,1,k).ge.0.0) then
          advy(1)=0.0
          else
          advy(1)=-v(i,1,k)*(a(i,2,k)-a(i,1,k))*dsi
          endif
          if (v(i,2,k).ge.0.0) then
          advy(2)=-v(i,2,k)*(a(i,2,k)-a(i,1,k))*dsi
          else
          advy(2)=-v(i,2,k)*(a(i,3,k)-a(i,2,k))*dsi
          endif
          if (v(i,3,k).ge.0.0) then
          advy(3)=-v(i,3,k)*(2.0*a(i,4,k)+3.0*a(i,3,k)-6.0*a(i,2,k)
     +         +a(i,1,k))*C1*dsi
          else
          advy(3)=-v(i,3,k)*(-2.0*a(i,2,k)-3.0*a(i,3,k)+6.0*a(i,4,k)
     +         -a(i,5,k))*C1*dsi
          endif
          if (v(i,JL,k).le.0.0) then
          advy(JL)=0.0
          else
          advy(JL)=-v(i,JL,k)*(a(i,JL,k)-a(i,JL-1,k))*dsi
          endif
          if (v(i,JL-1,k).ge.0.0) then
          advy(JL-1)=-v(i,JL-1,k)*(a(i,JL-1,k)-a(i,JL-2,k))*dsi
          else
          advy(JL-1)=-v(i,JL-1,k)*(a(i,JL,k)-a(i,JL-1,k))*dsi
          endif
          if (v(i,JL-2,k).ge.0.0) then
          advy(JL-2)=-v(i,JL-2,k)*(2.0*a(i,JL-1,k)+3.0*a(i,JL-2,k)-
     +         6.0*a(i,JL-3,k)+a(i,JL-4,k))*C1*dsi
          else
          advy(JL-2)=-v(i,JL-2,k)*(-2.0*a(i,JL-3,k)-3.0*a(i,JL-2,k)+
     +         6.0*a(i,JL-1,k)-a(i,JL,k))*C1*dsi
          endif
           do j=4,JL-3
           if(v(i,j,k).ge.0.0) then
            advy(j)=(-dsi*C2*v(i,j,k))*(-3.0*a(i,j+2,k)+30.0*
     +          a(i,j+1,k)+20.0*a(i,j,k)-60.0*a(i,j-1,k)+15.0
     +          *a(i,j-2,k)-2.0*a(i,j-3,k))
           else
            advy(j)=(-dsi*C2*v(i,j,k))*(2.0*a(i,j+3,k)-15.0
     +          *a(i,j+2,k)+60.0*a(i,j+1,k)-20.0*a(i,j,k)-30.0*
     +          a(i,j-1,k)+3.0*a(i,j-2,k))
           endif
         enddo
         do j=1,JL
           had(i,j,k)=had(i,j,k)+advy(j)
         enddo
c
         enddo
         enddo
c$omp end parallel do

      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine hdiff(ut,vt,u,v,akh,bw,lq,lp,km,ds,dt,nnest)
c
c---------------------------------------------------------------
      implicit none
      integer i,j,k,lq,lp,km,nnest,jm1,jp1
      real dsi2,ds,dt,akmax2,ag,hd00,hdcoe,dm1,dm2,dm,ah
      real ut(lq,lp,km),vt(lq,lp,km),u(lq,lp,km),v(lq,lp,km)
      real akh(lq,lp,km),bw(lq,lp)
c
      ah=0.5d0*0.4d0*ds
      hd00=3.0d-3*ds*ds/dt
      hdcoe=(ah*ah)/(2.0d0*ds)
      dsi2=1.0d0/(2.0d0*ds*ds)
      akmax2=ds*ds/(8.d0*dt)
c-------------------------------------------------------------
c  To include linear 4th-order horizontal diffusion term using 
c  the t-1 time level fields to ensure the numerical stability
c-------------------------------------------------------------
c
      if(nnest.eq.1) then
c
c$omp parallel do default(shared)
c$omp& private(j,i,k,jm1,jp1,dm1,dm2,dm)
      do 10 i=2,lp-1
      do 10 k=1,km
      do 10 j=1,lq
      jp1=j+1
      jm1=j-1
      if(j.eq.1) jm1=lq
      if(j.eq.lq) jp1=1
      dm1=u(jp1,i,k)-u(jm1,i,k)-v(j,i+1,k)+v(j,i-1,k)
      dm2=v(jp1,i,k)-v(jm1,i,k)+u(j,i+1,k)-u(j,i-1,k)
      dm=sqrt(dm1*dm1+dm2*dm2)
      akh(j,i,k)=min(akmax2,bw(j,i)*(hd00+hdcoe*dm))
 10   continue
c$omp end parallel do
c
c$omp parallel do default(shared)
c$omp& private(j,i,k)
      do j=1,lq
      do k=1,km
      akh(j,1,k)=akh(j,2,k)
      akh(j,lp,k)=akh(j,lp-1,k)
      enddo
      enddo
c$omp end parallel do

c$omp parallel do default(shared)
c$omp& private(j,i,k,jm1,jp1)
      do 20 i=2,lp-1
      do 20 k=1,km
      do 20 j=1,lq
        jp1=j+1
        jm1=j-1
      if(j.eq.1) jm1=lq
      if(j.eq.lq) jp1=1
      ut(j,i,k)=ut(j,i,k)+dsi2*((akh(jp1,i,k)+akh(j,i,k))*(u(jp1,i,k)
     +     -u(j,i,k))-(akh(j,i,k)+akh(jm1,i,k))*(u(j,i,k)-u(jm1,i,k))
     +     +(akh(j,i+1,k)+akh(j,i,k))*(u(j,i+1,k)-u(j,i,k))-(akh(j,i,
     +     k)+akh(j,i-1,k))*(u(j,i,k)-u(j,i-1,k)))
      vt(j,i,k)=vt(j,i,k)+dsi2*((akh(jp1,i,k)+akh(j,i,k))*(v(jp1,i,k)
     +     -v(j,i,k))-(akh(j,i,k)+akh(jm1,i,k))*(v(j,i,k)-v(jm1,i,k))
     +     +(akh(j,i+1,k)+akh(j,i,k))*(v(j,i+1,k)-v(j,i,k))-(akh(j,i,
     +     k)+akh(j,i-1,k))*(v(j,i,k)-v(j,i-1,k)))
 20   continue
c$omp end parallel do

      else 
c
c$omp parallel do default(shared)
c$omp& private(j,i,k,dm1,dm2,dm)
      do 30 i=2,lp-1
      do 30 k=1,km
      do 30 j=2,lq-1
      dm1=u(j+1,i,k)-u(j-1,i,k)-v(j,i+1,k)+v(j,i-1,k)
      dm2=v(j+1,i,k)-v(j-1,i,k)+u(j,i+1,k)-u(j,i-1,k)
      dm=sqrt(dm1*dm1+dm2*dm2)
      akh(j,i,k)=min(akmax2,bw(j,i)*(hd00+hdcoe*dm))
 30   continue
c$omp end parallel do
c
c$omp parallel do default(shared)
c$omp& private(j,i,k)
      do k=1,km
      do i=2,lp-1
      akh(1,i,k)=akh(2,i,k)
      akh(lq,i,k)=akh(lq-1,i,k)
      enddo
      do j=1,lq
      akh(j,1,k)=akh(j,2,k)
      akh(j,lp,k)=akh(j,lp-1,k)
      enddo
      enddo
c$omp end parallel do

c$omp parallel do default(shared)
c$omp& private(j,i,k)
      do 40 i=2,lp-1
      do 40 k=1,km
      do 40 j=2,lq-1
      ut(j,i,k)=ut(j,i,k)+dsi2*((akh(j+1,i,k)+akh(j,i,k))*(u(j+1,i,k)
     +     -u(j,i,k))-(akh(j,i,k)+akh(j-1,i,k))*(u(j,i,k)-u(j-1,i,k))
     +     +(akh(j,i+1,k)+akh(j,i,k))*(u(j,i+1,k)-u(j,i,k))-(akh(j,i,
     +     k)+akh(j,i-1,k))*(u(j,i,k)-u(j,i-1,k)))
      vt(j,i,k)=vt(j,i,k)+dsi2*((akh(j+1,i,k)+akh(j,i,k))*(v(j+1,i,k)
     +     -v(j,i,k))-(akh(j,i,k)+akh(j-1,i,k))*(v(j,i,k)-v(j-1,i,k))
     +     +(akh(j,i+1,k)+akh(j,i,k))*(v(j,i+1,k)-v(j,i,k))-(akh(j,i,
     +     k)+akh(j,i-1,k))*(v(j,i,k)-v(j,i-1,k)))
  40  continue
c$omp end parallel do
c
      endif

      return
      end
c
c
      subroutine vdiff(ut,vt,u,v,ak,IL,JL,km,dt,zz,dzz,dzu,top,land,cd)
        implicit none
        real al0, an
        parameter(al0=150.0,an=1.0e-6)
        integer i,j,k,IL,JL,km
        real ut(IL,JL,km),vt(IL,JL,km),u(IL,JL,km)
        real v(IL,JL,km),ak(IL,JL,km),zz(km+1), cd(IL,JL)
        real pu(JL,km),pv(JL,km),top(IL,JL)
        real dzz(km),dzu(km),akm(JL,km)
        real ax(JL,km),bx(JL,km),cx(JL,km)
        real xa(JL,km),xb(JL,km),xc(JL,km),dtdzz(km)
        real spd,al,al2,sh,z,dt,hh
        integer land(IL,JL)

        hh=zz(km+1)
        do k=1,km
          dtdzz(k)=dt/dzz(k)
        enddo
c$omp parallel do default(shared)
c$omp& private(i,j,k,spd,pu,pv,akm,sh,ax,bx,cx,xa,xb,xc,z,al,al2)
        do i=1,IL
          do k=2,km
          do j=1,JL
c         z=        top(i,j)+zz(k)*(1.0-top(i,j)/hh)
          z=min(hh, top(i,j)+zz(k)*(1.0-min(top(i,j)/hh,0.90)))
          al=0.4*z/(1.0+0.4*z/al0)
          al2=al*al
c    &       (dzu(k-1)*(1.0-top(i,j)/zz(km+1)))**2
          sh=((u(i,j,k)-u(i,j,k-1))**2+(v(i,j,k)-v(i,j,k-1))**2)/
     &       (dzu(k-1)*(1.0-min(top(i,j)/zz(km+1),0.90)))**2
          akm(j,k)=al2*sqrt(max(0.0,sh-an))
c    &         /(1.0-top(i,j)/hh)**2
          akm(j,k)=max(1.0e-3*dzu(k-1),akm(j,k))
     &         /(1.0-min(top(i,j)/hh,0.90))**2
          enddo
          enddo
          do j=1,JL
          spd=sqrt(u(i,j,1)*u(i,j,1)+v(i,j,1)*v(i,j,1))
          if(land(i,j).eq.0) then
            cd(i,j) = max(1.0, min(2.4,1.0+0.07*(spd-5.0)))*1.0e-3
          end if
c    &         /(1.0-top(i,j)/hh)**2
          akm(j,1)=1.4625*sqrt(cd(i,j))*spd
     &         /(1.0-min(top(i,j)/hh,0.90))**2
          ax(j,1)=-dtdzz(1)*akm(j,2)/dzu(1)
          cx(j,1)=0.0
          bx(j,1)=1.0-ax(j,1)+dtdzz(1)*cd(i,j)*spd
          ax(j,km)=0.0
          cx(j,km)=-dtdzz(km)*akm(j,km)/dzu(km-1)
          bx(j,km)=1.0-cx(j,km)
          do k=2,km-1
          ax(j,k)=-dtdzz(k)*akm(j,k+1)/dzu(k)
          cx(j,k)=-dtdzz(k)*akm(j,k)/dzu(k-1)
          bx(j,k)=1.0-ax(j,k)-cx(j,k)
          enddo
          enddo
          do k=1,km
          do j=1,JL
          pu(j,k)=u(i,j,k)
          pv(j,k)=v(i,j,k)
          xa(j,k)=ax(j,k)
          xb(j,k)=bx(j,k)
          xc(j,k)=cx(j,k)
          enddo
          enddo
          call tridiag(ax,bx,cx,pu,JL,km)
          call tridiag(xa,xb,xc,pv,JL,km)
          do k=1,km
          do j=1,JL
           ut(i,j,k)=ut(i,j,k)+(pu(j,k)-u(i,j,k))/dt
           vt(i,j,k)=vt(i,j,k)+(pv(j,k)-v(i,j,k))/dt
           ak(i,j,k)=akm(j,k)
          enddo
          enddo
        enddo
c$omp end parallel do
        return
        end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine tridiag(ax,bx,cx,yy,lq,km)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  To solve a tri-diagonal matrix using a direct solver
c--------------------------------------------------------------
      implicit none
      integer i,j,k,k1,lq,km
      real ap,yyk1
      real ax(lq,km),bx(lq,km),cx(lq,km),yy(lq,km)
c
      do 10 k=2,km
      k1=k-1
      do 10 j=1,lq
      ap=cx(j,k)/bx(j,k1)
      bx(j,k)=bx(j,k)-ap*ax(j,k1)
      ap=ap*yy(j,k1)
      yy(j,k)=yy(j,k)-ap
  10  continue
c
      do 20 j=1,lq
      yy(j,km)=yy(j,km)/bx(j,km)
 20   continue
c
      do 30 k=km-1,1,-1
      do 30 j=1,lq
      yyk1=yy(j,k+1)
      yy(j,k)=(yy(j,k)-ax(j,k)*yyk1)/bx(j,k)
  30  continue
c
      return
      end
c
      subroutine centre(u,v,lq,lp,km,vm1,um1,rm1,xcy,ycy,ds)

      implicit none
      integer ln, lm, r, it
      parameter(ln=144,lm=51,it=2)
      parameter(r=1)
c
      integer lq,lp,km, k7,k8,k9,j,i, itt
      real ux, a, b, x, phi, ycy0, xcy0, r1, vm,y
      real vm1,um1,rm1,xcy,ycy,ds, rm0, tm1,xxx,yyy,pi1, dd
      real u(lq,lp,km),v(lq,lp,km),us(lq,lp),vs(lq,lp)
      real uu(lm,ln),vv(lm,ln),umax(ln),vmax(ln)

      pi1=4.0*atan(1.0)/180.0
      dd=360.0/real(ln)
c------------------------------------------------------
      do j=1,lp
      do i=1,lq
      us(i,j)=u(i,j,8)
      vs(i,j)=v(i,j,8)
      enddo
      enddo

      do 100 itt=1,it
      xcy0=xcy
      ycy0=ycy
      vm1=0.0
      do 550 k7=1,10
      vm1=0.0
      do 500 k9=1,10
      do 500 k8=1,10
      do 110 j=1,ln
      phi=(j-1)*dd
      do 110 i=1,lm
      r1=r*(i-1)
      x=r1*cos(phi*pi1)
      y=r1*sin(phi*pi1)
      x=x+xcy
      y=y+ycy
      call SCINEX(x,y,us,uu(i,j),lq,lp)
      call SCINEX(x,y,vs,vv(i,j),lq,lp)
c.......................................................
c azimuthal
      a=vv(i,j)*cos(phi*pi1)-uu(i,j)*sin(phi*pi1)
c radial
      b=vv(i,j)*sin(phi*pi1)+uu(i,j)*cos(phi*pi1)
      uu(i,j)=a
      vv(i,j)=b
110   continue
c..............................................................
c calculate the max 
c..............................................................
        vm=0.0
         do i=1,lm
                umax(i)=0.
                vmax(i)=0.
                do j=1,ln
                umax(i)=umax(i)+uu(i,j)
                vmax(i)=vmax(i)+vv(i,j)
                enddo
                umax(i)=umax(i)/real(ln)
                vmax(i)=vmax(i)/real(ln)
                ux=sqrt(vmax(i)**2+umax(i)**2)
                if(umax(i).gt.vm) then
                        vm=umax(i)
                        rm0=i-1.0d0
                        um1=vmax(i)
                        tm1=ux
                endif
         enddo
          if(vm.gt.vm1) then
          vm1=vm
          xxx=xcy
          yyy=ycy
          rm1=rm0
          endif
          xcy=xcy0+(20.0-k9*4.)/real(k7)
          ycy=ycy0+(20.0-k8*4.)/real(k7)
500      continue
          xcy=xxx
          ycy=yyy
          xcy0=xxx
          ycy0=yyy
550       continue
c***********************************************************
100     continue
        rm1=rm1*ds/1000.0

        return 
        end

      SUBROUTINE SCINEX(GM,GN,SCALA,SCINTO,IL,JL)
C THIS SUBROUTINE PRODUCES THE VALUE SCINTO OF A SCALAR FIELD AT A POINT
C GM,GN BY INTERPOLATION OR EXTRAPOLATION OF THE FIELD SCALA  (2-DIRECTI
C BESSEL INTERPOLATION FORMULA). MMIN,MMAX AND NMIN,NMAX ARE THE BOUNDAR
C OF THE GRID ARRAY.
      REAL SCALA(IL,JL)
      MMIN=1
      NMIN=1
      MMAX=IL
      NMAX=JL
      IGM=int(GM)
      JGN=int(GN)
      FM=GM-IGM
      FN=GN-JGN
      IF(FM.LT.1.E-06)FM=0.
      IF(FN.LT.1.E-06)FN=0.
      MS=MMAX-1
      NS=NMAX-1
      MR=MMIN+1
      NR=NMIN+1
      IF(GM.LT.MMAX)GO TO 60
      IF(GN.LT.NMAX)GO TO 20
      E=GM-MMAX
      T1=E*(SCALA(MMAX,NMAX)-SCALA(MS,NMAX))
      E=GN-NMAX
      T2=E*(SCALA(MMAX,NMAX)-SCALA(MMAX,NS))
      SCINTO=SCALA(MMAX,NMAX)+T1+T2
      RETURN
   20 IF(GN.GE.NMIN)GO TO 40
      E=GM-MMAX
      T1=E*(SCALA(MMAX,NMIN)-SCALA(MS,NMIN))
      E=NMIN-GN
      T2=E*(SCALA(MMAX,NMIN)-SCALA(MMAX,NR))
      SCINTO=SCALA(MMAX,NMIN)+T1+T2
      RETURN
   40 P=SCALA(MMAX,JGN)+FN*(SCALA(MMAX,JGN+1)-SCALA(MMAX,JGN))
      H=SCALA(MS,JGN)+FN*(SCALA(MS,JGN+1)-SCALA(MS,JGN))
      E=GM-MMAX
      SCINTO=P+E*(P-H)
      RETURN
   60 IF(GM.GE.MMIN)GO TO 140
      IF(GN.LT.NMAX)GO TO 80
      E=GN-NMAX
      T2=E*(SCALA(MMIN,NMAX)-SCALA(MMIN,NS))
      E=MMIN-GM
      T1=E*(SCALA(MMIN,NMAX)-SCALA(MR,NMAX))
      SCINTO=SCALA(MMIN,NMAX)+T1+T2
      RETURN
   80 IF(GN.GE.NMIN)GO TO 100
      E=NMIN-GN
      T2=E*(SCALA(MMIN,NMIN)-SCALA(MMIN,NR))
      E=MMIN-GM
      T1=E*(SCALA(MMIN,NMIN)-SCALA(MR,NMIN))
      SCINTO=SCALA(MMIN,NMIN)+T1+T2
      RETURN
  100 E=MMIN-GM
      P=SCALA(MMIN,JGN)+FN*(SCALA(MMIN,JGN+1)-SCALA(MMIN,JGN))
      H=SCALA(MR,JGN)+FN*(SCALA(MR,JGN+1)-SCALA(MR,JGN))
      SCINTO=P+E*(P-H)
      RETURN
  120 E=GN-NMAX
      P=SCALA(IGM,NMAX)+FM*(SCALA(IGM+1,NMAX)-SCALA(IGM,NMAX))
      H=SCALA(IGM,NS)+FM*(SCALA(IGM+1,NS)-SCALA(IGM,NS))
      SCINTO=P+E*(P-H)
      RETURN
  140 IF(GN.GE.NMAX)GO TO 120
      IF(GN.GE.NMIN)GO TO 160
      E=NMIN-GN
      P=SCALA(IGM,NMIN)+FM*(SCALA(IGM+1,NMIN)-SCALA(IGM,NMIN))
      H=SCALA(IGM,NR)+FM*(SCALA(IGM+1,NR)-SCALA(IGM,NR))
      SCINTO=P+E*(P-H)
      RETURN
  160 IF(GM.LT.MS.AND.GM.GE.MR.AND.GN.LT.NS.AND.GN.GE.NR)GO TO 180
      P=SCALA(IGM+1,JGN)+FN*(SCALA(IGM+1,JGN+1)-SCALA(IGM+1,JGN))
         H=SCALA(IGM,JGN)+FN*(SCALA(IGM,JGN+1)-SCALA(IGM,JGN))
      SCINTO=H+FM*(P-H)
      RETURN
  180    FQ=0.25*(FM*FM-FM)
      A=SCALA(IGM,JGN-1)+FM*(SCALA(IGM+1,JGN-1)-SCALA(IGM,JGN-1))
     X+FQ*(SCALA(IGM+2,JGN-1)+SCALA(IGM-1,JGN-1)-SCALA(IGM+1,JGN-1)-
     XSCALA(IGM,JGN-1))
      B=SCALA(IGM,JGN)+FM*(SCALA(IGM+1,JGN)-SCALA(IGM,JGN))
     X+FQ*(SCALA(IGM+2,JGN)+SCALA(IGM-1,JGN)-SCALA(IGM+1,JGN)-
     XSCALA(IGM,JGN))
      C=SCALA(IGM,JGN+1)+FM*(SCALA(IGM+1,JGN+1)-SCALA(IGM,JGN+1))
     X+FQ*(SCALA(IGM+2,JGN+1)+SCALA(IGM-1,JGN+1)-SCALA(IGM+1,JGN+1)
     X-SCALA(IGM,JGN+1))
      D=SCALA(IGM,JGN+2)+FM*(SCALA(IGM+1,JGN+2)-SCALA(IGM,JGN+2))
     X+FQ*(SCALA(IGM+2,JGN+2)+SCALA(IGM-1,JGN+2)-SCALA(IGM+1,JGN+2)
     X-SCALA(IGM,JGN+2))
      SCINTO=B+FN*(C-B)+0.25*(FN*FN-FN)*(A+D-B-C)
      RETURN
      END
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine www(w,u,v,top,lq,lp,km,ds,zz,zu,dzz,nnest)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c---------------------------------------------------------------
c
      implicit none
      integer i,j,k,lq,lp,km,jm1,jp1,nnest
      real ds,tdsi,hh,fc,fc2,w1,w2,w3,term
      real u(lq,lp,km),v(lq,lp,km),w(lq,lp,km+1),dzz(km)
      real top(lq,lp),zz(km+1),zu(km)
c
      tdsi=1.d0/(2.0d0*ds)
      hh=zz(km+1)
c
c$omp parallel do default(shared)
c$omp& private(j,i,k,jm1,jp1,w1,w2,w3,term,fc,fc2)
      do 112 i=2,lp-1
      do 112 j=1,lq
      jm1=j-1
      jp1=j+1
      if(nnest.eq.1) then
      if(j.eq.1) jm1=lq
      if(j.eq.lq) jp1=1
      else
      if(j.eq.1) jm1=1
      if(j.eq.1) jp1=3
      if(j.eq.lq) jm1=lq-2
      if(j.eq.lq) jp1=lq
      endif
      w(j,i,1)=0.0
c Ming 0.90 - > 0.90 doesn't help
c     fc=1.0-top(j,i)/hh  
      fc=1.0-min(top(j,i)/hh, 0.90)
      fc2=fc*fc
      do 112 k=1,km
      w1=-tdsi*(u(jp1,i,k)-u(jm1,i,k)+v(j,i+1,k)
     &        -v(j,i-1,k))/fc2
      w2=(u(j,i,k)*(top(jp1,i)-top(jm1,i))+v(j,i,k)*
     &        (top(j,i+1)-top(j,i-1)))*tdsi/(hh*fc2)
      if(k.eq.1) then
      term=((u(j,i,k+1)-u(j,i,k))*(top(jp1,i)-top(jm1,i))+
     +     (v(j,i,k+1)-v(j,i,k))*(top(j,i+1)-top(j,i-1)))
     +     *tdsi/(zu(k+1)-zu(k))
      else if(k.eq.km) then
      term=((u(j,i,k)-u(j,i,k-1))*(top(jp1,i)-top(jm1,i))+
     +     (v(j,i,k)-v(j,i,k-1))*(top(j,i+1)-top(j,i-1)))
     +     *tdsi/(zu(k)-zu(k-1))
      else
      term=((u(j,i,k+1)-u(j,i,k-1))*(top(jp1,i)-top(jm1,i))+
     +     (v(j,i,k+1)-v(j,i,k-1))*(top(j,i+1)-top(j,i-1)))
     +     *tdsi/(zu(k+1)-zu(k-1))
      endif
c Ming
c     w3=-term*(1.0-zu(k)/hh)/(1.0-top(j,i)/hh)
      w3=-term*(1.0-zu(k)/hh)/(1.0-min(top(j,i)/hh,0.90))
      w(j,i,k+1)=w(j,i,k)+(w1+w2+w3)*dzz(k)
 112  continue
c$omp end parallel do
c
      return
      end 
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine bound(t,lq,lp,bnd,nbd,km,nnest)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  to damp the boundary noise for prognostic variable t using special
c  smoothing operation with open lateral boundary conditions
c--------------------------------------------------------------------
      implicit none
      integer i,j,k,lq,lp,km,nnest,nbd
      real t(lq,lp,km),to(lq,lp),bnd(nbd)
c
c$omp parallel do default(shared)
c$omp& private(j,i,k,to)
      do 9 k=1,km
c
      do 4 i=1,lp
      do 4 j=1,lq
      to(j,i)=t(j,i,k)
  4   continue
      if(nnest.eq.1) then
      call bound1c(to,lq,lp,bnd,nbd)
      else
      call bound1(to,lq,lp,bnd,nbd)
      endif
      do 6 i=1,lp
      do 6 j=1,lq
      t(j,i,k)=to(j,i)
  6   continue
c
  9   continue
c$omp end parallel do
c
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine bound1(p,lq,lp,bnd,nbd)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  to damp the boundary noise for p using a special smoothing
c  operation with open boundary conditions
c----------------------------------------------------------------
      implicit none
      integer i,j,n,n1,nq,np,lq,lp,nbd
      real sd,ed
      real p(lq,lp),pw(lp),pe(lp),p1(lq),pp(lq),bnd(nbd)
c
      do 9 n=1,nbd
c
      sd=bnd(n)
      ed=1.d0-4.d0*sd
c
      n1=n+1
      nq=lq-n
      np=lp-n
      do j=n1,nq
      p1(j)=p(j,n1)*ed+sd*(p(j-1,n1-1)+p(j-1,n1+1)+
     &                     p(j+1,n1-1)+p(j+1,n1+1))
      pp(j)=p(j,np)*ed+sd*(p(j-1,np-1)+p(j-1,np+1)+
     &                     p(j+1,np-1)+p(j+1,np+1))
      enddo
c
      n1=n+1
      nq=lq-n
      np=lp-n
      do i=n1+1,np-1
      pw(i)=p(n1,i)*ed+sd*(p(n1-1,i-1)+p(n1+1,i-1)+
     &                     p(n1-1,i+1)+p(n1+1,i+1))
      pe(i)=p(nq,i)*ed+sd*(p(nq-1,i-1)+p(nq+1,i-1)+
     &                     p(nq-1,i+1)+p(nq+1,i+1))
      enddo
c
      do j=n1,nq
      p(j,n1)=p1(j)
      p(j,np)=pp(j)
      enddo
      do i=n1+1,np-1
      p(n1,i)=pw(i)
      p(nq,i)=pe(i)
      enddo
c
  9   continue
c
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine bound1c(p,lq,lp,bnd,nbd)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  to damp the boundary noise for p using a special smoothing
c  operation with open boundary conditions
c----------------------------------------------------------------
      implicit none
      integer i,j,lq,lp,n,n1,np,jm1,jp1,nbd
      real p(lq,lp),p1(lq),pp(lq),bnd(nbd),sd,ed
c
      do j=1,lq
      jm1=j-1
      jp1=j+1
      if(j.eq.1) jm1=lq
      if(j.eq.lq) jp1=1
      p1(j)=0.25d0*(p(jm1,1)+2.d0*p(j,1)+p(jp1,1))
      pp(j)=0.25d0*(p(jm1,lp)+2.d0*p(j,lp)+p(jp1,lp))
      enddo
      do j=1,lq
      p(j,1)=p1(j)
      p(j,lp)=pp(j)
      enddo
c
      do 9 n=1,nbd
c
      sd=bnd(n)
      ed=1.d0-4.d0*sd
      n1=n+1
      np=lp-n
c
      do j=1,lq
      jm1=j-1
      jp1=j+1
      if(j.eq.1) jm1=lq
      if(j.eq.lq) jp1=1
      p1(j)=p(j,n1)*ed+sd*(p(jm1,n1-1)+p(jm1,n1+1)+p(jp1,n1-1)+
     &                     p(jp1,n1+1))
      pp(j)=p(j,np)*ed+sd*(p(jm1,np-1)+p(jm1,np+1)+p(jp1,np-1)+
     &                     p(jp1,np+1))
      enddo
c
      do j=1,lq
      p(j,n1)=p1(j)
      p(j,np)=pp(j)
      enddo
c
  9   continue
c
      return
      end

       subroutine H10(u,v,f_2d,IL,JL,km,ds,vm,rm,txc,tyc,r34,pc,pe,
     +      ux,vy)

       implicit none
       integer i,j,k,IL,JL,km
       real tx,ty,rr,vt,vm,rm,ds,txc,tyc,fac,ux,vy,vm_sfc_sym
       real ee, pe, pc, r34, bs, term1, xx_it, xxterm, term1_34
       real u(IL,JL,km), v(IL,JL,km),f_2d(IL,JL)
       real corr_fac, spd_tc, correction_vt

c      pe = environment pressure value is set to climatology
c      ee = base of natural logarithms
c      sfc density = 1.13 (at 980mb)
c      eqn 7 in Holland 2010. 
c      Assume Vmax and Pmin are reliablly observed
c      It is not critical and a constant value of density can e used.
c       pe = 1010.   
c      fac is a reduction factor from gradient winf to 10-m wind 
c      Ming Ge: March 2018. adjust Vmax by removing forward motion of storm
c      Vmax_sfc_sym = vmax - 0.5*TC_speed
c       vm_sfc_sym = vm - 0.5*sqrt(ux*ux+vy*vy)
       vm_sfc_sym = vm
       ee = 2.718281
       bs = (vm_sfc_sym*vm_sfc_sym*1.13*ee)/(100.*(pe-pc))
       term1_34 = (rm/r34)**bs

c JAMES v2
c Change the reduction factor to 0.76
c       fac = 0.85
c Change reduction factor to 0.71
c      fac = 0.76
       fac = 0.71
c JAMES v2

c  calculate xx_it in eq. 10  at R34, v=17
c  when rr = r34, xxterm = xx_it = 1
      xx_it = 0.5
 10   vt = vm_sfc_sym*(term1_34*exp(1 - term1_34))**xx_it
      if(vt.gt.16.and.xx_it.le.2.0) then
         xx_it = xx_it + 0.02
         go to 10
      end if

      spd_tc = sqrt(ux*ux+vy*vy)

c$omp parallel do default(shared)
c$omp& private(j,i,k,tx,ty,rr,term1,xxterm,vt,corr_fac)
      do j=1,JL
         ty = (real(j)-tyc)*ds
         do i=1,IL
           tx = (real(i)-txc)*ds
           rr = max(1.0,sqrt(ty*ty+tx*tx))
           term1 = (rm/rr)**bs
           if(rr.le.rm) then
             xxterm = 0.5
           else
             xxterm = 0.5 + (rr - rm)*(xx_it-0.5)/(r34 - rm)
           end if

           vt = (vm_sfc_sym/rr)*(term1*exp(1 - term1))**xxterm/fac


c JAMES v2
c Decided to add the assymetry in the post-processing
           corr_fac = 0.
c           corr_fac = 2.0*(rm*rr)/(rm**2 + rr**2)
c     +                 *1.173*spd_tc**0.63/spd_tc
c JAMES v2


c JAMES v2
c add a correction to get the correct sign of the U and V winds in each hemisphere
           if(f_2d(i,j).le.0) then
             do k=1,km
               u(i,j,k)=+vt*ty + corr_fac*ux
               v(i,j,k)=-vt*tx + corr_fac*vy
             enddo
           else
             do k=1,km
               u(i,j,k)=-vt*ty + corr_fac*ux
               v(i,j,k)=+vt*tx + corr_fac*vy
             enddo
           end if
c JAMES v2

c JAMES v2
c remove the special treatment of the lowest four levels.
c           do k=1,4
c             u(i,j,k)=-vt*ty*(fac+(k-1)*(1.-fac)/4.0)
c     +                 + corr_fac*ux
c             v(i,j,k)= vt*tx*(fac+(k-1)*(1.-fac)/4.0)
c     +                 + corr_fac*vy
c           end do
c JAMES v2

         enddo
       enddo
c$omp end parallel do

       return
       end
c
       subroutine H10_slp(slp,IL,JL,ds,vm,rm,txc,tyc,pc,pe)

       implicit none
       integer i,j,IL,JL
       real slp_tmp,ds,tx,ty,rr,vm,rm,txc,tyc,fac,vm_sfc_sym
       real ee, pe, pc, bs, b, term1
       real slp(IL,JL)

c      Using eqn 2 in Holland 2010:                                                                                                                                                  
c      slp = pc + (pe - pc)*exp(-1.*(rm/r)**b)
c      ee = base of natural logarithms                                                                                                                                          
c      sfc density = 1.13 (at 980mb)  
       ee = 2.718281
       vm_sfc_sym = vm

       bs = (vm_sfc_sym*vm_sfc_sym*1.13*ee)/(100.*(pe-pc))

c      fac is a reduction factor from gradient wind to 10-m wind                                                                                                                      
       fac = 0.71

c      calculate b from bs from Holland et al. (2010)
       b = bs*(1/(fac*fac))

cc      find TC center and set pc
c       ic = int(txc)
c       jc = int(tyc)
c       slp(ic, jc) = pc


c$omp parallel do default(shared)                                                                                                                                                     
c$omp& private(j,i,tx,ty,rr,term1,slp_tmp)                                                                                                                             
      do j=1,JL
         ty = (real(j)-tyc)*ds
         do i=1,IL
           tx = (real(i)-txc)*ds
c           rr = max(1.0,sqrt(ty*ty+tx*tx))
           rr = sqrt(ty*ty+tx*tx)

           term1 = (rm/rr)**b

           slp(i,j) = pc + (pe - pc)*exp(-1.*term1)

         enddo
       enddo
c$omp end parallel do                                                                                                                                                              
       return
       end

c
      subroutine prxy_H10(u, v, prx,pry,f_2d,IL,JL,km, ds,vm,rm,r34,
     +                     txc, tyc, pc, pe, ux, vy)
       implicit none
       integer i,j,k,IL,JL,km
       real tx,ty,rr,vt,vm,rm,ds,txc,tyc,fac,ux,vy,vm_sfc_sym, spd_tc
       real ee, pe, pc, r34, bs, term1, xx_it, xxterm, term1_34
       real prx(IL,JL), pry(IL,JL),f_2d(IL,JL),u(IL,JL,km), v(IL,JL, km)
       real corr_fac, correction_vt

c      pe = environment pressure value is set to climatology
c      ee = base of natural logarithms
c      sfc density = 1.13 (at 980mb)
c      eqn 7 in Holland 2010. 
c       pe = 1010.
c      fac is a reduction factor from gradient winf to 10-m wind 
c      vm_sfc_sym = vm - 0.5*sqrt(ux*ux+vy*vy)
       vm_sfc_sym = vm
       ee = 2.718281
       bs = (vm_sfc_sym*vm_sfc_sym*1.13*ee)/(100.*(pe-pc))
       term1_34 = (rm/r34)**bs

c JAMES v2
c Change the reduction factor to 0.76
c       fac = 0.85
c Change the reduction factor to 0.71
c       fac = 0.76
       fac = 0.71
c JAMES v2



c  calculate xx_it in eq. 10  at R34, v=17
c  when rr = r34, xxterm = xx_it
      xx_it = 0.5
 10   vt = vm_sfc_sym*(term1_34*exp(1 - term1_34))**xx_it
      if(vt.gt.16.and.xx_it.le.2.0) then
         xx_it = xx_it + 0.02
         go to 10
      end if
c
      spd_tc = sqrt(ux*ux+vy*vy)
c
c$omp parallel do default(shared)
c$omp& private(j,i,k,tx,ty,rr,term1,xxterm,vt,corr_fac )
      do j=1,JL
         ty = (real(j)-tyc)*ds
         do i=1,IL
           tx = (real(i)-txc)*ds
           rr = max(1.0,sqrt(ty*ty+tx*tx))
           term1 = (rm/rr)**bs
           if(rr.le.rm) then
             xxterm = 0.5
           else
             xxterm = 0.5 + (rr - rm)*(xx_it-0.5)/(r34 - rm)
           end if
           vt = (vm_sfc_sym/rr)*(term1*exp(1 - term1))**xxterm/fac
           corr_fac = 0.
c           corr_fac = 2.0*(rm*rr)/(rm**2 + rr**2)
c     +               *1.173*spd_tc**0.63/spd_tc
           prx(i,j)=(f_2d(i,j)+vt)*vt*tx +corr_fac*f_2d(i,j)*vy
           pry(i,j)=(f_2d(i,j)+vt)*vt*ty -corr_fac*f_2d(i,j)*ux
         enddo
       enddo
c$omp end parallel do

       return
       end
c
       subroutine vortex(u,v,prx,pry,f_2d,IL,JL,km,ds,vm,rm,ux,vy,
     +                   txc,tyc)

       implicit none
       integer i,j,k,IL,JL,km
       real b,r0,rd,tx,ty,r,rr,vt,vm,rm,ux,vy,ds,txc,tyc,fac
       real u(IL,JL,km),v(IL,JL,km),prx(IL,JL),pry(IL,JL),f_2d(IL,JL)

       b=0.6
       r0=1.0e6
       fac=0.85

       do j=1,JL
       ty=(real(j)-tyc)*ds
       do i=1,IL
       tx=(real(i)-txc)*ds
       rr=sqrt(ty*ty+tx*tx)
       rd=rr/rm
       vt=vm/rm*exp((1.d0-rd**b)/b) 
       if(rd.lt.1.0) vt=vt*0.2*(1.0+4.0*sqrt(rd))
       if(rr.ge.r0) vt=vt*exp(-2.0*(rr-r0)/r0)
       vt=vt/fac
       prx(i,j)=(f_2d(i,j)+vt)*vt*tx+f_2d(i,j)*vy
       pry(i,j)=(f_2d(i,j)+vt)*vt*ty-f_2d(i,j)*ux
       do k=1,km
       u(i,j,k)=-vt*ty+ux
       v(i,j,k)=+vt*tx+vy
       enddo
       enddo
       enddo

       return
       end

c  to get estimated symmetric surface-level tangential winds from gradient-levels winds, 
c  uses a method from Knaff et al. (2011).
c  This method uses a reduction factor of 0.90 up to a radius(rr0) of 100 km
c  a reduction factor of 0.75 for any radius 700 km or greater 
c  and a linear decreasing reduction factor for any radius between those two radius values 
c  If the point is over land (true for any county), this reduction factor is further reduced by 20%
       real function wind_sfc_sym(spd_gl, rr_dist, is_land)
         
       implicit none
       integer is_land
       real reduction_fac, reduction_fac0, reduction_fac1
       real rr0, rr1, rr_dist, spd_gl

       rr0 = 100.0e3
       rr1 = 700.0e3
       reduction_fac0 = 0.9
       reduction_fac1 = 0.75

       if (rr_dist.le.rr0) then
         reduction_fac = reduction_fac0
       else if(rr_dist.ge.rr1) then
         reduction_fac = reduction_fac1
       else
         reduction_fac = reduction_fac0-
     +         (rr_dist-rr0)*((reduction_fac0-reduction_fac1)/(rr1-rr0))
       end if

       if(is_land.eq.1) then
         reduction_fac = reduction_fac*0.8
       end if

       wind_sfc_sym = spd_gl*reduction_fac

       end function wind_sfc_sym

c under construction
c calculate gradient wind from symmetric surface wind
      real function spd_gl(vmax, is_land)
      real vmax, reduction_factor
      integer is_land

      reduction_factor = 0.9
      if(is_land.eq.1) then
         reduction_factor = reduction_factor*.8
      end if
      spd_gl = vmax/reduction_factor
  
      end
c Calculate the direction of gradient winds at each location
      real function dir_gl(u,v)
      
      implicit none
      real, parameter :: pi = 3.14159265 
      real u,v
    
      if (abs(v) > .0001) then
         dir_gl =atan(u/v)
      else
         if (u > 0.) then
            dir_gl = -pi/2.0
         else
            dir_gl = pi/2.0
         endif
      endif
      dir_gl = dir_gl*360./(2.0*pi)
      if (v   > 0.) dir_gl = dir_gl + 180.
      if (dir_gl < 0.) dir_gl = dir_gl + 360. 
      end function dir_gl

c Calculate the surface wind direction from gradient winds
c the function adds an inflow angle to the gradient wind direction.
c because surface friction changes the wind direction 
c The inflow angle is calculated as a function of the distance from the storm center 
c to the location and the storm's Rmax at that observation point (eq. 11a--c, Phadke et al. 2003):

      real function dir_sfc(dir_gl, rr_dist, is_land, rmax)
      implicit none
      real dir_gl, rr_dist, is_land, rmax,  inflow_angle 

c     Calculate inflow angle over water based on radius of location from storm
c     center in comparison to radius of maximum winds (Phadke et al. 2003)
      if(rr_dist.lt.rmax) then      
         inflow_angle = 10 + (1 + (rr_dist/rmax))
      else if(rr_dist.ge.rmax.and.rr_dist.le.1.2*rmax) then
         inflow_angle = 20 + 25 * ((rr_dist/rmax) - 1)
      else
        inflow_angle = 25 
      end if

c Add 20 degrees to inflow angle since location is over land, not water
      if(is_land.eq.1) then
        inflow_angle = inflow_angle + 20
      end if      
      
      dir_sfc = dir_gl + inflow_angle
      if (dir_sfc > 360.) dir_sfc = dir_sfc - 360.
      end function dir_sfc

c add back in wind from forward speed of TC
      subroutine  wind_sfc(u_sfc,v_sfc,spd_sfc_sym,dir_sfc,rmax,rr_dist
     +            ,tcspd_u, tcspd_v)
      implicit none
      real*4 u_sfc,v_sfc
      real deg_2_rad, spd_sfc_sym, dir_sfc
      real correction_factor, rmax, rr_dist, tcspd_u, tcspd_v
 
c     3.14/180      
      deg_2_rad = 0.01745239 
      u_sfc = -spd_sfc_sym*sin(dir_sfc*deg_2_rad)
      v_sfc = -spd_sfc_sym*cos(dir_sfc*deg_2_rad)
 
c     Add back in component from forward motion of the storm
      correction_factor = (rmax * rr_dist)/(rmax*rmax + rr_dist*rr_dist)
      
c     Add tangential and forward speed components and calculate
c     magnitude of this total wind
      u_sfc = u_sfc + correction_factor * tcspd_u      
      v_sfc = v_sfc + correction_factor * tcspd_v  
      
      end subroutine wind_sfc 

c calculate all the inputs needed for the Willoughby wind profile   
c equation 10 a-c
c lat_tc: Latitude, in decimal degrees 
      subroutine get_wind_param_eq10(X1, NN, AA, vm_gl, lat_tc)

c JAMES v2
c Convert latitude to absolute values of latitude so the parameters are the same in each hemisphere

       implicit none
       real X1, NN, AA, vm_gl, lat_tc

       ! exponential decay length(km) in outer vortex (equation 10a)
       X1 = 317.1 - 2.026 * vm_gl + 1.915 * abs(lat_tc)
c JAMES test
       ! Irma Test to get better outer winds.
       ! X1 = 380. is good for control simulation of Irma
c       X1 = 380
       ! Adding 10% according to climate change
c       X1 = 418.
c JAMES test
       ! km ->m  
       X1 = X1*1.0e3
       ! Power-law exponent for porfile inside the eye (equation 10b)
       NN = 0.4067 + 0.0144 * vm_gl - 0.0038 * abs(lat_tc)
       ! Fraction of outer profile represented  (equation 10c)
c  JAMES test
       AA = 0.0696 + 0.0049*vm_gl - 0.0064 * abs(lat_tc)
c       AA = 0.05
c JAMES test
       if(AA.lt.0) AA = 0.


       Print*, "X1 = ", X1
c       Print*, "NN = ", NN
       Print*, "AA = ", AA
c JAMES v2

       return
       end

c calculate all the inputs needed for the Willoughby wind profile   
c equation 11 a-c
c lat_tc: Latitude, in decimal degrees 
c rm: radius of maximum wind. 
      subroutine get_wind_param_eq11(X1, NN, AA, vm_gl, lat_tc, rm)

c JAMES v2
c Convert latitude to absolute values of latitude so the parameters are the same in each hemisphere

       implicit none
       real X1, NN, AA, vm_gl, lat_tc, rm

       ! exponential decay length(km) in outer vortex (equation 10a)
       X1 = 287.6 - 1.942 * vm_gl + 7.799*log(rm)+ 1.819 * abs(lat_tc)
       ! km ->m  
       X1 = X1*1.0e3
       ! Power-law exponent for porfile inside the eye (equation 10b)
       NN = 2.1340 + 0.0077* vm_gl - 0.4522*log(rm) -0.0038*abs(lat_tc)
       ! Fraction of outer profile represented  (equation 10c)
       AA = 0.5913 + 0.0029*vm_gl -0.1361*log(rm)- 0.0042 * abs(lat_tc)
       if(AA.lt.0) AA = 0.

c       Print*, "X1 = ", X1
c       Print*, "NN = ", NN
c       Print*, "AA = ", AA
c JAMES v2

       return
       end

c calculate lower(R1) and upper(R2) boundary of the transition zone (km)
c assumes that  R2 − R1 (the width of the transition region) is 25 kilometers 
c when Rmax is larger than 20 kilometers and 15 kilometers otherwise.
      subroutine get_R1R2(R1, R2, rm, X1, NN, AA)

      implicit none
      real R1, R2, WID, rm, NN, AA, X1, X2, xx, wt

      X2 = 25.e3
      ! Compute weight at radius of maximum wind  dual exponentials
      wt = 1.0/(1.0 + (rm/NN)*((1.0-AA)/X1 + AA/X2))

      ! Find nondimensional coordinate of that weight
      ! by interval bisection.
       call findx(wt, xx)

      if (rm.gt.20.e3) then
        WID = 25.e3
      else
        WID = 15.e3
      end if

      !Approximate nondimensional coordinate using asyptotic cubic parabola
      ! x2 = 1.0 - (0.1*(1-w))**0.333333
      !Compute inner boundary
      R1 = rm - xx*WID

      !Compute outer boundary
      R2 = R1 + WID
      
      return
      end subroutine get_R1R2

c Polynomial ramp function
c when x, the nondimensional argument, is  <= zero, wgt=0
c when x >= 1 wgt =1, 
c First and second derivitives are zero at both ends, C2 continuity.
      real function wgt(x)

      implicit none
      real x

      if (x .le. 0)then
        wgt = 0.0                !Zero before start of ramp
      else if (x .lt. 1.0) then  !Plynomial ramp form zero to one
        wgt = x**5*(126.0-x*(420.0-x*(540.0-x*(315.0-70.0*x)))) !Forth
      else                       !One after end of ramp
        wgt = 1.0
      endif

      end function wgt
c Uses interval bisection to find coordinate XX of a given
c weighting WW on the polynomial ramp function
c Returns: Nondimensional coordinate XX=[(r-R1)/(R2-R1)] where
c the polynomial ramp function, wgt, is equal to WW
c This subroutine is adapted from Hugh Willoughby
      subroutine findx(WW,XX)
      implicit none
      real eps, WW, XX, X1, X2, X3, wgt, W1, W2, W3
      parameter(eps=1.0e-4)  !Convergence criterion
      logical done           !Logical flag for convergence

       !Initialize interval
       done = .FALSE.
       X1=0.0             !Lower bound of interval
       W1 = wgt(X1) !Value there
       X2=1.0             !Upper bound of interval
       W2 = wgt(X2) !Value there

        !Iterate until "done"
       do while (.not.done)
        !Find value at midpoint of interval
         X3=0.5*(X1+X2)
         W3 = wgt(X3)

        !Which half contains zero crossing?
         if(W3.lt.WW)then  !Upper interval
           X1 = X3               !Discard lower
           W1 = W3
         else                            !Lower interval
           X2 = X3         !Dicard upper
           W2 = W3
         endif
         XX=X3             !Return midpoint

        !Converged yet?
         done = abs(X1-X2).lt.eps
       enddo               !End of iteration
       return
       end subroutine findx

      real function correction_vt(vt_tc) result(val)
      ! 2.0 for Vt=2kts, 1.0 for Vt=10kts  0.5 for Vt=20kts
      implicit none
      real vt_tc

      if(vt_tc.lt.1) then
        val = 2.0
      else if(vt_tc.lt.5) then
        val = 1.0 + 0.25*(5.0-vt_tc)
      else if(vt_tc.ge.5.and.vt_tc.lt.10) then
        val = 1.0 - 0.1*(vt_tc-5.0)
      else
         val = 0.5
      end if

      end function  correction_vt

c       calculates Rmax using Eq. 7a from Willoughby et al. (2006):
c        if(rm_1d(k).ge.800.0e3.or.rm_1d(k).le.1.0e3) then
c          rm_1d(k) = 46.4 * exp(-0.0155 * vm_1d(k) + 0.0169 * lon_1d(k))
c        end if
c --------------------------------------------
      subroutine willouby(u,v,f_2d,IL,JL,km,ds,vm,rm,txc,tyc,ux,vy,lat
     +                     ,land)

      implicit none
      integer i,j,k,IL,JL,km, land(IL,JL)
      real rr,vt,vm,rm,ds,txc,tyc,ux,vy, lat
      real tx, ty, VofR2, correction_factor
      real vm_gl, vm_sfc_sym,reduction_factor, reduction_land
      real X1, X2, NN, AA, R1, R2, WID,fac
      real spd_tc, correction_vt
      real u(IL,JL,km), v(IL,JL,km),f_2d(IL,JL)

      ! adjust Vmax by removing forward motion of storm
      ! (Phadke et al. 2003)
      ! Vmax_sfc_sym = vmax - 0.5*TC_speed
      !vm_sfc_sym = vm - 0.5*sqrt(ux*ux+vy*vy) 
      vm_sfc_sym = vm

c JAMES v2
c Changed the reduction factor to 0.76
c      reduction_factor = 0.85
c Changed the reduction factor to 0.71
c      reduction_factor = 0.76
      reduction_factor = 0.71
c JAMES v2


      ! origional reduction_land   = 0.8
      ! Dec. 3 2018. 
      !The reduction_land factor in Willoughby is causing us some problems with a strange speed up of winds near the coast.
      reduction_land   = 1.0 

      vm_gl = vm_sfc_sym/reduction_factor

      X2 = 25.e3

      if (rm.gt.20.e3) then
        WID = 25.e3
      else
        WID = 15.e3
      end if

      call get_wind_param_eq10(X1, NN, AA, vm_gl, lat)
      !call get_wind_param_eq11(X1, NN, AA, vm_gl, lat, rm*0.001)
 
      call get_R1R2(R1, R2, rm, X1, NN, AA)

      spd_tc = sqrt(ux*ux+vy*vy)
c      
c$omp parallel do default(shared)
c$omp& private(j,i,k,tx,ty,rr,vt, fac, correction_factor)
      ! two-exponential profile
      do j=1,JL
        ty = (real(j)-tyc)*ds
        do i=1,IL
          tx = (real(i)-txc)*ds
          rr = max(1.0,sqrt(ty*ty+tx*tx))
          !Compute winds 
          if(land(i,j).eq.1) then
            fac = reduction_factor*reduction_land
          else
            fac = reduction_factor
          end if

c JAMES v2
c Remove the double counting of the reduction factor 'fac'. 
c The reduction factor was applied earlier in this subroutine.
c          vt = VofR2(rr,vm_gl, rm, R1,R2, WID, NN, X1, X2, AA)/fac
          vt = VofR2(rr,vm_gl, rm, R1,R2, WID, NN, X1, X2, AA)
c JAMES v2


c JAMES v2
c Decided to add the asymmetry i the post-processing 
          correction_factor = 0.0
c          correction_factor = 2.0*(rm*rr)/(rm**2 + rr**2)
c     +                        *1.173*spd_tc**0.63/spd_tc
c JAMES v2

c JAMES v2
c Added a correction to get the correct sign of the U and V winds in each hemisphere

          if(f_2d(i,j).le.0) then
            do k=1,km
              u(i,j,k)= vt*ty/rr + correction_factor*ux
              v(i,j,k)=-vt*tx/rr + correction_factor*vy
            enddo
          else
            do k=1,km
              u(i,j,k)=-vt*ty/rr + correction_factor*ux
              v(i,j,k)= vt*tx/rr + correction_factor*vy
            enddo
          endif
c JAMES v2

c JAMES v2
c Remove this different treatment of the lowest levels.
c           do k=1,4
c             u(i,j,k)=-vt*ty*(fac+(k-1)*(1.-fac)/4.0)
c     +                 /rr + correction_factor*ux
c             v(i,j,k)= vt*tx*(fac+(k-1)*(1.-fac)/4.0)
c     +                 /rr + correction_factor*vy
c           end do
c JAMES v2
        end do
      end do
c$omp end parallel do

      return
      end subroutine willouby

c --------------------------------------------
      subroutine prxy_willouby(prx,pry,f_2d,IL,JL,km,ds,vm,rm,txc,tyc,
     +                     ux,vy,lat,land)

      implicit none
      integer i,j,k,IL,JL,km, land(IL,JL)
      real rr,vt,vm,rm,ds,txc,tyc,ux,vy, lat
      real tx, ty, VofR2,fac, correction_factor
      real vm_gl, vm_sfc_sym,reduction_factor, reduction_land
      real X1, X2, NN, AA, R1, R2, WID
      real spd_tc, correction_vt
      real f_2d(IL,JL), prx(IL,JL),pry(IL,JL)

      ! adjust Vmax by removing forward motion of storm
      ! (Phadke et al. 2003)
      ! Vmax_sfc_sym = vmax - 0.5*TC_speed
      !vm_sfc_sym = vm - 0.5*sqrt(ux*ux+vy*vy)
      vm_sfc_sym = vm


c JAMES v2
c Changed the reduction factor to 0.76
c      reduction_factor = 0.85
c Changed the reduction factor to 0.71
c      reduction_factor = 0.76
      reduction_factor = 0.71
c Justified in KW01 showed it ranges from .9 to .7 around the storm, 
c and that is the factor from gradient to 22 m. 
c We are going to 10 m so it'll be even more justified at the .7 level.  
c Vickerey shows a factor of about 0.72 going from the supergradient jet to sfc. 
c We don't allow enough time for jet to form so we artifically put it in.   
c In anycase, 0.76 is well within 1sd of the mean (0.81) of 150 dropsode obs of Powell.
c JAMES v2


      ! original reduction_land   = 0.8
      reduction_land   = 1.0

      vm_gl = vm_sfc_sym/reduction_factor

      X2 = 25.e3

      if (rm.gt.20.e3) then
        WID = 25.e3
      else
        WID = 15.e3
      end if

      call get_wind_param_eq10(X1, NN, AA, vm_gl, lat)
      !call get_wind_param_eq11(X1, NN, AA, vm_gl, lat, rm*0.001)
      call get_R1R2(R1, R2, rm, X1, NN, AA)

      spd_tc = sqrt(ux*ux+vy*vy)
 
c$omp parallel do default(shared)
c$omp& private(j,i,k,tx,ty,rr,vt,fac,correction_factor)
      ! two-exponential profile
      do j=1,JL
        ty = (real(j)-tyc)*ds
        do i=1,IL
          tx = (real(i)-txc)*ds
          rr = max(1.0,sqrt(ty*ty+tx*tx))
          !Compute winds 
          if(land(i,j).eq.1) then
            fac = reduction_factor*reduction_land
          else
            fac = reduction_factor
          end if

c JAMES v2
c Remove the double counting of the reduction factor 'fac'
c The reduction factor was applied earlier in this subroutine.
c          vt = VofR2(rr,vm_gl, rm, R1,R2, WID, NN, X1, X2, AA)/fac
          vt = VofR2(rr,vm_gl, rm, R1,R2, WID, NN, X1, X2, AA)
c JAMES v2

c JAMES v2
c Decided to add the asymmetry in the post-processing
          correction_factor = 0.0
c          correction_factor = 2.0*(rm*rr)/(rm**2 + rr**2) 
c     +                       *1.173*spd_tc**0.63/spd_tc
c JAMES v2
          prx(i,j)=(f_2d(i,j)+vt/rr)*vt*tx/rr
     +             +correction_factor*f_2d(i,j)*vy
          pry(i,j)=(f_2d(i,j)+vt/rr)*vt*ty/rr
     +             -correction_factor*f_2d(i,j)*ux

        end do
      end do
c$omp end parallel do
      return
      end subroutine prxy_willouby   

c --------------------------------------------
c Sectionally continuous wind-profile function with a
c power law inside the eye and a double-length exponential
c decay outside. There is a smooth polynimial transition
c across the radius of maximum wind.
c Returns: Double exponent wind profile as a function of radius
c Adapted from Hugh Willoughby
      real function VofR2(RR,Vmax,Rmax,R1,R2,WID,ENI,XL1,XL2,AA)

      implicit none
      real RR,Vmax,Rmax,R1,R2,WID,ENI,XL1,XL2,AA, Vin, Vout, w, wgt

      ! Inside the eye, use power law
      if(RR.le.R1) then
         VofR2 = Vin(RR,Vmax,Rmax,ENI)

      !In transition, Compute weighted combination of inside and outside profiles
      else if(RR.le.R2) then
        w = wgt((RR-R1)/WID)
        VofR2 =(1.0-w)*          Vin(RR,Vmax,Rmax,ENI)
     1            + w*((1.0-AA)*Vout(RR,Vmax,Rmax,XL1)
     2                 + AA *Vout(RR,Vmax,Rmax,XL2))
      !Outisde eye and transition,Compute exponential profiles
      else
        VofR2 = (1.0-AA)*Vout(RR,Vmax,Rmax,XL1)
     1             + AA *Vout(RR,Vmax,Rmax,XL2)
        endif
      return
      end function VofR2

c --------------------------------------------
c Sectionally continuous Willoughby wind-profile function with a
c power law inside the eye and a hyperbolic cosecant
c decay outside. There is a smooth polynimial transition
c across the radius of maximum wind.
c Returns: Swirling wind as a function or radius
      real function VofRH(RR,Vmax,Rmax,R1,R2,WID,ENI,XL1,AA)
      implicit none
      real RR,Vmax,Rmax,R1,R2,WID,ENI,XL1,AA,w, Vin, VH, wgt

      if(RR.le.R1) then
        ! Inside the eye, use power law 
        VofRH = Vin(RR,Vmax,Rmax,ENI)
      else if(RR.le.R2) then
        !In transition area, 
        w = wgt((RR-R1)/WID)

        !Compute weighted combination of inside and outside profiles
        VofRH = (1.0-w)*Vin(RR,Vmax,Rmax,ENI)
     1           + w* VH(RR,Vmax,Rmax,XL1,AA)
      else
        !Outisde eye and transition,Compute hyperbolic profile
        VofRH = VH(RR,Vmax,Rmax,XL1,AA)
      endif
      return
      end function VofRH

c --------------------------------------------
c Power-law wind profile for inside the eye
      real function Vin(RR,Vmax,Rmax,ENI)
      implicit none
      real RR,Vmax,Rmax,ENI

      Vin = Vmax*(RR/Rmax)**ENI
      return
      end function Vin

c --------------------------------------------
c Single-length exponential for outisde the eye
      real function Vout(RR,Vmax,Rmax,Xlen)
      implicit none
      real RR,Vmax,Rmax,Xlen

      Vout = Vmax*exp(-(RR - Rmax)/Xlen)
      return
      end function Vout

c --------------------------------------------
c Hyperbolic cosecant outer profile
      real function VH(RR,Vmax,Rmax,Xlen,AA)
      implicit none
      real RR,Vmax,Rmax,Xlen,AA, x

      x = exp((RR-Rmax)/Xlen)
      VH = Vmax*(1.0-AA)/(x - AA/x)
      return
      end function VH
c --------------------------------------------
c       calculates Rmax using Eq. 7a from Willoughby et al. (2006):
c        if(rm_1d(k).ge.800.0e3.or.rm_1d(k).le.1.0e3) then
c          rm_1d(k) = 46.4 * exp(-0.0155 * vm_1d(k) + 0.0169 * lon_1d(k))
c        end if


