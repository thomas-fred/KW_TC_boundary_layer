;==================================================
; select, process, and plot all TC tracks within 200km of land
; input files:  rmax_ib.nc, vmax_ib.nc, IBTrACS.ALL.v04r00.nc
;	        geo.nc
; output files: bdy_10min.txt, lat_2d.dat, topo.dat, landuse.dat
; Ming Ge, March 2020 and James Done March 2021
;==================================================

begin

  region_s = "PR"

  ; do we want to write out the domain data files needed for the footprint model ( landuse.dat  lat_2d.dat  topo.dat )?
  is_write = True
  ; do we want to write out the domain data files just once for only the first storm encountered and then stop the script?
  is_dat   = False

  ; set working directories
  dir_in = "/glade/work/done/footprints/raincoat/input_data/"
  dir_code = "/glade/work/done/footprints/raincoat/code/"
  dir_out = "/glade/work/done/footprints/raincoat/output_data/"

  ; set missing values
  missing_i = -99
  missing_r = -99.

  ; set constants
  delim        = " "
  km_2_nm      = 0.539956803         ; 1 kilometer = 0.539956803 nm
  min10_per_hr = 6
  kt2mps       = 0.514
  nm2m         = 1852

  ; read in the domain geo file
  flnm = dir_in + "geo_" + str_lower(region_s) + ".nc"

  ; read in domain data
  print("read in:" + flnm)
  ff = addfile(flnm, "r")
  lat_2d = ff->XLAT_M(0,:,:)
  lon_2d = ff->XLONG_M(0,:,:)
  topo   = ff->HGT_M(0,:,:)
  landuse = toint(ff->LU_INDEX(0,:,:))

  print("Setting TOPO to zero for ocean points")
  topo = where(landuse.eq.17,0.,topo)

  ; Perform some smoothing of terrain for numerical stability
  print("TOPO before smooth:" + max(topo) + " the maxmum will be 2000")
  if(max(topo).gt.2000.)then
    topo = 2000.0/max(topo)*topo
  end if

  ; used to cut TC tracks close to domain boundaries
  dis_bdy_s = 0.3         ; deg
  dis_bdy_n = 0.3
  dis_bdy_e = 0.3
  dis_bdy_w = 0.3

  lat_s = min(lat_2d)
  lat_n = max(lat_2d)
  lon_w = min(lon_2d)
  lon_e = max(lon_2d)

  dims  = dimsizes(lat_2d)
  n_lat = dims(0)
  n_lon = dims(1)
  delete(dims)

  print("For wind field boundary layer model set-up:  n_lat(lp)=" + n_lat + " n_lon(lq)=" + n_lon)

  yct = tointeger(0.5*(n_lat + 1)) - 1
  xct = tointeger(0.5*(n_lon + 1)) - 1
  print("center point xct=" + xct + " yct=" +  yct)

  ; Read in historical track data from IBTrACS
  flnm_ib = dir_in + "rmax_ib.nc"
  f_ib    = addfile(flnm_ib,   "r")
  rmax_ib = f_ib->rmax_ib
  print("---- read in: " + flnm_ib )
  print(dimsizes(rmax_ib))

  flnm_ib = dir_in + "vmax_ib.nc"
  f_ib    = addfile(flnm_ib,   "r")
  spd_ib = f_ib->vmax_ib	;(storm, date_time)
  print("---- read in: " + flnm_ib )
  print(dimsizes(spd_ib))

  flnm_ib = dir_in + "IBTrACS.ALL.v04r00.nc"
  print("---- read in: " + flnm_ib )

  f_ib     = addfile(flnm_ib,   "r")
  basin    = tostring(f_ib->basin)
  time     = f_ib->time	 ;v4  double time(storm, date_time) days since 1858-11-17 00:00:00
  print(dimsizes(time))

  name_ib  = tostring(f_ib->name)
  dist2land = toint(f_ib->landfall)
  lat_ib   = f_ib->lat   ; (storm, date_time)
  lon_ib   = f_ib->lon
  pe_ib    = f_ib->bom_poci
  r34_ib   = dim_avg_n(f_ib->usa_r34, 2)
  pmin_ib  = f_ib->usa_pres
  vt_ib    = short2flt(f_ib->storm_speed)  ; short(storm, date_time) kts
  dir_ib   = short2flt(f_ib->storm_dir)  ; (storm, date_time)

; Prevent a 360 rotation when the direction changes from 359 to 1 degree in one timestep.
  dir_ib   = where(dir_ib.lt.180., dir_ib+360., dir_ib)
  spd_ib   = where(spd_ib.gt.0., spd_ib, spd_ib@_FillValue)

  dims     = dimsizes(lat_ib)
  n_tc_ib  = dims(0)
  n_date_time = dims(1)
  delete(dims)

  utc_date  = ut_calendar(time, 0)
  year_ib   = toint(utc_date(:, :, 0))
  mm_ib     = toint(utc_date(:, :, 1))
  dd_ib     = toint(utc_date(:, :, 2))
  hh_ib     = toint(utc_date(:, :, 3))
  mmddhh_ib = toint(mm_ib*10000 + dd_ib*100 + hh_ib)

  juli_ib   = greg2jul(year_ib, mm_ib, dd_ib, hh_ib)

  ; NOTE: for different regions, make sure longitudes are correct.
  if(region_s.eq."GULF") then   ; negative lon
    lon_ib = where(lon_ib.gt.180, lon_ib-360, lon_ib)
  end if

  ; Start to look for all histroical TC tracks the cross the geographic domain

  ; note: variables without missing value: time, vt_ib, dir_ib
  ;       last one is missing: dist2land
  ; there are missing lat_ib/lon_ib in the middle of obs.
  ; cut data to model domain
  lat_ib = where(lat_ib.lt.lat_s+dis_bdy_s.or.lat_ib.gt.lat_n-dis_bdy_n.or.	\
 		 lon_ib.lt.lon_w+dis_bdy_w.or.lon_ib.gt.lon_e-dis_bdy_e,	\
                 lat_ib@_FillValue, lat_ib)
  lon_ib = where(.not.ismissing(lat_ib), lon_ib, lon_ib@_FillValue)

  ; ##########################
  ; for checking data change False to True here.
  if(False) then
  ct = 0
  do nt = 0, n_tc_ib - 1
    if(.not.all(ismissing(lat_ib(nt,:)))) then
      ind_all = ind(.not.ismissing(lat_ib(nt,:)))
      ind_s = min(ind_all)
      ind_e = max(ind_all)
      ; change variable
      do ni = ind_s, ind_e
        if(ismissing(lat_ib(nt,ni))) then
          print(year_ib(nt,ni)+""+mmddhh_ib(nt,ni)+" "+ dir_ib(nt,ni)+" "+ vt_ib(nt,ni) + " " + name_ib(nt)+" "+ni+"-"+ind_e )

          ;print(year_ib(nt,ind_s:ind_e) + "" + mmddhh_ib(nt,ind_s:ind_e) + " " + dist2land(nt,ind_s:ind_e) + " " + name_ib(nt))
          ct = ct + 1
        end if
      end do
      delete([/ind_all/])
    end if
  end do
  print(ct + "")
  end if
  ; ####################
  year_ib   = where(lat_ib.lt.lat_s+dis_bdy_s.or.lat_ib.gt.lat_n-dis_bdy_n.or. \
                    lon_ib.lt.lon_w+dis_bdy_w.or.lon_ib.gt.lon_e-dis_bdy_e,    \
                    year_ib@_FillValue, year_ib)

  ; Optional: Only keep tracks from a specified year
;  year_ib   = where(year_ib.eq.2007,year_ib,year_ib@_FillValue)

  mmddhh_ib = where(.not.ismissing(year_ib), mmddhh_ib, mmddhh_ib@_FillValue)
  rmax_ib   = where(.not.ismissing(year_ib(0:13359,:)), rmax_ib,   rmax_ib@_FillValue)
  spd_ib    = where(.not.ismissing(year_ib(0:13359,:)), spd_ib,    spd_ib@_FillValue)
  vt_ib     = where(.not.ismissing(year_ib), vt_ib,     vt_ib@_FillValue)
  dir_ib    = where(.not.ismissing(year_ib), dir_ib,    dir_ib@_FillValue)
  dist2land = where(.not.ismissing(year_ib), dist2land, dist2land@_FillValue)

  print("----IB data info ----")
  print("lat_track : "   + min(lat_ib) + "/" + max(lat_ib)  + "   lon: " + min(lon_ib)  + "/" + max(lon_ib))
  print("lat_domain: "   + min(lat_s) + "/" + max(lat_n)  + "   lon: " + min(lon_w)  + "/" + max(lon_e))
  print("vmax: "  + min(spd_ib) + "/" + max(spd_ib))
  print("dist2land(km): "  + min(dist2land) + "/" + max(dist2land))
  print("vt(kt): "  + min(vt_ib) + "/" + max(vt_ib))
  print("month: " + min(mmddhh_ib) + "/" + max(mmddhh_ib) + "year: "  + min(year_ib) + "/" + max(year_ib))

  delete([/time, utc_date, mm_ib, dd_ib, hh_ib/])

  ; plot a map of all the identified TC tracks
  flnm = dir_out + "tracks_landfall_50"
  wks  = gsn_open_wks("png", flnm )
  gsn_define_colormap(wks, "MPL_Greys")

  res                        = True
  res@gsnDraw                = False
  res@gsnFrame               = False
  res@gsnMaximize            = True
  res@mpOutlineBoundarySets  = "GeophysicalAndUSStates"

  rest                       = True
  rest@gsLineThicknessF      = 1.5
  rest@gsLineColor           = "Black"

  first                      = True
  first@gsMarkerIndex        = 1
  first@gsMarkerColor        = "black"
  first@gsMarkerThicknessF   = 6.0
  first@gsMarkerSizeF	     = 0.015

  res@mpMinLatF = lat_s
  res@mpMaxLatF = lat_n
  res@mpMinLonF = lon_w
  res@mpMaxLonF = lon_e
  res@mpCenterLonF = 0.5*(lon_w + lon_e)

  txres = True
  txres@txFontHeightF = 0.007
  txres@txFontColor   = "black"
  txres@txAngleF      = 0

  map = gsn_csm_map_ce(wks, res)
  text_st = new(1000, graphic)
  dum     = new(1000, graphic)

  ; for interpolation
  opt = True
  opt@MAP_PROJ  = ff@MAP_PROJ
  opt@TRUELAT1  = ff@TRUELAT1
  opt@TRUELAT2  = ff@TRUELAT2
  opt@STAND_LON = ff@STAND_LON
  opt@DX        = ff@DX
  opt@DY        = ff@DY
  opt@REF_LAT = lat_2d(0,0)
  opt@REF_LON = lon_2d(0,0)
  opt@KNOWNI  = 1.0
  opt@KNOWNJ  = 1.0

  ; search and keep only tracks with valid data and that come within 200km of land
  ct = 0
  do ns = 0,  n_tc_ib - 1
     ; at least 2 valid rmax points & 200km closer to land
     if(any(dist2land(ns,:).le.200).and.dimsizes(ind(ismissing(rmax_ib(ns,:)))).lt.n_date_time - 1) then
      ind_all = ind(.not.ismissing(dist2land(ns,:)))
      ind_s = min(ind_all)
      ind_e = max(ind_all)
      delete(ind_all)

      ; get track starting point ind_0
      ind_0 = 999
      do nd = ind_s, ind_e
        ; NOTE: by doing following, it will miss the landing after the first
        ;       that is one TC one count
        ;if(.not.ismissing(rmax_ib(ns,nd)).and.dist2land(ns,nd).le.200 \
	;   .and.all(dist2land(ns,nd:ind_e).le.200)) then
        if(.not.ismissing(rmax_ib(ns,nd)).and.dist2land(ns,nd).le.1000) then
          ind_0 = nd
          break
        end if
      end do

      if((ind_0.eq.999).or.(region_s.eq."EP".and.lat_ib(ns,ind_0).ge.18.and.lon_ib(ns,ind_0).ge.-98) \
         .or.(region_s.eq."EP_E".and.lat_ib(ns,ind_0).ge.16)) then
        continue
      end if

      ; get track ending point ind_1
      ind_1 = 999
      do nd = ind_0 + 2, ind_e
       ind_c = ind_e+ind_0 - nd + 2
       dist2land_tmp = where(dist2land(ns,:).gt.1000, dist2land@_FillValue, dist2land(ns,:))
       if(.not.ismissing(dist2land_tmp(ind_c-1))) then
          ind_1 = ind_c - 1
          break
        end if
      end do

      if(ind_1.eq.999) then
        continue
      end if

      print("JAMES ind_1 is " + ind_1)

      n_time = ind_1 - ind_0 + 1
      print("JAMES lat_ib " + lat_ib(ns,ind_0::))

      lat_ib_1d  = lat_ib (ns,ind_0:ind_1)
      lon_ib_1d  = lon_ib (ns,ind_0:ind_1)
      juli_ib_1d = juli_ib(ns,ind_0:ind_1)

      ; monotonically increasing or decreasing
      do nn = 1, n_time - 1
        if(abs(juli_ib_1d(nn)-juli_ib_1d(nn-1)).le.0.01) then
           juli_ib_1d(nn) = juli_ib_1d(nn) + 0.01
        end if
      end do

      vt_ib_1d   = 1.173*(vt_ib(ns, ind_0:ind_1)*kt2mps)^0.63
      dir_ib_1d  = dir_ib (ns, ind_0:ind_1)
      vmax_ib_1d = spd_ib (ns, ind_0:ind_1)
      rmax_ib_1d = rmax_ib(ns, ind_0:ind_1)
      r34_ib_1d  = r34_ib (ns, ind_0:ind_1)
      pmin_ib_1d = pmin_ib(ns, ind_0:ind_1)
      pe_ib_1d   = pe_ib  (ns, ind_0:ind_1)
      pe_ib_1d   = where(ismissing(pe_ib_1d), toshort(1005.), pe_ib_1d)

      ; interpolate 6-hourly data to 10minute data
      n_int = toint((max(juli_ib_1d) - min(juli_ib_1d))*24*min10_per_hr) + 1

      xi      = juli_ib_1d
      xo      = fspan(min(xi), max(xi), n_int)
      xo_juli = fspan(min(xi), max(xi), (n_int-1)/min10_per_hr+1)
      greg    = jul2greg(xo_juli)

      mmddhh_ib_1d = greg(:,1)*10000 + greg(:,2)*100 + greg(:,3)

      mmddhhmm_1 = new(n_int, integer)
      mmddhhmm_s = new(n_int, string)
      do nn = 0, n_int-1
        mmddhhmm_1(nn) = mmddhh_ib_1d(nn/min10_per_hr)*100 + mod(nn,min10_per_hr)*10
      end do

      mmddhhmm_s = where(mmddhhmm_1.gt.10000000, tostring(mmddhhmm_1), "0"+tostring(mmddhhmm_1))

      ind_r = ind(rmax_ib_1d)

      lat_1 = linint1 (xi, lat_ib_1d,  False, xo, 0)
      lon_1 = linint1 (xi, lon_ib_1d,  False, xo, 0)
      vm_1  = linint1 (xi, vmax_ib_1d*kt2mps- vt_ib_1d, False, xo, 0)
      rm_1  = linint1 (xi(ind_r), rmax_ib_1d(ind_r), False, xo, 0)*nm2m
      r34_1 = linint1 (xi, r34_ib_1d,  False, xo, 0)*nm2m
      pc_1  = linint1 (xi, pmin_ib_1d, False, xo, 0)
      pe_1  = linint1 (xi, pe_ib_1d,   False, xo, 0)
      vt_1  = linint1 (xi, vt_ib_1d,   False, xo, 0)
      dir_1 = linint1 (xi, dir_ib_1d,  False, xo, 0)

; set missing rmax rows to last positive value
      non_missing_rm_1 = ind(.not.ismissing(rm_1))
      rm_1 = where(ismissing(rm_1),rm_1(max(non_missing_rm_1)),rm_1)

;; set missing r34 rows to last positive value
;      non_missing_r34_1 = ind(.not.ismissing(r34_1))
;      r34_1 = where(ismissing(r34_1),r34_1(max(non_missing_r34_1)),r34_1)

      loc = wrf_ll_to_ij(lon_1, lat_1, opt)
      typ = loc(0,:)
      txp = loc(1,:)

      yyyy_1 = mmddhhmm_1
      yyyy_1 = toint(year_ib(ns, ind_0))
      yyyy_1 = where(mmddhhmm_1.lt.mmddhhmm_1(0), yyyy_1+1, yyyy_1)

      name_1 = new(dimsizes(yyyy_1), string)
      name_1 = name_ib(ns)

      if(is_write) then
        ; write out track data
        line_s = new(n_int+1, string)
         line_s(0)  = " jt0(x lon)  it0(y lat)   vmax     pmin       rmax      r34       pe        lat       lon  storm speed storm dir yyyymmddhhmm name"
         line_s(1:) = sprintf("%10.2f", typ)  + sprintf("%10.2f", txp)   +  \
                      sprintf("%10.1f", vm_1) + sprintf("%10.1f", pc_1)  +  \
                      sprintf("%10.1f", rm_1) + sprintf("%10.1f", r34_1) +  \
                      sprintf("%10.1f", pe_1) + sprintf("%10.2f", lat_1) + sprintf("%10.2f", lon_1) + \
                      sprintf("%10.1f", vt_1) + sprintf("%10.1f", dir_1)  +  \
                     "  " + yyyy_1(0) + "" + mmddhhmm_s + "  " + name_1

        dir_pbl = dir_out + name_ib(ns) + "_" + year_ib(ns,ind_0)
        print(dir_pbl + " out" )
        flnm = dir_pbl + "/bdy_10min.txt"
        print("output file:" + flnm)
        mkdir_cmd = "csh -c 'if (! -d " + dir_pbl + ") then ;  mkdir " + dir_pbl + "; endif'"
        system(mkdir_cmd)
        asciiwrite(flnm, line_s)

        if(is_dat) then
          filo = dir_out + "landuse.dat"       ; output file
          print("output file: "+ filo)
          fbinwrite (filo, landuse)

          filo = dir_out + "topo.dat"       ; output file
          print("output file: "+ filo)
          fbinwrite (filo, topo)

          filo = dir_out + "lat_2d.dat"       ; output file
          print("output file: "+ filo)
          fbinwrite (filo, lat_2d)

          flnm =  dir_out + "topo.nc"
  	  print("output file: "+ flnm)
  	  system("/bin/rm -f " + flnm)
  	  ncdf = addfile( flnm, "c")

  	  fAtt             = True
  	  fAtt@script      = "James Done: preprocess_historical.ncl"
  	  fileattdef(ncdf, fAtt )

  	  ncdf->topo   = topo
  	  ncdf->landuse = landuse
  	  ncdf->lat_2d = lat_2d
  	  ncdf->lon_2d = lon_2d

          print("set is_dat=False, it is only used once")
          exit
        end if
        delete(line_s)
      end if
      print(lon_ib_1d + " lon/lat " + lat_ib_1d)
      dum(ct) = gsn_add_polyline(wks, map, lon_ib_1d(ind(lon_ib_1d)), lat_ib_1d(ind(lon_ib_1d)), rest)
      text_st(ct) = gsn_add_text(wks, map, name_ib(ns), lon_ib_1d(0), lat_ib_1d(0), txres)
      ct =  ct + 1
      print(ct + " ct")
      delete([/yyyy_1, name_1, mmddhhmm_s, mmddhhmm_1, non_missing_rm_1/])
      delete([/xo, xo_juli, xi, mmddhh_ib_1d, greg, lat_1, lon_1, vm_1, rm_1, r34_1, pc_1, pe_1, loc, txp, typ/])
      delete([/vt_ib_1d, dir_ib_1d, dir_1, vt_1, ind_r/])
      delete([/ vmax_ib_1d, rmax_ib_1d, r34_ib_1d , pmin_ib_1d, pe_ib_1d/])
      delete([/lat_ib_1d, lon_ib_1d,  juli_ib_1d /])
    end if
  end do	; end of ns-loop

  draw(map)
  frame(wks)

  print(ct + " TC")
end

