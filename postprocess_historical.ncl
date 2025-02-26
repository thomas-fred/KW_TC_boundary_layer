; plot and writen out wind swath from PBL model 
; 10-min input PBL model output data is interpolated 
; UPDATE: is_vt  = True         ; read vt data from bdy_10min.txt file 
; Ming Ge Feb 2020 and modified by James Done Apr 2020

load "substring.ncl"
load "rd_IBTrACS.ncl" 
load "vt_latlon_min.ncl"

; JAMES. The fix for the too-high (gradient-like) winds at time zero is to write over the first time:
; cp maria_WILLOUBY_001.d maria_WILLOUBY_000.d


begin

; select the wind profile model used
  pro_wind = "WILLOUBY"
;  pro_wind = "H10"

; select the storm
  sens_case = "MARIA_2017/"
;  sens_case = "IRENE_2011/"
;  sens_case = "IRMA_2017/"

; set dir
  dir_in = "./"

; read vt from bdy_10min.txt file if it is available
  is_vt = True          

  bdy_file = "bdy_10min.txt"

  vmin = 16
  rad  = 4.*atan(1.0)/180.    ; degress to radians

  ; read in domain file
  flnm = "geo_pr.nc"

  print("input lat-lon file: " + flnm)
  ff     = addfile(flnm, "r")
  topo   = ff->HGT_M(0,:,:)
  land   = ff->LANDMASK(0,:,:)	; 1=land, 0=water
  lu     = ff->LU_INDEX(0,:,:)

  lat_2d = ff->XLAT_M(0,:,:)
  lon_2d = ff->XLONG_M(0,:,:)
  lon_2d = where(lon_2d.le.0, lon_2d+360, lon_2d)

  ds     = ff@DX

  dims   = dimsizes(lat_2d)
  
  topo@description = ""
  topo@units       = ""

  ny = dims(0)      ; lp
  nx = dims(1)      ; lq

  print("nx=" + nx + " ny=" + ny)
  
;  nk = 18
; changed this because the fortran code is no longer outputting every vertical level of pbl model
  nk = 1

  swath_wind           = topo
  swath_wind@units     = "m/s"
  swath_slp       = topo
  swath_slp@units = "hPa"

  res                       = True
  res@gsnDraw               = False            
  res@gsnFrame              = False
  res@gsnMaximize           = True
  res@gsnAddCyclic          = False       ; regional data
  res@tfDoNDCOverlay        = True        ; do not transform data, plot the data natively
  res@cnLinesOn             = True
  res@cnFillOn              = True
  res@cnLevelSelectionMode  = "ExplicitLevels"
;  res@cnLevels              = ispan(vmin, 49, 3)
  res@cnLevels              = ispan(vmin, 67, 3)
  res@pmTickMarkDisplayMode = "Always" 
  res@mpFillOn              = False 
  res@pmLabelBarDisplayMode = "Always" 
  res@pmLabelBarHeightF     = 0.1   
  res@mpProjection          = "LambertConformal"
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
  res@mpLambertMeridianF    = ff@CEN_LON
  res@mpLambertParallel1F   = ff@TRUELAT1
  res@mpLambertParallel2F   = ff@TRUELAT2
  res@lbAutoManage          = False             ; we control label bar
  res@mpLimitMode           = "Corners" 
  res@mpLeftCornerLatF      = lat_2d(0,0) 
  res@mpLeftCornerLonF      = lon_2d(0,0)
  res@mpRightCornerLatF     = lat_2d(dims(0)-1, dims(1)-1)
  res@mpRightCornerLonF     = lon_2d(dims(0)-1, dims(1)-1)  
 
  sres                      = True
  sres@gsnDraw              = False             ; do not draw the plot
  sres@gsnFrame             = False
  sres@gsnAddCyclic         = False 
  sres@tfDoNDCOverlay       = True 
  sres@cnLinesOn            = True
  sres@cnLevelSelectionMode = "ExplicitLevels"  
  sres@cnLineLabelsOn       = False
  sres@cnLineColor          = "white"
  sres@cnLineThicknessF     = 1.0
  sres@cnInfoLabelOn        = False
  sres@cnLevels 	    = ispan(800, 2400, 800)

  rest                       = True
  rest@gsLineThicknessF      = 3.5
  rest@gsLineColor           = "Black"

  cmap = read_colormap_file("wind_17lev")
;  cmap = read_colormap_file("precip2_17lev")
 

  mkres               = True
  mkres@gsMarkerIndex = 16     ; Filled circle
  mkres@gsMarkerSizeF = 0.012

  ;all_dir = systemfunc("ls " + dir_in + " |grep _20")
  ;all_dir = systemfunc("ls " + dir_in + " |grep _19;  ls " + dir_in + " |grep _20")
  all_dir = new((/1/),"string")
  all_dir(0) = dir_in + sens_case
  print(all_dir)
;  n_dir = dimsizes(all_dir) - 1
  n_dir = dimsizes(all_dir) 
  print( "Number of dir is " + n_dir)
  delim  = " "	; for bdy_10min.txt
  do nd = 0, n_dir - 1
    name_d = all_dir(nd)
    dir_out = dir_in + sens_case

    ; TC name is the 1st part of directory name
    name_s = str_get_field(sens_case, 1, "_")
    year_s = str_get_field(sens_case, 2, "_")
    print(" " + dir_out)

    ; read in pbl model output data
    ; read in boundry data to determin the map boundray of plots
    flnm_bdy = dir_out + bdy_file
    print("read in: " + flnm_bdy)

    delim  = " "
    data_t = asciiread(flnm_bdy, -1, "string")
    data_s = data_t(1:)
    n_hr   = dimsizes(data_s)
    n_field = str_fields_count(data_s(0), delim)
    ; IBTrACS tracks has 13 field (with storm speed and direction), EB 11
    print(n_field + " " + data_s(0))

;    yyyy_s = str_get_cols(str_get_field(data_s, 11, delim), 0, 3)
;    mm_s   = str_get_cols(str_get_field(data_s, 11, delim), 4, 5)
;    dd_s   = str_get_cols(str_get_field(data_s, 11, delim), 6, 7)
;    hh_s   = str_get_cols(str_get_field(data_s, 11, delim), 8, 9)
    x_lon  = stringtofloat(str_get_field(data_s,  1, delim))
    y_lat  = stringtofloat(str_get_field(data_s,  2, delim))
    lat_trk = stringtofloat(str_get_field(data_s, 8, delim))
    lon_trk = stringtofloat(str_get_field(data_s, 9, delim))
    rmax    = stringtofloat(str_get_field(data_s, 5, delim))
    vmax    = stringtofloat(str_get_field(data_s, 3, delim)) 

    if(n_field.eq.13)
      vt_spd  = stringtofloat(str_get_field(data_s, 10, delim)) ; ie reading the reduced Vt from bdy_10min.txt
      vt_dir  = stringtofloat(str_get_field(data_s, 11, delim))
      ; u      = -wspd*sin(wdir*rad)  ; u component (zonal)
      ; v      = -wspd*cos(wdir*rad)  ; v component (meridional)
      u_vt = vt_spd*sin(vt_dir*rad)  
      v_vt = vt_spd*cos(vt_dir*rad) 
 
      yyyymmddhhmm_s = str_get_field(data_s, 12, delim)
      yyyy_s = str_get_cols(str_get_field(data_s, 12, delim), 0, 3)
      mm_s   = str_get_cols(str_get_field(data_s, 12, delim), 4, 5)
      dd_s   = str_get_cols(str_get_field(data_s, 12, delim), 6, 7)
      hh_s   = str_get_cols(str_get_field(data_s, 12, delim), 8, 9)
      mn_s   = str_get_cols(str_get_field(data_s, 12, delim), 10, 11)

   else
      ; vt_latlon(lat_int, lon_int, n)  n: n minute data
      uv = vt_latlon_min(lat_trk, lon_trk, 10)
      u_vt  = uv[0]
      v_vt  = uv[1]
      vt_spd = sqrt(u_vt^2 + v_vt^2) ; ie calculating the full vt from the x and y track locations of the vortex center.
      vt_dir = vt_spd

      yyyymmddhhmm_s = str_get_field(data_s, 10, delim)
      yyyy_s = str_get_cols(str_get_field(data_s, 10, delim), 0, 3)
      mm_s   = str_get_cols(str_get_field(data_s, 10, delim), 4, 5)
      dd_s   = str_get_cols(str_get_field(data_s, 10, delim), 6, 7)
      hh_s   = str_get_cols(str_get_field(data_s, 10, delim), 8, 9)
      mn_s   = str_get_cols(str_get_field(data_s, 10, delim), 10, 11)

    end if

    ;  Translation speed correction
    ;  vt_eb = 1.173*vt_eb^0.63
    ; previously Ming applied the vt reduction for both reading the reduced vt from bdy_10min file and calculating the vt from x and y locations
    if (n_field.eq.13)   
    print("We read in the reduced Vt from bdy file so no need to further reduce it")
    else
    v_vt = where(v_vt.gt.0, 1.173*v_vt^0.63, -1.173*abs(v_vt)^0.63)  
    u_vt = where(u_vt.gt.0, 1.173*u_vt^0.63, -1.173*abs(u_vt)^0.63)
    end if

    rr_3d = new((/n_hr,ny,nx/), float)
    ; calculate track Vmax by adding back the Vt effect
    ; previously Ming has the next line for both reading the reduced vt from bdy_10min file and calculating the vt from x and y locations
    ; vmax = vmax + 1.173*vt_spd^0.63
    ; James corrected it to

    if (n_field.eq.13)
      vmax = vmax + vt_spd   ; ie adding back the reduced vt that was taken off prior to pbl run
    else
      vmax = vmax + 1.173*vt_spd^0.63 ; ie adding back the reduced vt as calc by x and y locations.
    end if


;    print("vt(m/s)=" + sqrt(u_vt*u_vt +v_vt*v_vt) + " vmax= " + vmax +  " ddhh=" + dd_s + "" + hh_s )

    do ii = 0, ny - 1
    do jj = 0, nx - 1
      rr = ds*sqrt((ii-y_lat)*(ii-y_lat)+(jj-x_lon)*(jj-x_lon))
      rr_3d(:, ii, jj) = 2.0*rmax*rr/(rr*rr+rmax*rmax)
    end do
    end do
  
    allfile_d = systemfunc("ls " + dir_out + "/*_" + pro_wind + "_???.d")
    print(allfile_d)
;    n_time = dimsizes(allfile_d) - 1
   n_time = dimsizes(allfile_d) 

    print(n_time)
; for code testing
;    n_time = 10  

    uu   = new((/n_time, ny, nx/), float)
    uu!0 = "time"
    uu!1 = "Y"
    uu!2 = "X"

    vv  = uu
    slp = uu
    prx = uu
    pry = uu

    uu@_FillValue = -99.
    vv@_FillValue = -99.
    slp@_FillValue = -99.

    ;rr_3d = 0
    ; jt0, it0  txp, typ

    do nh = 0, n_time - 1
;    do nh = 0 , 5 - 1

;; If time = zero, set values to missing.
;      if (nh.eq.0) then
;       uu(nh,:,:)=uu@_FillValue
;       vv(nh,:,:)=vv@_FillValue
;       slp(nh,:,:)=slp@_FillValue
;      else

      ct = 0
      flnm = allfile_d(nh)
      print("process: " + flnm + nx + " " + ny)
  
      uu(nh, :,:) = fbinrecread(flnm, ct, (/ny, nx/), "float") + u_vt(nh)*rr_3d(nh,:,:)
      ct = ct + nk
      vv(nh, :,:) = fbinrecread(flnm, ct, (/ny, nx/), "float") + v_vt(nh)*rr_3d(nh,:,:)
      ct = ct + 3
      slp(nh, :,:) = fbinrecread(flnm, ct, (/ny, nx/), "float") 
;      ct = ct + 1
;      prx(nh, :,:) = fbinrecread(flnm, ct, (/ny, nx/), "float")
;      ct = ct + 1
;      pry(nh, :,:) = fbinrecread(flnm, ct, (/ny, nx/), "float")

      print("ntime:" + nh + " max_speed:" + max(sqrt(uu(nh, :,:)*uu(nh, :,:)+ vv(nh, :,:)*vv(nh, :,:))))
      print("ntime:" + nh + " min_slp: " + min(slp) + " max_slp: " + max(slp) )
;      end if
    end do

; create time variable
    yyyy = stringtointeger(yyyy_s)
    mm  = stringtointeger(mm_s)
    dd   = stringtointeger(dd_s)
    hh   = stringtointeger(hh_s)
    mn   = stringtointeger(mn_s)

; calc seconds since 1900-01-01 00:00:00
   sc  = new((/n_time/), integer)
   sc  = 0

  units  = "seconds since 1900-01-01 00:00:00"
  fptime   = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,units, 0)

  print(fptime)


    spd_tmp = sqrt(uu*uu+vv*vv)

    swath_wind  = 0.
; note that for some reason the swath_wind and swath_slp will stay as missing if we include the first time in the simulation.
    swath_slp   = 1005.
    do nh = 0, n_time - 1
      swath_wind     = where(swath_wind.lt.spd_tmp(nh,:,:), spd_tmp(nh,:,:), swath_wind)
      swath_slp      = where(swath_slp.ge.slp(nh,:,:), slp(nh,:,:), swath_slp)
    end do	; end of nh-loop n_time

;    swath_wind = where(swath.le.tofloat(vmin-3), swath@_FillValue,  swath)
    swath_wind@_FillValue = -99.
    swath_slp@_FillValue = -99.

    res@tiMainString = name_s

    flnm_plot = dir_out + str_lower(name_s) + "_" + pro_wind 
    wks = gsn_open_wks("pdf", flnm_plot)
    gsn_define_colormap(wks, "wind_17lev")
   
    plot    = gsn_csm_contour_map(wks, swath_wind, res)
    plot_ov = gsn_csm_contour    (wks,  topo, sres)
    overlay(plot, plot_ov)

  ; read IBTrACS lat, lon, Vmax for plotting
    ;track = rd_IBTrACS(name_s, toint(yyyy_s(0)), min(lat_trk), max(lat_trk), min(lon_trk), max(lon_trk))
    ;lat_t = track[0]
    ;lon_t = track[1]
    ;wind_t = track[2]
 
    lat_t = lat_trk(::12)
    lon_t = lon_trk(::12)
    wind_t = vmax(::12) 

    ; plot trajectory and wind as marker
    trajectory = gsn_add_polyline(wks, plot, lon_t, lat_t, rest)
    marker = new(dimsizes(lat_t), graphic)
    do n = 0, dimsizes(lat_t)-1
      v_tmp = toint(wind_t(n) - vmin)/3
      v_tmp = where(v_tmp.gt.17, 17, v_tmp)
;      v_tmp = where(v_tmp.gt.11, 11, v_tmp)
      v_tmp = where(v_tmp.lt.1,  1,  v_tmp)
      mkres@gsMarkerColor  =  cmap(v_tmp ,:)
      marker(n) = gsn_add_polymarker(wks, plot, lon_t(n), lat_t(n), mkres)
    end do

    draw(plot)
    frame(wks)
 
    ; write out
    line_s = new((/nx*ny + 1/), string) 
    line_s(0) = "wind(m/s),    lat,     lon"
    ct = 0
    
    do nyy = 0, ny-1
    do nxx = 0, nx-1
      if(.not.ismissing(swath_wind(nyy, nxx))) then
        ct =  ct + 1
        line_s(ct) = sprintf("%8.2f", swath_wind(nyy, nxx)) + "," + sprintf("%8.2f", lat_2d(nyy, nxx)) + "," +	\
	     	   sprintf("%8.2f", lon_2d(nyy, nxx)) 
      end if
    end do
    end do

    flnm = dir_out + str_lower(name_s) + "_" + pro_wind + ".txt"
    print("write to " + flnm  )
    asciiwrite(flnm, line_s(:ct))

    flnm = dir_out + str_lower(name_s) + "_" + pro_wind + ".nc"
    print("write to " + flnm)
    system("/bin/rm -f " + flnm)
    ncdf = addfile(flnm, "c")

    dim_names = (/ "time", "south_north",  "west_east" /)
    dim_sizes = (/ n_time ,  ny ,  nx  /)
    dimUnlim = (/ True , False,  False  /)
    filedimdef( ncdf, dim_names, dim_sizes, dimUnlim )

    fAtt               = True            ; assign file attributes
    fAtt@title         = "Wind and sea level pressure data from the NCAR Tropical Cyclone model"
    fAtt@creation_date = systemfunc ("date")
    fileattdef(ncdf, fAtt )            ; copy file attributes

    swath_wind@description = "storm lifetime max wind speed"
    swath_slp@description = "storm lifetime min sea level pressure"

    fptime@long_name = "model time" ;
    fptime@standard_name = "time" ;
    fptime@base_date = "1900-01-01 00:00:00 UTC" 
    fptime@_FillValue = -9999
    fptime!0 = "time"

    uu@units = "meters per second"
    uu@description = "one-minute sustained U-wind at 10 meters above the surface"
    uu@_FillValue = -99.
    uu!0 = "time"
    uu!1 = "south_north"
    uu!2 = "west_east"

    vv@units = "meters per second"
    vv@description = "one-minute sustained V-wind at 10 meters above the surface"
    vv@_FillValue = -99.
    vv!0 = "time"
    vv!1 = "south_north"
    vv!2 = "west_east"

    slp@units = "hPa"
    slp@description = "sea level pressure"
    slp@_FillValue = -99.
    slp!0 = "time"
    slp!1 = "south_north"
    slp!2 = "west_east"

    ncdf->swath_wind = swath_wind
;    ncdf->swath_slp  = swath_slp
    ncdf->lat_2d     = lat_2d
    ncdf->lon_2d     = lon_2d  
;    ncdf->U          = uu(:,:,:)    
;    ncdf->V          = vv(:,:,:)
;    ncdf->SLP        = slp(:,:,:)
;    ncdf->time       = fptime

    delete([/data_t, data_s, yyyymmddhhmm_s, yyyy_s, mm_s, dd_s, hh_s, mn_s, x_lon, y_lat, lat_trk, lon_trk, rmax, vmax/])
;    delete([/data_t, data_s, yyyymmddhhmm, x_lon, y_lat, lat_trk, lon_trk, rmax, vmax/])
;    delete([/uu,vv, spd_tmp, slp, prx, pry, u_vt, v_vt, vt_spd, rr, rr_3d/])
    delete([/uu,vv, spd_tmp, slp, u_vt, v_vt, vt_spd, rr, rr_3d/])
    delete([/allfile_d,line_s, marker, lat_t, lon_t, wind_t/])
    delete([/vt_dir/])
  end do
end
