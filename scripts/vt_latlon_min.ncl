; calculate transational speed U and V (m/s) of TC from n-minute lat-lon position
; Ming Ge 2019

undef ("vt_latlon_min")

function vt_latlon_min(lat_trk:numeric, lon_trk:numeric, n_min)

local lat_trk, lon_trk, n_data, u_vt, v_vt, ds_u, ds_v, n_min

begin

  n_data = dimsizes(lat_trk)
  u_vt = new(n_data, float)
  v_vt = new(n_data, float)
  n_sec_n_min  = 60*n_min
  n_sec_2n_min = 2*n_sec_n_min

  do nn = 1, n_data - 2
    ; 3 in meter
    ds_u = gc_latlon(lat_trk(nn),   lon_trk(nn+1), lat_trk(nn), lon_trk(nn-1), 2, 3)
    ds_v = gc_latlon(lat_trk(nn+1), lon_trk(nn),   lat_trk(nn-1), lon_trk(nn), 2, 3)
    if(lon_trk(nn+1).gt.lon_trk(nn-1)) then
      u_vt(nn) = ds_u/n_sec_2n_min
    else
      u_vt(nn) = -ds_u/n_sec_2n_min
    end if

    if(lat_trk(nn+1).gt.lat_trk(nn-1)) then
      v_vt(nn) = ds_v/n_sec_2n_min
    else
      v_vt(nn) = -ds_v/n_sec_2n_min
    end if
  end do

  ds_v = gc_latlon(lat_trk(1), lon_trk(0), lat_trk(0), lon_trk(0), 2, 3)
  ds_u = gc_latlon(lat_trk(0), lon_trk(1), lat_trk(0), lon_trk(0), 2, 3)
  if(lon_trk(1).gt.lon_trk(0)) then
    u_vt(0) = ds_u/n_sec_n_min
  else
    u_vt(0) = -ds_u/n_sec_n_min
  end if

  if(lat_trk(1).gt.lat_trk(0)) then
    v_vt(0) = ds_v/n_sec_n_min
  else
    v_vt(0) = -ds_v/n_sec_n_min
  end if

  ds_u = gc_latlon(lat_trk(n_data-1), lon_trk(n_data-1), lat_trk(n_data-1), lon_trk(n_data-2), 2, 3)
  ds_v = gc_latlon(lat_trk(n_data-1), lon_trk(n_data-1), lat_trk(n_data-2), lon_trk(n_data-1), 2, 3)
  if(lon_trk(n_data-1).gt.lon_trk(n_data-2)) then
    u_vt(n_data-1) =  ds_u/n_sec_n_min
  else
    u_vt(n_data-1) = -ds_u/n_sec_n_min
  end if

  if(lat_trk(n_data-1).gt.lat_trk(n_data-2)) then
    v_vt(n_data-1) = ds_v/n_sec_n_min
  else
    v_vt(n_data-1) = -ds_v/n_sec_n_min
  end if

   return([/u_vt, v_vt/])
end
