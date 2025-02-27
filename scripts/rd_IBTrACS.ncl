; get IBTrACS track info. lat, lon, wind
; given storm name and year, and lat-lon range
; Ming Ge Dec. 2019

undef("rd_IBTrACS")

function rd_IBTrACS(name_s:string, year_i:numeric, lat_1:numeric, lat_2:numeric, lon_1:numeric, lon_2)

local flnm_ib, f_ib, time, lat_ib, lon_ib, name_ib, wind_ib, dims, n_storm, utc_date, year_ib,	\
      lat, lon, wind_trk, ind_t, lat_t, lon_t, wind_t, lat_1, lat_2, lon_1, lon_2

begin
  flnm_ib = "/glade/work/mingge/TC_DBS/IBTrACS.ALL.v04r00.nc"
  print("function rd_IBTrACS read: "  + flnm_ib)

  f_ib    = addfile(flnm_ib,   "r")
  time    = f_ib->time               ; (storm, time) days since 1858-11-17 00:00:00
  lat_ib  = f_ib->lat   ; (storm, time)
  lon_ib  = f_ib->lon
  name_ib = f_ib->name
  ;wind_ib = f_ib->usa_wind
  dims    = dimsizes(lat_ib)
  n_storm = dims(0)

  utc_date  = ut_calendar(time, 0)      ;  (storm, time, 6
  year_ib   = toint(utc_date(:, :, 0))

  flnm_ib = "/glade/work/mingge/TC_DBS/vmax_ib.nc"
  f_ib    = addfile(flnm_ib,   "r")
  wind_ib = f_ib->vmax_ib        ;(storm, date_time)
  print("---- read in: " + flnm_ib )

  do ns = 0, n_storm - 1
    if(year_ib(ns, 0).eq.year_i.and.chartostring(name_ib(ns,:)).eq.name_s.and.\
       any(abs(lat_ib(ns,:)-lat_1).le.0.5)) then
      lat    = lat_ib(ns,:)
      lon    = lon_ib(ns,:)
      wind_trk = wind_ib(ns,:)
      ind_t = ind(lat.ge.lat_1.and.lat.le.lat_2.and.      \
                  lon.ge.lon_1.and.lat.le.lon_2.and.      \
                  .not.ismissing(wind_trk))
      lat_t = lat(ind_t)
      lon_t = lon(ind_t)
      wind_t = wind_trk(ind_t)
      return([/lat_t, lon_t, wind_t/])
    end if
  end do
end

; for a quick test
;load "rd_IBTrACS.ncl"

;begin

;  name_s = "JEBI"
;  year_i = 2018
;  track = rd_IBTrACS(name_s, year_i, 25, 60, 130, 160)
;  lat_t = track[0]
;  lon_t = track[1]
;  wind_t = track[2]
;  print(lat_t + " " + lon_t + " " + wind_t)
;end
