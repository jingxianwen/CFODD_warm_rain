program prog_nicam_cfodd_jsim

!-------------------------------------
! This program does two things:
! 1. compute the PDF of cloud top Re
! 2. compute CFODD statistics 
!   (COT profiles from adiabatic growth
!    calculation)
!-------------------------------------

implicit none

!parameters
integer, parameter :: nlat=128
integer, parameter :: nlon=256
integer, parameter :: nlev=40
integer, parameter :: ntmax=1459
integer, parameter :: nsub=25

integer, parameter :: nmax_dbz = 25
!!integer, parameter :: nmax_reff = 6, nmax_tauc = 15
integer, parameter :: nmax_temp = 20

!parameters for Re pdf.
integer, parameter :: num_rebnd = 200
real :: rebnd(num_rebnd+1)  ! bands for pdf calculation

!parameters for num_aerosol range.
integer, parameter :: num_naerbnd = 5
real :: naerbnd(num_naerbnd+1)  ! bands for aerosol number ranges.
character(5) :: naerbnd_str(num_naerbnd)  ! bands for aerosol number ranges.
integer :: n

!input variables
integer :: landsea_mask(nlon,nlat)
real :: tem(nlev,nlon,nlat) ! air temperature (K)
real :: qc(nlev,nlon,nlat)  ! mass mixing ration of liq cloud. (kg/kg)
real :: qi(nlev,nlon,nlat)  ! mass mixing ration of ice cloud. (kg/kg)
real :: recld(nlev,nlon,nlat)  ! reff of cloud (m)
real :: gdzm(nlev,nlon,nlat)  ! height of model interface level (m)
real :: tauliq(nlev,nlon,nlat)  ! cld optical depth profile
real :: tauadb(nlev,nlon,nlat)  ! cld optical depth profile, downward cumulated
real :: tauadbsgl(nlev,nlon,nlat)  !cld optical depth profile, sigle layer
real :: taumodis(nlon,nlat)  ! modis cld optical depth
real :: numaes(nlon,nlat)   ! vertically accumulated aerosol number 
real :: zetot(nlev,nlon,nlat) ! 
real :: frout(nlev,nlon,nlat) ! 

!output variables
real :: cntbnd_re(num_rebnd)  !count, for pdf calculation
real :: cntbnd_naer(num_naerbnd)  !count, for pdf calculation
real :: bndmid(num_rebnd)  !middle of bnd, for pdf calculation
!!integer :: cnt_cfodd(num_naerbnd,nmax_reff,nmax_tauc,nmax_dbz)
integer :: cnt_cftd(num_naerbnd,nmax_temp,nmax_dbz)

!local variables
real :: cldthk 
real :: cldhpos
real :: qliq(nlev)  ! mass mixing ration of liq cloud+rain. (kg/kg)
real :: qice(nlev)  ! mass mixing ration of ice cloud+rain. (kg/kg)
real :: qcld(nlev)  ! mass mixing ration of total cloud+rain. (kg/kg)
real :: tctop  ! cloud top temperature. (K) 
real :: tau0 
real :: re  ! effective radius
real :: nc  ! number of liquid cloud droplets

!!real :: reff_min, reff_max, reff_del !, reff_bin( nmax_reff ) 
!!real :: tauc_min, tauc_max, tauc_del !, tauc_bin( nmax_tauc )
real :: temp_min, temp_max, temp_del !, temp_bin( nmax_tauc )
real :: dbz_min, dbz_max, dbz_del    !, dbz_bin( nmax_dbz )

real :: latt(nlat) !latitude
real :: long(nlon) !longitude
real :: long_bg(2), long_ed(2) ! longitude begin/end for data sample by local time.
real :: latt_bg(2), latt_ed(2) ! latitude begin/end for data sample by local time.
logical,parameter :: sample_region=.true.
logical,parameter :: sample_time=.true.
real :: latt_rng(2) ! sample latitude range, only if sample_region=.true.
real :: long_rng(2) ! sample longitude range, only if sample_region=.true.

!!integer :: icotbin, irebin, izebin, inaerbin
integer :: izebin, inaerbin, itempbin
real :: tau
integer :: cbase_idx, ctop_idx
logical :: cld_idx
logical :: sgl_idx

!settings
character(50) :: data_path
character(80) :: out_path
character(15) :: filenm 
!!character(10) :: rerng(nmax_reff)
real, parameter :: qcld_thr=1.0e-5 !(kg/kg)

character(3) :: idcol 
character(2) :: ctime
character(4) :: hhmm !hour/minute in the day (GMT)
character(4) :: hhmm_all(4) !hour/minute in the day (GMT at which the data are saved)

!others
integer :: bidx
integer :: it, iz,ilat,ilon,ilev, irg,isub, ix, iy, kk, i, j, k, iaer
integer :: idrec,idrec2


!----set data and output paths----
  data_path="../inputdata/test_br_cosp2/y2502/grd/"
  out_path="../cfodd_output/cftd_aero_catego_br/"
 
!----set Re, Tau, and dbz bins for CFODD----
  !!reff_min = 0.0
  !!reff_max = 30.0
  !!reff_del = ( reff_max-reff_min )/nmax_reff
  !!rerng=(/"0-5  ","5-10 ","10-15","15-20","20-25","25-30"/)
  
  !!tauc_min = 0.0
  !!tauc_max = 60.0
  !!tauc_del = ( tauc_max-tauc_min )/nmax_tauc
  
  temp_min = 185.0
  temp_max = 305.0
  temp_del = ( temp_max-temp_min )/nmax_temp

  dbz_min = -30.0
  dbz_max =  20.0
  dbz_del = ( dbz_max-dbz_min )/nmax_dbz
  
!----set Re bins for PDF calculation----
  !!rebnd(:)=0.0
  !!bndmid(:)=0.0
  !!do i=2,num_rebnd
  !! rebnd(i)=rebnd(i-1)+0.25
  !! bndmid(i-1)=(rebnd(i)+rebnd(i-1))*0.5
  !!end do
  
!----set num_aerosol bins for CFODD stratification----
  !naerbnd(:)=(/1.e+5,2.0e+11,5.0e+11,8.0e+11,2.0e+12,1.0e+15/)
  naerbnd(:)=(/1.e+5,5.0e+12,1.5e+13,3.0e+13,4.5e+13,1.0e+15/)
  naerbnd_str(:)=(/"naer1","naer2","naer3","naer4","naer5"/)
  !bndmid(:)=0.0
  !do i=2,num_rebnd
  ! rebnd(i)=rebnd(i-1)+0.25
  ! bndmid(i-1)=(rebnd(i)+rebnd(i-1))*0.5
  !end do
  
!----read land_sea mask data----
  open(10,file="../utils/landsea_mask.txt")
  do iy=1,nlat
    read(10,*) landsea_mask(:,iy)
  end do
  close(10)
!----read latitude----
  open(10,file="../utils/lat_t85.txt")
  read(10,*) latt(:)
  close(10)
!----set longitude----
  do ix=1,nlon
    long(ix)=1.40625*(ix-1)
  end do

  hhmm_all=(/"0000","0600","1200","1800"/)

!---set sample region if required---   
  if (sample_region) then 
   latt_rng(:)=(/0.,45. /) 
   long_rng(:)=(/60.,140. /)
  end if 

!----set input data----
  open (11,file=trim(data_path)//"qliqlcp.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (12,file=trim(data_path)//"qicelcp.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (13,file=trim(data_path)//"T.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (14,file=trim(data_path)//"reff1cp.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (15,file=trim(data_path)//"dtauscp.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (16,file=trim(data_path)//"tauwmodis.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (17,file=trim(data_path)//"gdzmcp.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)
  open (180,file=trim(data_path)//"numaes_sprt.grd",form = "unformatted", access="direct", &
                status="old",recl=nlat*nlon*4)

 !----initialize counters----
  !cnt_cfodd(:,:,:)=0
  cnt_cftd(:,:,:)=0
  !!cntbnd_re(:)=0
  cntbnd_naer(:)=0

 !----start loop over time steps----
  do it=1,ntmax
    if (sample_time) then 
      hhmm=hhmm_all(mod(it,4)+1)
    !..sample by local time
     !only 01:30 and 13:30 local time are used, approximately the A-Train overpass time.
      if (hhmm .eq. "0000") then
        long_bg(:)=(/15., 195./)
        long_ed(:)=(/30., 210./)
      else if (hhmm .eq. "0600") then
        long_bg(:)=(/285., 105./)
        long_ed(:)=(/300., 120./)
      else if (hhmm .eq. "1200") then
        long_bg(:)=(/195., 15./)
        long_ed(:)=(/210., 30./)
      else if (hhmm .eq. "1800") then
        long_bg(:)=(/105.,285./)
        long_ed(:)=(/120.,300./)
      end if
    end if ! sample_time

    read(16,rec=it) taumodis(:,:)
    read(180,rec=it) numaes(:,:)
    do iz=1,nlev
      idrec=(it-1)*nlev+iz
      read(11,rec=idrec) qc(iz,:,:)
      read(12,rec=idrec) qi(iz,:,:)
      read(13,rec=idrec) tem(iz,:,:)
      read(14,rec=idrec) recld(iz,:,:)
      read(15,rec=idrec) tauliq(iz,:,:)
      read(17,rec=idrec) gdzm(iz,:,:)
    end do
    tem(:,:,:)=tem(:,:,:)*0.5
 !----compute tau profile according to adiabatic growth model----
    !!tauadb(:,:,:)=0.0
    !!tauadbsgl(:,:,:)=0.0
    !!do ilat=1,nlat
    !! do ilon=1,nlon
    !!   ctop_idx=-1
    !!   cbase_idx=-1
    !!   if (taumodis(ilon,ilat) .gt.0) then
    !!     ! where is cloud base
    !!     do ilev=1,nlev
    !!        if (qc(ilev,ilon,ilat) .gt. qliq_thr) then
    !!           cbase_idx=ilev
    !!           exit
    !!        end if
    !!     end do
    !!     ! where is cloud top
    !!     do ilev=nlev,1,-1
    !!        if (qc(ilev,ilon,ilat) .gt. qliq_thr) then
    !!           ctop_idx=ilev
    !!           exit
    !!        end if
    !!     end do
    !!     ! is single layer?
    !!     if (any(qc(cbase_idx:ctop_idx,ilon,ilat) .lt. qliq_thr )) then
    !!       continue
    !!     else
    !!       cldthk=gdzm(ctop_idx+1,ilon,ilat)-gdzm(cbase_idx,ilon,ilat)
    !!       do ilev=ctop_idx,cbase_idx,-1 
    !!         cldhpos=cldthk - (gdzm(ctop_idx+1,ilon,ilat)-gdzm(ilev,ilon,ilat)) 
    !!         tauadb(ilev,ilon,ilat)=taumodis(ilon,ilat)*(1.0-(cldhpos/cldthk)**(5.0/3.0))
    !!         tauadbsgl(ilev,ilon,ilat)=tauadb(ilev,ilon,ilat)-tauadb(ilev+1,ilon,ilat)
    !!       end do 
    !!     end if
    !!   end if
    !! end do !ilon
    !!end do !ilat

  !----loop through sub-columns----
   do isub=1,nsub
     print*,"it=",it," isub=",isub

     write(idcol,"(i3)") isub

     filenm=trim(adjustl(idcol))//".grd"

     open (18,file=trim(data_path)//"dbze"//trim(filenm),form = "unformatted",& 
              status="old",access="direct", recl=nlat*nlon*4)
     open (19,file=trim(data_path)//"cfrac"//trim(filenm),form = "unformatted",& 
              status="old",access="direct", recl=nlat*nlon*4)

     do iz=1,nlev
       idrec2=(it-1)*nlev+iz
       read(18,rec=idrec2) zetot(iz,:,:)
       read(19,rec=idrec2) frout(iz,:,:)
     end do !iz

     do i=1,nlon 
      !only 01:30 and 13:30 local time are used
       if (sample_time) then
         if ((long(i) .ge. long_bg(1) .and. long(i) .le. long_ed(1)) .or.  &
             (long(i) .ge. long_bg(2) .and. long(i) .le. long_ed(2))) then
           continue
         else
           cycle
         end if
       end if ! sample time 

       do j=1,nlat 

         if (sample_region) then
           if ((long(i) .ge. long_rng(1) .and. long(i) .le. long_rng(2)) .or.  &
               (latt(j) .ge. latt_rng(1) .and. latt(j) .le. latt_rng(2))) then
             continue
           else
             cycle
           end if
         end if ! sample_region

        !over ocean only
         if (landsea_mask(i,j) .eq. 0) then !0:ocean. 1:land.
           !!qliq(:)=qc(:,i,j) !*rho(:,i,j)
           !!qice(:)=qi(:,i,j) !*rho(:,i,j)
           qcld(:)=qc(:,i,j)+qi(:,i,j)
           !!tctop=-999.9
           !!do k=nlev,1,-1
           !! if ((qliq(k)+qice(k)).ge.qliq_thr) then 
           !!   tctop=tem(k,i,j)  
           !!   exit
           !! end if
           !!end do
          !warm clouds only
           !!if (tctop .gt. 273.15) then 
               cld_idx=.false.
               do k=1,nlev
                 !!if (qliq(k) .ge. qliq_thr .and. zetot(k,i,j) .ge. -30.) then
                 if (qcld(k) .ge. qcld_thr .and. zetot(k,i,j) .ge. -30.) then
                   cbase_idx=k
                   exit
                 end if
               end do
               do k=nlev,1,-1
                 !!if (qliq(k) .ge. qliq_thr .and. zetot(k,i,j) .ge. -30.) then
                 if (qcld(k) .ge. qcld_thr .and. zetot(k,i,j) .ge. -30. ) then
                    ctop_idx=k
                    cld_idx=.true.
                    exit
                 end if
               end do

               if (cld_idx .eqv. .true.) then
               !single layer cloud only
                sgl_idx=.true.
                if (any(qcld(cbase_idx:ctop_idx) .lt. qcld_thr ) .or. &
                    any(zetot(cbase_idx:ctop_idx,i,j) .lt. -30.)) then
                   sgl_idx=.false.
                end if

                if (sgl_idx .eqv. .true.)  then
                  !determin Reff bin for CFODD
                  !!irebin=int( (recld(ctop_idx,i,j)*1.0e6-reff_min)/reff_del ) + 1
                  !determin naer bin for CFODD
                  inaerbin=-1
                  if ( numaes(i,j) .gt. naerbnd(1) .and. numaes(i,j) .le. naerbnd(2)) then
                     inaerbin=1
                     cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  else if ( numaes(i,j) .gt. naerbnd(2) .and. numaes(i,j) .le. naerbnd(3)) then
                     inaerbin=2
                     cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  else if ( numaes(i,j) .gt. naerbnd(3) .and. numaes(i,j) .le. naerbnd(4)) then
                     inaerbin=3
                     cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  else if ( numaes(i,j) .gt. naerbnd(4) .and. numaes(i,j) .le. naerbnd(5)) then 
                     inaerbin=4
                     cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  else if ( numaes(i,j) .gt. naerbnd(5) .and. numaes(i,j) .le. naerbnd(6)) then
                     inaerbin=5
                     cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  end if
                  !!cntbnd_naer(inaerbin)=cntbnd_naer(inaerbin)+1.0
                  !print*,numaes(i,j),inaerbin
                  !pause
                  !Re PDF:
                  !!if (recld(ctop_idx,i,j) .gt. 0.0) then
                  !!   bidx=maxloc(rebnd,1,rebnd.lt.recld(ctop_idx,i,j)*1.0e6)
                  !!   cntbnd_re(bidx)=cntbnd_re(bidx)+1.0
                  !!end if

                  !!tau0=0.0
                  do k=ctop_idx,cbase_idx,-1
                  !cot profile
                     !!tau=tauadbsgl(k,i,j) +tau0
                     !!tau0=tau
                  !determin cot bin for CFODD
                    !!icotbin=int( (tau-tauc_min)/tauc_del ) + 1
                  !determin temperature bin for CFODD
                    itempbin=int( (tem(k,i,j)-temp_min)/temp_del ) + 1
                    !print*,tem(k,i,j),itempbin
                  !determin dbze bin for CFODD
                    izebin=int( (zetot(k,i,j)-dbz_min)/dbz_del ) + 1
  
                  !count for CFODD
                    if (inaerbin .ge.1 .and. inaerbin .le. num_naerbnd .and. &
                        itempbin .ge. 1 .and. itempbin .le. nmax_temp .and. &
                        izebin .ge. 1 .and. izebin .le. nmax_dbz) then
                      !cnt_cfodd(irebin,icotbin,izebin)=cnt_cfodd(irebin,icotbin,izebin)+1
                      cnt_cftd(inaerbin,itempbin,izebin)=cnt_cftd(inaerbin,itempbin,izebin)+1
                    end if
                  end do !k
                 end if !sgl_idx
               end if !cld_idx
           !!end if !tctop
         end if !landsea_mask
       end do !lat
     end do !lon
   end do !isub

  end do !it

!----write out counter for CFODD
  do iaer=1,num_naerbnd
   !do i=1,nmax_reff
    open (50,file=trim(out_path)//'cftd_'//trim(naerbnd_str(iaer))//".txt", &
         form='formatted')
     do j=1,nmax_temp
       write (50,*) cnt_cftd(iaer,j,:)
     end do
    close (50)
   !end do
  end do
  
  !-----write out Re PDF----
  !!open(60,file=trim(out_path)//'pdf_cnt_re.txt',form='formatted',status="replace")
  !!open(61,file=trim(out_path)//'pdf_center_re.txt',form='formatted',status="replace")
  open(62,file=trim(out_path)//'pdf_cnt_naer.txt',form='formatted',status="replace")
  !!write(60,*) cntbnd_re(:)
  !!write(61,*) bndmid(:)
  write(62,*) cntbnd_naer(:)
 
end program
