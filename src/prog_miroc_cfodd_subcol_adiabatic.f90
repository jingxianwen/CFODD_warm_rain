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
integer, parameter :: nmax_reff = 6, nmax_tauc = 15

!parameters for Re pdf.
integer, parameter :: num_rebnd = 200
real :: rebnd(num_rebnd+1)  ! bands for pdf calculation

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
real :: zetot(nlev,nlon,nlat) ! 
real :: frout(nlev,nlon,nlat) ! 

!output variables
real :: cntbnd(num_rebnd)  !count, for pdf calculation
real :: bndmid(num_rebnd)  !middle of bnd, for pdf calculation
integer :: cnt_cfodd(nmax_reff,nmax_tauc,nmax_dbz)

!local variables
real :: cldthk 
real :: cldhpos
real :: qliq(nlev)  ! mass mixing ration of liq cloud+rain. (kg/kg)
real :: qice(nlev)  ! mass mixing ration of ice cloud+rain. (kg/kg)
real :: tctop  ! cloud top temperature. (K) 
real :: tau0 
real :: re  ! effective radius
real :: nc  ! number of liquid cloud droplets

real :: reff_min, reff_max, reff_del !, reff_bin( nmax_reff ) 
real :: tauc_min, tauc_max, tauc_del !, tauc_bin( nmax_tauc )
real :: dbz_min, dbz_max, dbz_del    !, dbz_bin( nmax_dbz )

real :: long(nlon) !longitude
real :: long_bg(2), long_ed(2) ! longitude begin/end for data sample, defined by A-Train overpass time.

integer :: icotbin, irebin, izebin
real :: tau
integer :: cbase_idx, ctop_idx
logical :: cld_idx
logical :: sgl_idx

!settings
character(50) :: data_path
character(80) :: out_path
character(15) :: filenm 
character(10) :: rerng(nmax_reff)
real, parameter :: qliq_thr=1.0e-4 !(kg/kg)

character(3) :: idcol 
character(2) :: ctime
character(4) :: hhmm !hour/minute in the day (GMT)
character(4) :: hhmm_all(4) !hour/minute in the day (GMT at which the data are saved)

!others
integer :: bidx
integer :: it, iz,ilat,ilon,ilev, irg,isub, ix, iy, kk, i, j, k
integer :: idrec,idrec2


!----set data and output paths----
  data_path="../inputdata/test_br_cosp2/y2502/grd/"
  out_path="../cfodd_output/cfodd_br_cosp2/"
 
!----set Re, Tau, and dbz bins for CFODD----
  reff_min = 0.0
  reff_max = 30.0
  reff_del = ( reff_max-reff_min )/nmax_reff
  rerng=(/"0-5  ","5-10 ","10-15","15-20","20-25","25-30"/)
  
  tauc_min = 0.0
  tauc_max = 60.0
  tauc_del = ( tauc_max-tauc_min )/nmax_tauc
  
  dbz_min = -30.0
  dbz_max =  20.0
  dbz_del = ( dbz_max-dbz_min )/nmax_dbz
  
!----set Re bins for PDF calculation----
  rebnd(:)=0.0
  bndmid(:)=0.0
  do i=2,num_rebnd
   rebnd(i)=rebnd(i-1)+0.25
   bndmid(i-1)=(rebnd(i)+rebnd(i-1))*0.5
  end do
  
!----read land_sea mask data----
  open(10,file="../utils/landsea_mask.txt")
  do iy=1,nlat
    read(10,*) landsea_mask(:,iy)
  end do
  close(10)
!----set longitude----
  do ix=1,nlon
    long(ix)=1.40625*(ix-1)
  end do
  hhmm_all=(/"0000","0600","1200","1800"/)
  
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

 !----initialize counters----
  cnt_cfodd(:,:,:)=0
  cntbnd(:)=0

 !----start loop over time steps----
  do it=1,ntmax
    
    hhmm=hhmm_all(mod(it,4)+1)
  !..sample by local time 
  !..only 01:30 and 13:30 local time are used
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

    read(16,rec=it) taumodis(:,:)
    do iz=1,nlev
      idrec=(it-1)*nlev+iz
      read(11,rec=idrec) qc(iz,:,:)
      read(12,rec=idrec) qi(iz,:,:)
      read(13,rec=idrec) tem(iz,:,:)
      read(14,rec=idrec) recld(iz,:,:)
      read(15,rec=idrec) tauliq(iz,:,:)
      read(17,rec=idrec) gdzm(iz,:,:)
    end do
 !----compute tau profile according to adiabatic growth model----
    tauadb(:,:,:)=0.0
    tauadbsgl(:,:,:)=0.0
    do ilat=1,nlat
     do ilon=1,nlon
       ctop_idx=-1
       cbase_idx=-1
       if (taumodis(ilon,ilat) .gt.0) then
         ! where is cloud base
         do ilev=1,nlev
            if (qc(ilev,ilon,ilat) .gt. qliq_thr) then
               cbase_idx=ilev
               exit
            end if
         end do
         ! where is cloud top
         do ilev=nlev,1,-1
            if (qc(ilev,ilon,ilat) .gt. qliq_thr) then
               ctop_idx=ilev
               exit
            end if
         end do
         ! is single layer?
         if (any(qc(cbase_idx:ctop_idx,ilon,ilat) .lt. qliq_thr )) then
           continue
         else
           cldthk=gdzm(ctop_idx+1,ilon,ilat)-gdzm(cbase_idx,ilon,ilat)
           do ilev=ctop_idx,cbase_idx,-1 
             cldhpos=cldthk - (gdzm(ctop_idx+1,ilon,ilat)-gdzm(ilev,ilon,ilat)) 
             tauadb(ilev,ilon,ilat)=taumodis(ilon,ilat)*(1.0-(cldhpos/cldthk)**(5.0/3.0))
             tauadbsgl(ilev,ilon,ilat)=tauadb(ilev,ilon,ilat)-tauadb(ilev+1,ilon,ilat)
           end do 
         end if
       end if
     end do !ilon
    end do !ilat

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
       if ((long(i) .ge. long_bg(1) .and. long(i) .le. long_ed(1)) .or.  &
           (long(i) .ge. long_bg(2) .and. long(i) .le. long_ed(2))) then
         continue
       else
         cycle
       end if

       do j=1,nlat 
        !over ocean only
         if (landsea_mask(i,j) .eq. 0) then !0:ocean. 1:land.
           qliq(:)=qc(:,i,j) !*rho(:,i,j)
           qice(:)=qi(:,i,j) !*rho(:,i,j)
           tctop=-999.9
           do k=nlev,1,-1
            if ((qliq(k)+qice(k)).ge.qliq_thr) then 
              tctop=tem(k,i,j)  
              exit
            end if
           end do
          !warm clouds only
           if (tctop .gt. 273.15) then 
               cld_idx=.false.
               do k=1,nlev
                 if (qliq(k) .ge. qliq_thr .and. zetot(k,i,j) .ge. -30.) then
                   cbase_idx=k
                   exit
                 end if
               end do
               do k=nlev,1,-1
                if (qliq(k) .ge. qliq_thr .and. zetot(k,i,j) .ge. -30. ) then
                    ctop_idx=k
                    cld_idx=.true.
                    exit
                 end if
               end do

               if (cld_idx .eqv. .true.) then
               !single layer cloud only
                sgl_idx=.true.
                if (any(qliq(cbase_idx:ctop_idx) .lt. qliq_thr ) .or. &
                    any(zetot(cbase_idx:ctop_idx,i,j) .lt. -30.) .or. &
                    any(tauadbsgl(cbase_idx:ctop_idx,i,j) .le. 0)) then
                   sgl_idx=.false.
                end if

                if (sgl_idx .eqv. .true.)  then
                  !determin Reff bin for CFODD
                  irebin=int( (recld(ctop_idx,i,j)*1.0e6-reff_min)/reff_del ) + 1
                  !Re PDF:
                  if (recld(ctop_idx,i,j) .gt. 0.0) then
                     bidx=maxloc(rebnd,1,rebnd.lt.recld(ctop_idx,i,j)*1.0e6)
                     cntbnd(bidx)=cntbnd(bidx)+1.0
                  end if

                  tau0=0.0
                  do k=ctop_idx,cbase_idx,-1
                  !cot profile
                     tau=tauadbsgl(k,i,j) +tau0
                     tau0=tau
                  !determin cot bin for CFODD
                    icotbin=int( (tau-tauc_min)/tauc_del ) + 1

                  !determin dbze bin for CFODD
                    izebin=int( (zetot(k,i,j)-dbz_min)/dbz_del ) + 1
  
                  !count for CFODD
                    if (irebin .ge.1 .and. irebin .le. nmax_reff .and. &
                        icotbin .ge. 1 .and. icotbin .le. nmax_tauc .and. &
                        izebin .ge. 1 .and. izebin .le. nmax_dbz) then
                      cnt_cfodd(irebin,icotbin,izebin)=cnt_cfodd(irebin,icotbin,izebin)+1
                    end if
                  end do !k
                 end if !sgl_idx
               end if !cld_idx
           end if !tctop
         end if !landsea_mask
       end do !lat
     end do !lon
   end do !isub

  end do !it

!----write out counter for CFODD
  do i=1,nmax_reff
   open (50,file=trim(out_path)//'cfodd_'//trim(rerng(i))//".txt", &
        form='formatted')
    do j=1,nmax_tauc
      write (50,*) cnt_cfodd(i,j,:)
    end do
   close (50)
  end do
  
  !-----write out Re PDF----
  open(60,file=trim(out_path)//'pdf_cnt.txt',form='formatted',status="replace")
  open(61,file=trim(out_path)//'pdf_center.txt',form='formatted',status="replace")
  write(60,*) cntbnd(:)
  write(61,*) bndmid(:)
 
end program
