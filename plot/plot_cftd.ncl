;******************************************************************
; plot_cfodd.ncl
;
; Concepts illustrated:
;   - plot CFODD using numbers of occurence at each temp&dBZ pixel.
;
;*****************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
;
;************************************************
begin
 
  temp_min=195.0
  temp_max=305.0
  num_tempbin=20

  ze_min=-30.
  ze_max=20.
  num_zebin=25

 ;--temp bin setup--
  tempbnd=new((/num_tempbin+1/),"float")
  tempbnd_center=new((/num_tempbin/),"float")
  tempbnd=fspan(temp_min,temp_max,num_tempbin+1)
  do ic=0,num_tempbin-1
     tempbnd_center(ic)=0.5*(tempbnd(ic)+tempbnd(ic+1))
  end do
 ;--ze bin setup--
  zebnd=new((/num_zebin+1/),"float")
  zebnd_center=new((/num_zebin/),"float")
  zebnd=fspan(ze_min,ze_max,num_zebin+1)
  do iz=0,num_zebin-1
     zebnd_center(iz)=0.5*(zebnd(iz)+zebnd(iz+1))
  end do

  filpath="../cfodd_output/cftd_aero_catego_br/"
  naer="naer4"
  pdfcz1=asciiread(filpath+"cftd_naer1.txt",(/num_tempbin,num_zebin/),"float")
  pdfcz2=asciiread(filpath+"cftd_naer2.txt",(/num_tempbin,num_zebin/),"float")
  pdfcz3=asciiread(filpath+"cftd_naer3.txt",(/num_tempbin,num_zebin/),"float")
  pdfcz4=asciiread(filpath+"cftd_naer4.txt",(/num_tempbin,num_zebin/),"float")
  pdfcz5=asciiread(filpath+"cftd_naer5.txt",(/num_tempbin,num_zebin/),"float")

  ze_del=(ze_max-ze_min)/num_zebin ;ze bin width
  do itemp=0,num_tempbin-1
    ncnt=sum(pdfcz1(itemp,:))       ;number of smaples in each temp bin
   ;---convert numbers to PDF---
    if (ncnt.gt.0) then
    pdfcz1(itemp,:)=pdfcz1(itemp,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz2(itemp,:))
    if (ncnt.gt.0) then
    pdfcz2(itemp,:)=pdfcz2(itemp,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz3(itemp,:))
    if (ncnt.gt.0) then
    pdfcz3(itemp,:)=pdfcz3(itemp,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz4(itemp,:))
    if (ncnt.gt.0) then
    pdfcz4(itemp,:)=pdfcz4(itemp,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz5(itemp,:))
    if (ncnt.gt.0) then
    pdfcz5(itemp,:)=pdfcz5(itemp,:)/ncnt/ze_del*100.
    end if
    ;ncnt=sum(pdfcz6(itemp,:))
    ;if (ncnt.gt.0) then
    ;pdfcz6(itemp,:)=pdfcz6(itemp,:)/ncnt/ze_del*100.
    ;end if
  end do

  pdfcz1!0="temp"
  pdfcz1!1="ze"
  pdfcz1&temp=tempbnd_center
  pdfcz1&ze=zebnd_center
 
  copy_VarCoords(pdfcz1,pdfcz2)
  copy_VarCoords(pdfcz1,pdfcz3)
  copy_VarCoords(pdfcz1,pdfcz4)
  copy_VarCoords(pdfcz1,pdfcz5)

;---plot figure---

  wks  = gsn_open_wks("x11","cfodd_br68_modistau") 
  res  = True
  res@gsnDraw=False
  res@gsnFrame=False
  res@tiXAxisString   = "Reflectivity [dBZ]"
  res@tiYAxisString   = "Optical Depth"
  res@trXMinF = -30
  res@trXMaxF = 20
  res@trYMinF = 0
  res@trYMaxF = 60
  res@trYReverse = True

  res@cnLinesOn = False
  res@cnFillOn = True
  res@cnFillPalette = "NCV_jet"
  ;res@gsnSpreadColors = True
  res@gsnSpreadColorStart = 2
  res@gsnSpreadColorEnd = 241

  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = fspan(1,10,10)

  res@lbLabelBarOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False

  res@tmXBLabelFontHeightF = 0.03
  res@tmYLLabelFontHeightF = 0.03

  plot=new((/5/),"graphic")

  res@tiXAxisFontHeightF = 0.04
  res@tiYAxisFontHeightF = 0.04
  res@gsnLeftString = "(a)"
  res@gsnCenterString = "Naer range 1"
  res@tiXAxisString = ""
  res@tiYAxisString = "Temperature"
  res@gsnLeftStringFontHeightF = 0.04
  res@gsnCenterStringFontHeightF = 0.04
  res@tmXBMinorOn=False
  res@tmYLMinorOn=False
  res@tmXTOn=False
  res@tmYROn=False

  plot(0) = gsn_csm_contour (wks,pdfcz1, res)
  res@tiYAxisString = ""
  res@tiXAxisString   = "Reflectivity (dBZ)"
  res@gsnLeftString = "(b)"
  res@gsnCenterString = "Naer range 2"
  res@tiMainString = "MIROC5(BR68)"
  res@tiMainFontHeightF = 0.04
  ;res@tiMainOffsetYF = -0.01
  res@tiXAxisString   = "Reflectivity (dBZ)"
  plot(1) = gsn_csm_contour (wks,pdfcz2, res)
  res@gsnLeftString = "(c)"
  res@gsnCenterString = "Naer range 3"
  res@tiYAxisString = ""
  res@tiXAxisString   = ""
  res@tiMainString = ""
  plot(2) = gsn_csm_contour (wks,pdfcz3, res)
  res@gsnLeftString = "(d)"
  res@gsnCenterString = "Naer range 4"
  res@tiYAxisString = ""
  res@tiXAxisString   = ""
  res@tiMainString = ""
  plot(3) = gsn_csm_contour (wks,pdfcz4, res)
  res@gsnLeftString = "(e)"
  res@gsnCenterString = "Naer range 5"
  res@tiYAxisString = ""
  res@tiXAxisString   = ""
  res@tiMainString = ""
  plot(4) = gsn_csm_contour (wks,pdfcz5, res)

  respn=True
  respn@gsnPanelLabelBar = True
  respn@lbAutoManage          = False
  respn@lbBoxLinesOn     = False
  respn@pmLabelBarWidthF = 0.08
  respn@pmLabelBarHeightF = 0.23
  respn@lbOrientation = "Vertical"
  respn@lbLabelStride = 1
  respn@lbLabelFontHeightF = 0.013
  ;respn@pmLabelBarOrthogonalPosF = 0.0
  respn@pmLabelBarParallelPosF = 0.015

  gsn_panel(wks, plot(:), (/2,3/),respn)
end
