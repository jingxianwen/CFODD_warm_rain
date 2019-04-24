;******************************************************************
; plot_cfodd.ncl
;
; Concepts illustrated:
;   - plot CFODD using numbers of occurence at each COT&dBZ pixel.
;
;*****************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
;
;************************************************
begin
 
  cot_min=0.0
  cot_max=60.0
  num_cotbin=15

  ze_min=-30.
  ze_max=20.
  num_zebin=25

 ;--cot bin setup--
  cotbnd=new((/num_cotbin+1/),"float")
  cotbnd_center=new((/num_cotbin/),"float")
  cotbnd=fspan(cot_min,cot_max,num_cotbin+1)
  do ic=0,num_cotbin-1
     cotbnd_center(ic)=0.5*(cotbnd(ic)+cotbnd(ic+1))
  end do
 ;--ze bin setup--
  zebnd=new((/num_zebin+1/),"float")
  zebnd_center=new((/num_zebin/),"float")
  zebnd=fspan(ze_min,ze_max,num_zebin+1)
  do iz=0,num_zebin-1
     zebnd_center(iz)=0.5*(zebnd(iz)+zebnd(iz+1))
  end do

  filpath="../cfodd_output/cfodd_br_cosp2/"

  pdfcz1=asciiread(filpath+"cfodd_0-5.txt",(/num_cotbin,num_zebin/),"float")
  pdfcz2=asciiread(filpath+"cfodd_5-10.txt",(/num_cotbin,num_zebin/),"float")
  pdfcz3=asciiread(filpath+"cfodd_10-15.txt",(/num_cotbin,num_zebin/),"float")
  pdfcz4=asciiread(filpath+"cfodd_15-20.txt",(/num_cotbin,num_zebin/),"float")
  pdfcz5=asciiread(filpath+"cfodd_20-25.txt",(/num_cotbin,num_zebin/),"float")
  pdfcz6=asciiread(filpath+"cfodd_25-30.txt",(/num_cotbin,num_zebin/),"float")

  ze_del=(ze_max-ze_min)/num_zebin ;ze bin width
  do icot=0,num_cotbin-1
    ncnt=sum(pdfcz1(icot,:))       ;number of smaples in each COT bin
   ;---convert numbers to PDF---
    if (ncnt.gt.0) then
    pdfcz1(icot,:)=pdfcz1(icot,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz2(icot,:))
    if (ncnt.gt.0) then
    pdfcz2(icot,:)=pdfcz2(icot,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz3(icot,:))
    if (ncnt.gt.0) then
    pdfcz3(icot,:)=pdfcz3(icot,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz4(icot,:))
    if (ncnt.gt.0) then
    pdfcz4(icot,:)=pdfcz4(icot,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz5(icot,:))
    if (ncnt.gt.0) then
    pdfcz5(icot,:)=pdfcz5(icot,:)/ncnt/ze_del*100.
    end if
    ncnt=sum(pdfcz6(icot,:))
    if (ncnt.gt.0) then
    pdfcz16(icot,:)=pdfcz16(icot,:)/ncnt/ze_del*100.
    end if
  end do

  pdfcz1!0="cot"
  pdfcz1!1="ze"
  pdfcz1&cot=cotbnd_center
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

  plot=new((/3/),"graphic")

  res@tiXAxisFontHeightF = 0.04
  res@tiYAxisFontHeightF = 0.04
  res@gsnLeftString = "(a)"
  res@gsnCenterString = "Re=5-10 ~F33~m~F21~m"
  res@tiXAxisString = ""
  res@tiYAxisString = "Optical Depth"
  res@gsnLeftStringFontHeightF = 0.04
  res@gsnCenterStringFontHeightF = 0.04
  res@tmXBMinorOn=False
  res@tmYLMinorOn=False
  res@tmXTOn=False
  res@tmYROn=False

  plot(0) = gsn_csm_contour (wks,pdfcz2, res)
  res@tiYAxisString = ""
  res@tiXAxisString   = "Reflectivity (dBZ)"
  res@gsnLeftString = "(b)"
  res@gsnCenterString = "Re=10-15 ~F33~m~F21~m"
  res@tiMainString = "MIROC5(BR68)"
  res@tiMainFontHeightF = 0.04
  ;res@tiMainOffsetYF = -0.01
  res@tiXAxisString   = "Reflectivity (dBZ)"
  plot(1) = gsn_csm_contour (wks,pdfcz3, res)
  res@gsnLeftString = "(c)"
  res@gsnCenterString = "Re=15-20 ~F33~m~F21~m"
  res@tiYAxisString = ""
  res@tiXAxisString   = ""
  res@tiMainString = ""
  plot(2) = gsn_csm_contour (wks,pdfcz4, res)

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

  gsn_panel(wks, plot(:), (/1,3/),respn)
end
