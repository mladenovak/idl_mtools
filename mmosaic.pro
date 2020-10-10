
pro mosaic, outfile, pathsP, pathsPB, $
	sizeX=sizeX, sizeY=sizeY, raCenter=raCenter, decCenter=decCenter, $
	noregrid=noregrid, linregrid=linregrid, $
	calcRMS=calcRMS, calcEffFreq=calcEffFreq, rmsMap=rmsMap,$
	addweights=addweights

	nnn=N_ELEMENTS(pathsP) ;number of images to mosaic

	if nnn eq 0 then begin
		print, 'No images to mosaic!'
		return
	endif

	if KEYWORD_SET(calcEffFreq) then doEffFreq=1 else doEffFreq=0
	if not KEYWORD_SET(pathsPB) then begin
		pathsPB=strarr(nnn)
		pathsPB[*]=''
	endif
	if not KEYWORD_SET(addweights) then begin
		addweights=fltarr(nnn,/nozero)
		addweights[*]=1.0
	endif

	if N_ELEMENTS(pathsPB) ne nnn or N_ELEMENTS(addweights) ne nnn then begin
		print, 'Mismatch of input fields dimensions. Exiting.'
		return
	endif

	print, 'Final output will be ' + outfile

	centralFreq=0.
	centralFreqW=0.

	initialize=0
	for idx=0, N_ELEMENTS(pathsP)-1 do begin

		print, 'Mosaicing ' + str(idx+1) + ' out of ' + str(N_ELEMENTS(pathsP))

		pathP=pathsP[idx]

		fitsexists = FILE_TEST(pathP)
		if not fitsexists then begin
			print,"File doesn't exist: "+pathP
			continue
		endif

		print, 'pathP: ', pathP
		;read pointing
		dataP = readfits(pathP, headP)
		extast, headP, astrP

		freq=SXPAR(headP, 'CRVAL3')*1e-9

		;initialize data from first pointing
		if(initialize eq 0) then begin
			pathPB='$'

			; create PB response instead of reading it
			;if(not dospws) then begin
			;	response=getPointingResponse(headP,psize=getPointingSize(dataP,headP), feather=3)
			;	;response=readfits(pathsPB)
			;endif

			;take header of first pointing for the mosaic
			headM=headP

			; change headers according to new parameters
			if KEYWORD_SET(sizeX) then begin
				SXADDPAR, headM, 'NAXIS1', sizeX, format='I'
				SXADDPAR, headM, 'CRPIX1', sizeX/2+1, format='E19.12'
			endif else sizeX=(size(dataP))[1]

			if  KEYWORD_SET(sizeY) then begin
				SXADDPAR, headM, 'NAXIS2', sizeY, format='I'
				SXADDPAR, headM, 'CRPIX2', sizeY/2+1, format='E19.12'
			endif else sizeY=(size(dataP))[2]

			if  KEYWORD_SET(raCenter) then $
				SXADDPAR, headM, 'CRVAL1', raCenter, format='E19.12' $
			else raCenter=SXPAR(headP, 'CRVAL1')

			if  KEYWORD_SET(decCenter) then $
				SXADDPAR, headM, 'CRVAL2', decCenter, format='E19.12' $
			else decCenter=SXPAR(headP, 'CRVAL2')

			;initialize arrays for mosaic, weights and eff freq
			dataM=fltarr(sizeX, sizeY)
			dataMw=dataM
			if doEffFreq then dataFreq=dataM

			SXADDPAR, headM, 'OBJECT', 'Mosaic'
			SXADDPAR, headM, 'HISTORY', timenow()+' IDL start'

			extast, headM, astrM

			initialize=1
		endif


		SXADDPAR, headM, 'HISTORY', 'Adding: '+ FILE_BASENAME(pathP, '.fits')

		; get center
		centerPpix=astrP.CRPIX-1 ;fits convention iz 1-based
		xyad, headP, centerPpix[0],centerPpix[1],x,y
		centerPradec=[x,y]

		;get pixel position in mosaic for current pointing
		adxy, headM, centerPradec[0], centerPradec[1], pxm, pym
		;print, 'center in new coord:',[pxm, pym]
		print, 'center in float coord:',pxm,pym


		;PB
		newPB=pathsPB[idx]
		if(newPB ne pathPB) then begin
			pathPB=newPB

			if pathPB eq '' then begin
				print, 'Will not use any PB scaling (no 2nd input).'
				response=dataP*0+1.
			endif else begin
				response=readfits(pathPB, headW)
			endelse

			; create a pb response, works only for 3 GHz
			;response=getPointingResponse(headP,psize=getPointingSize(dataP,headP), feather=3) $
		endif

		addw=addweights[idx]
		if(addw ne 1.0) then $
			print, 'Additional weight factor = ' + str(addw)

		rms=1.
		; Noise weights: rmsmap or rms/pb
		if KEYWORD_SET(rmsMap) then begin
			print, 'RMS from second input map'
			weight=addw*(1./response)^2
			
			
			rms=mean(response[where(finite(response) eq 1 and response gt 0)])
			;rms=min(response[where(finite(response) eq 1)])

			print, 'mean RMS = '+str(rms)
			
		endif else begin
			if pathPB ne '' then $
				print, 'Weights scaled with the PB response'
			;RMS
			if KEYWORD_SET(calcRMS) then begin
				; fits a histogram on pb uncorrected image that has uniform noise
				; a bit dubious for widebandpbcorr images
				; maybe limit calc to center parts
				rmsGauss, dataP*response, rms,/ORIGUNIT,binsize=0.1e-6,/NOMIRROR
					
				print, 'RMS = '+str(rms)
			endif


			weight=addw*(response/rms)^2

		endelse

		; No fancy regridding, all images should have same pixel size and rot

		if KEYWORD_SET(noregrid) then begin
			;round the center pixel to nearest one in mosaic
			print, 'Nearest neighbour gridding.'
			CenterMpix = ROUND([pxm, pym],/L64) ;long
		endif else begin
			print, 'Regridding data:'
			;SXADDPAR, headM, 'HISTORY', 'Regridding data'
			if KEYWORD_SET(linregrid) then begin
				regrid, dataP, pxm, pym
				regrid, weight, pxm, pym
			endif else begin
				regrid, dataP, pxm, pym,/cubic
				regrid, weight, pxm, pym,/cubic
			endelse
			;must truncate position
			CenterMpix = FIX([pxm, pym],type=3) ;long
		endelse

		print, 'center in new coord:',CenterMpix
		SXADDPAR, headM, 'HISTORY', 'Centered at: ' + str(CenterMpix[0]) + ', '+str(CenterMpix[1])


		dataP=dataP*weight

		if doEffFreq then freqP=freq*weight

		; central frequency after stacking
		centralFreq+=freq*addw/rms^2
		centralFreqW+=1.*addw/rms^2


		sizeX=(size(dataP))[1]
		sizeY=(size(dataP))[2]


		; consider edges of the final image
		left=centerMpix[0]-centerPpix[0]
		bottom=centerMpix[1]-centerPpix[1]
		right=left+sizeX-1
		top=bottom+sizeY-1

		cropRight=min([astrM.naxis[0]-1 - right,0])
		cropTop=min([astrM.naxis[1]-1 - top,0])
		cropLeft=-min([ left,0])
		cropBottom=-min([ bottom,0])

		;print,	cropLeft,cropRight	,cropBottom,cropTop
		notInMosaic = abs(cropLeft) ge sizeX or abs(cropRight) ge sizeX $
			or abs(cropBottom) ge sizeY or abs(cropTop) ge sizeY

		;inf='Pointing['+str(cropLeft)+':'+str(cropRight-1)+', '+ str(cropBottom)+':'+str(cropTop-1)+']'+$
		;	' -> Mosaic['+str(left+cropLeft)+':'+str(right+cropRight)+', '+ str(bottom+cropBottom)+':'+str(top+cropTop)+']'


		if(not notInMosaic)then begin

			print, 'Copying data to mosaic'

			partP=dataP[cropLeft:(cropRight-1), cropBottom:(cropTop-1)]
			w=where(finite(partP) eq 0)

			partP[w]=0
			if doEffFreq then begin
				partFreq=freqP[cropLeft:(cropRight-1), cropBottom:(cropTop-1)]
				partFreq[w]=0			; NAN values have weight 0
			endif

			partW=weight[cropLeft:(cropRight-1), cropBottom:(cropTop-1)]
			partW[w]=0 			; NAN values have weight 0

			; add data to mosaic and weight map
			dataM[(left+cropLeft):(right+cropRight),(bottom+cropBottom):(top+cropTop)]+=partP
			dataMw[(left+cropLeft):(right+cropRight),(bottom+cropBottom):(top+cropTop)]+=partW

			if doEffFreq then begin
				dataFreq[(left+cropLeft):(right+cropRight),(bottom+cropBottom):(top+cropTop)]+=partFreq
			endif

			print, 'OK!'
		endif else begin
			SXADDPAR, headM, 'HISTORY', 'Not inside image!'
			print, 'Not inside image!'
		endelse
	endfor

	if initialize eq 0 then begin
		print, 'Could not initialize the mosaic. Exiting.'
		return
	endif

	dataM/=dataMw
	if doEffFreq then  dataFreq/=dataMw
	;Add central frequency in the header, in GHz
	SXADDPAR, headM, 'CRVAL3', centralFreq/centralFreqW*1e9, format='E19.12'

	;add min max values to header
	;SXADDPAR, headM, 'DATAMAX', max(dataM[where(finite(dataM))]), BEFORE='HISTORY'
	;SXADDPAR, headM, 'DATAMIN', min(dataM[where(finite(dataM))]), BEFORE='HISTORY'


	SXADDPAR, headM, 'HISTORY', timenow()+' IDL end'
	SXADDPAR, headM, 'DATE', timenow()

	print, 'Writing mosaic fits file.'
	writefits, outfile+'.fits', dataM, headM
	writefits, outfile+'.err.fits', 1./sqrt(dataMw), headM
	if doEffFreq then  writefits, outfile+'.freq.fits', dataFreq, headM

	print, 'Fits file written in ',outfile+'.fits'
end

pro regrid, im, xpix, ypix, cubic=cubic ;linear

	;xpix, ypicx - pixel position of new image in the underlaying grid

	; decimal part of pixel coordinate
	;shift ammount
	tx=0
	ty=0

	tx=1.-float(xpix-fix(xpix, type=3))
	if(KEYWORD_SET(ypix)) then ty=1.-float(ypix-fix(ypix, type=3))

	sizedim=(size(im))[0]

	sizeX=(size(im))[1]
	if(sizedim lt 2) then sizeY=0 else sizeY=(size(im))[2]



	;print,'sizex',sizex
	;print,'sizey',sizey

	;craete buffer image
	imBuff = im*0.


	;cubic interpol
	if(KEYWORD_SET(cubic)) then begin
		print, 'Regrid with CUBIC interpolation.'
		;create buffer image
		imBuff = im*0.

		c_1=0.5*(-tx^3+2.*tx^2-tx)
		c0=0.5*(3.*tx^3-5.*tx^2+2.)
		c1=0.5*(-3.*tx^3+4.*tx^2+tx)
		c2=0.5*(tx^3-tx^2)

		;mirror pixels on the edge to create the needed 4x4 kernel
		;
		;shift all columns
		for x = 0L, sizeX-1 do begin
			if x eq 0 then imBuff[x,*]=c_1*im[x+1,*]+c0*im[x,*]+c1*im[x,*]+c2*im[x+1,*]
			if x eq 1 then imBuff[x,*]=c_1*im[x-1,*]+c0*im[x-1,*]+c1*im[x,*]+c2*im[x+1,*]
			if x ge 2 and x lt sizeX-1 then imBuff[x,*]=c_1*im[x-2,*]+c0*im[x-1,*]+c1*im[x,*]+c2*im[x+1,*]
			if x eq sizeX-1 then imBuff[x,*]=c_1*im[x-2,*]+c0*im[x-1,*]+c1*im[x,*]+c2*im[x,*]
		endfor

		;for one dimensional, just copy from buffer
		if(sizeY eq 0) then im[*]=imbuff[*] else im[*]=0.

		c_1=0.5*(-ty^3+2.*ty^2-ty)
		c0=0.5*(3.*ty^3-5.*ty^2+2.)
		c1=0.5*(-3.*ty^3+4.*ty^2+ty)
		c2=0.5*(ty^3-ty^2)

		;shift all rows
		for y = 0L, sizeY-1 do begin
			if y eq 0 then im[*,y]=c_1*imBuff[*,y+1]+c0*imBuff[*,y]+c1*imBuff[*,y]+c2*imBuff[*,y+1]
			if y eq 1 then im[*,y]=c_1*imBuff[*,y-1]+c0*imBuff[*,y-1]+c1*imBuff[*,y]+c2*imBuff[*,y+1]
			if y ge 2 and y lt sizeY-1 then im[*,y]=c_1*imBuff[*,y-2]+c0*imBuff[*,y-1]+c1*imBuff[*,y]+c2*imBuff[*,y+1]
			if y eq sizeY-1 then im[*,y]=c_1*imBuff[*,y-2]+c0*imBuff[*,y-1]+c1*imBuff[*,y]+c2*imBuff[*,y]
		endfor

		;linear interpol
	endif else begin
		print, 'Regrid with LINEAR interpolation.'
		;shift all columns
		for x = 0L, sizeX-1 do begin
			if x eq 0 then begin
				;mirror first column
				imBuff[x,*]=im[x,*]
			endif else begin
				imBuff[x,*]=im[x-1,*]*(1.-tx)+im[x,*]*tx
			endelse
		endfor
		;for one dimensional, just copy from buffer
		if(sizeY eq 0) then im[*]=imbuff[*] else im[*]=0.
		;shift all rows
		for y = 0L, sizeY-1 do begin
			if y eq 0 then begin
				im[*,y]=imBuff[*,y]
			endif else begin
				im[*,y]=imBuff[*,y-1]*(1.-ty)+imBuff[*,y]*ty
			endelse
		endfor
	endelse
end

pro rmsGauss, im, rms, ra, dec, radius, head=head, plotX=ooo, plotYhisto=h, plotYgauss=gaus, nomirror=nomirror, origunit=origunit, binsize=binsize
	;ra&dec in deg
	;radius in arcmin, if 0 use all image
	;if no head, ra, dec and radius in pixels
	;rms is returned in muJy
	
	rms=0
	
	;print,'headdefined',KEYWORD_SET(head)
	if(KEYWORD_SET(radius)) then begin
		s=radius
		px=ra
		py=dec
	endif else s=0
		
	if(KEYWORD_SET(head)) then begin
		extast, head, astr
		adxy,head,ra,dec,px,py
		s=radius/60./astr.cdelt[1]
		;print,'s',s
	endif
	
	
	
	if s gt 0 then	begin
	
		l=round(px-s)
		r=round(px+s)
		b=round(py-s)
		t=round(py+s)
		
		maxx=(size(im))[1]-1
		maxy=(size(im))[2]-1
		
		l=min([max([0, l]), maxx])
		r=min([max([0, r]), maxx])
		b=min([max([0, b]), maxy])
		t=min([max([0, t]), maxy])
		
		
		aQuad=im[l:r,b:t]
	endif else aQuad=im
	
	
	aQuad=aQuad(where(finite(aQuad))) ; makni NaN 
	
	if not KEYWORD_SET(origunit) then aQuad*=1e6 ; Jy pretvori u mikroJy
	
	if(KEYWORD_SET(nomirror)) then begin
		z=aQuad
	endif else begin
		;take only negative noise and mirror it
		neg=where(aQuad lt 0,num)
		zer=where(aQuad eq 0,num2)
		
		z=[aQuad(neg), -aQuad(neg), aQuad(zer)]
	endelse
	
	if N_ELEMENTS(z) gt 5 then begin
		
		
		if(not KEYWORD_SET(binsize)) then binsize=0.1
	
		;print,'binsize',binsize

	
		h=hist1d(z, binsize=binsize,obin=ooo) ; histogram, obin su lokacije binova
		
		gaus = GAUSSFIT(ooo, h, param, NTERMS=3)
		rms=param[2]
	endif
end

function rmsMap, im, boxsize, increment

	;copy image
	mapa=im
	mapa[*]=!VALUES.F_NAN

	print,'boxsize',boxsize
	print,'increment',increment
	;return

	maxx=(size(im))[1]-1
	maxy=(size(im))[2]-1

	;needValidPix=0.05*boxsize^2
	needValidPix=3

	for i = 0L+increment/2, maxx, increment do begin

		for j = 0L+increment/2, maxy, increment do begin

			rms=!VALUES.F_NAN

			;using histogram + gauss method
			;rmsGauss,im,rms,i,j,boxsize/2,/ORIGUNIT,binsize=1e-6,/NOMIRROR

			l=min([max([0, i-boxsize/2]), maxx])
			r=min([max([0, i+boxsize/2]), maxx])
			b=min([max([0, j-boxsize/2]), maxy])
			t=min([max([0, j+boxsize/2]), maxy])

			;print,'l,r,b,t',l,r,b,t

			boxpix=im[l:r,b:t]
			w=where(finite(boxpix) eq 1, num)
			boxpix=boxpix[w]

			if(num lt needValidPix) then continue;

			for n=0, 30 do begin
				;print,'n',n

				rms=stddev(boxpix)
				;print,'rms',rms

				if(finite(rms) eq 0) then break

				mean=mean(boxpix)
				;w=where(abs(boxpix) lt abs(mean+3*rms), num)
				;w=where(boxpix lt mean+3*rms, num)

				w=where(abs(boxpix - mean) lt 3*rms, num)

				;print,'num',num
				if(N_ELEMENTS(boxpix) eq num) then break

				boxpix=boxpix[w]
			endfor
			;print,'n,rms',n,rms

			;print,'rms,sigma',rms,sigma

			l=min([max([0, i-increment/2]), maxx])
			r=min([max([0, i+increment/2]), maxx])
			b=min([max([0, j-increment/2]), maxy])
			t=min([max([0, j+increment/2]), maxy])

			mapa[l:r,b:t]=rms

		endfor
		print, 'Progress: '+string(float(i+1)/(maxx)*100, format='(i4)')+'%'
	endfor



	;interpolate
	print,'interpolate x'
	for i = 0L+increment/2, maxx-increment, increment do begin
		for sub=1.0, increment-1 do begin
			factor=sub/increment
			mapa[i+sub,*]=mapa[i,*]+factor*(mapa[i+increment,*]-mapa[i,*])
		endfor
	endfor

	print,'interpolate y'
	for j = 0L+increment/2, maxy-increment, increment do begin
		for sub=1.0, increment-1 do begin
			factor=sub/increment
			mapa[*,j+sub]=mapa[*,j]+factor*(mapa[*,j+increment]-mapa[*,j])
		endfor
	endfor

	;put NANs where if the original map have them
	mapa[where( finite(im) eq 0)] = !VALUES.F_NAN

	return, mapa


end
