
pro mtools
	; test
	racen=150.1191667d
	decen=2.2058333d
	r=0.5*sqrt(2)
	print,sphere_area_deg2(racen-r,racen+r,decen-r,decen+r)

end

; compute radio spectrum at specific frequency assuming a power law with a spectral index of alpha
function radiospec, flux0, freq0, freq1, alpha
	return, flux0*(freq1/freq0)^alpha
end

; writes columns of data into a file, adapted function
pro writecol_ml, file, v1, v2, v3, v4, v5, v6, v7, v8, v9, $
	v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20,v21,v22,v23,v24,$
	v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,$
	FMT=fmt, FILNUM=filnum, head=head

	; writecol -- Writes a 2 column ascii file

	if (N_params() LT 2) then begin
		print,'Syntax - ' + $
			'writecol, file, v1, v2, [v3-v19] FMT=, FILNUM= '
		return
	endif

	flgvn = N_params()-1
	if not keyword_set( FMT ) then    flgfmt    = 0 else begin
		flgfmt = 1
		fmt = fmt[0]
	endelse

	if not keyword_set(FILNUM) then begin
		filnum = 91
		close, filnum
		openw, filnum, file, WIDTH=500
		flg_fil = 91
	endif
	
	if  keyword_set(HEAD) then begin
		printf, filnum, '# ' + head
	endif
	
	for i=0L,n_elements(v1)-1 do begin
		case flgvn of
			34: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i],v30[i],v31[i],v32[i],v33[i],v34[i]
			33: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i],v30[i],v31[i],v32[i],v33[i]
			32: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i],v30[i],v31[i],v32[i]
			31: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i],v30[i],v31[i]
			30: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i],v30[i]
			29: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i],v29[i]
			28: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i],v28[i]
			27: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i],v27[i]
			26: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i],v26[i]			
			25: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i],$
				v25[i]
			24: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i], v23[i],v24[i]
			23: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i],v23[i]
			22: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i],v22[i]
			21: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i], v21[i]
			20: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i], v20[i]
			19: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i], v19[i]
			18: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], $
				v16[i],v17[i], v18[i]
			17: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i]
			16: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i]
			15: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i]
			14: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i]
			13: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i],v13[i]
			12: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i],v12[i]
			11: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i],v11[i]
			10: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i],v10[i]
			9: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
				v8[i],v9[i]
			8: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i]
			7: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i]
			6: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i]
			5: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i], v5[i]
			4: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i]
			3: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i]
			2: printf, filnum, FORMAT=fmt, v1[i], v2[i]
			1: printf, filnum, FORMAT=fmt, v1[i]
			else: stop
		endcase
	endfor
	if keyword_set(FLG_FIL) then close, filnum
	return
end

; compute area in square degres
function area_deg2, ramin,ramax,decmin,decmax
	return,(ramax-ramin)*( sin(decmax*!pi/180.) - sin(decmin*!pi/180.) )*180./!pi
end

; Gaussian function, not normalized
function gauss, x, gausparams	
	A=gausparams[0]
	mu=gausparams[1]
	sigma=gausparams[2]
	
	; FWHM = 2*sqrt(ln(4))*sigma
	return, A*exp(-0.5*((x-mu)/sigma)^2 )
end

; Gaussian function, normalized
function gaussnorm, x, mu, sig

	; FWHM = 2*sqrt(ln(4))*sigma
	return, 1./(sig*sqrt(2*!pi))*exp(-0.5*((x-mu)/sig)^2 )
end

; generate equaly spaced numbers between two extreme values (including them as well)
function genrange, first, last, num, bin=bin
	; optional: define bin instead of total numbers

	if KEYWORD_SET(bin) then begin
		num=fix((last-first)/bin)+1
		return, findgen(num)*bin+first
	endif else begin
		if num eq 1 then return, first
		;if num eq 2 then return, [first,last]
		return, findgen(num)/float(num-1)*(last-first)+first
	endelse

	return,-1
end

; generate histograms (h) using user defined bin edges (bins)
pro histomla, bins, data, h, plotx, ploty, bincen, logcen=logcen
	; output histogram line for plotting (each bin reaches 0)
	; bincen are centers of bins (for fitting)
	; use logcen if geometrical mean for bincenter is prefered to arithmetical
	  
	; data at bin limits are handled as  a <= xxx < b
	h=histogram(value_locate(bins, data), min=0, max=N_ELEMENTS(bins)-2)

	plotx=fltarr((N_ELEMENTS(bins)-1)*4)
	ploty=fltarr((N_ELEMENTS(bins)-1)*4)

	for i = 0L, N_ELEMENTS(bins)-2 do begin
		plotx[i*4+0]=bins[i]
		plotx[i*4+1]=bins[i]
		plotx[i*4+2]=bins[i+1]
		plotx[i*4+3]=bins[i+1]

		ploty[i*4+0]=0
		ploty[i*4+1]=h[i]
		ploty[i*4+2]=h[i]
		ploty[i*4+3]=0
	endfor
	
	if KEYWORD_SET(logcen) then $
		bincen=10^((alog10(bins[0:-2])+alog10(bins[1:-1]))*0.5) $
	else bincen=(bins[0:-2]+bins[1:-1])*0.5

end

; calculate specific percentile of the input array
function percentile, arr, perc, gauss=gauss, interq=interq

	if(not KEYWORD_SET(perc)) then perc=[0,0.5,1]
	if(KEYWORD_SET(gauss)) then perc=[0.16,0.5,0.84]
	if(KEYWORD_SET(interq)) then perc=[0.25,0.5,0.75]

	n=N_ELEMENTS(arr)
	arr_pos=findgen(n)/(n-1.)
	perc=(arr[sort(arr)])[VALUE_LOCATE(arr_pos, perc)]
	
	return, perc
end

; new line string
function newline
	if (!D.NAME eq 'WIN') then nl = string([13B, 10B]) else nl = string(10B)
	return, nl
end

; tab string
function tab
	return,STRING(9B)
end

; convert to string function
function str, input, f=f
	;float and double will have god output now
	if( input eq !NULL)then return, ''
	
	if(KEYWORD_SET(f))then return, strtrim(string(input,format='(F)'),2) $
	else return, strtrim(input,2)
end

; convert to string, add a sign in front
function strsign,num
	if num lt 0 or num eq 0 then return, string(num, FORMAT = '(F0.1)')
	if num gt 0 then return ,'+'+string(num, FORMAT = '(F0.1)')
end

; get timestamp
function timeNow
	time = Systime(UTC=Keyword_Set(utc))
	
	;weekday = Strmid(time, 0, 3)
	day = String(StrMid(time, 8, 2), Format='(I2.2)')
	month = Strmid(time, 4, 3)
	year = Strmid(time, 20, 4)
	;timestring = Strmid(time, 11, 8)
	hour = Strmid(time, 11, 2)
	minute = Strmid(time, 14, 2)
	second = Strmid(time, 17, 2)
	months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
	m = (Where(months EQ StrUpCase(month))) + 1
	
	time=TIMESTAMP(DAY=day, MONTH=m, YEAR=year, HOUR=hour, MINUTE=minute, SECOND=second)
	
	time=Strmid(time,0,19)
	;if KEYWORD_SET(casa) then time+='.000000'
	
	return, time
end

; start plotting
pro openplot, filename, ps, xs, ys
	!p.multi=0
	if(not KEYWORD_SET(ps))then ps=0
	
	if(not KEYWORD_SET(xs))then xs=16
	if(not KEYWORD_SET(ys))then ys=16

	if(ps eq 1)then begin
		; bad greek letters problem?
		;!p.font=0  ;vektorski font
		set_plot, 'ps'
		
		!p.font=1
		device, /tt_font
		
		device, file=filename,/isolatin,/helv,/color,$
			xsize=xs,ysize=ys, font_size = 11, encapsulated=1, /portrait
	endif else begin
		set_plot, 'x'
		!p.font=0  ;vektorski font
		device,true_color=24,decomposed=0
		device, retain=2             ; window
		
		if xs ge ys then Window, XSIZE=512, YSIZE=512/xs*ys
		if xs lt ys then Window, XSIZE=512*xs/ys, YSIZE=512
	endelse
end

; end plotting
pro closeplot, ps
	;if(not KEYWORD_SET(ps))then ps=-1
	if (ps eq 1) then device, /close
	;openplot
end

; color presets
pro colors,black,white,gray,red,magenta,purple,blue,cyan,green,darkgreen,pink
	loadct,12, /silent
	
	black=0
	white=240
	gray=224
	red=192
	magenta=132
	purple=112
	blue=96
	cyan=80
	green=32
	darkgreen=16
	pink=208
end

; read a file into array
function readTxt, path
	OPENR, 1, path
	; Read one line at a time, saving the result into array
	array = ''
	line = ''
	WHILE NOT EOF(1) DO BEGIN & $
		READF, 1, line & $
		array = [array, line] & $
	ENDWHILE
	; Close the file and free the file unit
	close,1

	return, array
end

; write out region file for DS9, also a mask file for CASA
pro ds9reg,path, ras,decs, ids=ids, colors=colors, sizes=sizes, casa=casa
	num=N_ELEMENTS(ras)

	if(not KEYWORD_SET(colors)) then begin
		colors=strarr(num)
		colors[*]='green'
	endif else if (size(colors))[0] eq 0 then begin
		val=colors
		colors=strarr(num)
		colors[*]=val
	endif
	
	if(not KEYWORD_SET(ids)) then begin
		ids=strarr(num)
		ids[*]=''
	endif else if (size(ids))[0] eq 0 then begin
		val=ids
		ids=strarr(num)
		ids[*]=val
	endif
	
	if(not KEYWORD_SET(sizes)) then begin
		sizes=strarr(num)
		sizes[*]=5
	endif else if (size(sizes))[0] eq 0 then begin
		val=sizes
		sizes=fltarr(num)
		sizes[*]=val
	endif

	filename1=path
	openw,1,filename1
	printf,1,'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
	printf,1,'fk5'

	if(KEYWORD_SET(casa)) then begin
		filename2=filename1+'.casamask'
		openw,2,filename2
		printf,2,'#CRTFv0' ;CASA magic word
	endif

	for i = 0L, N_ELEMENTS(ras)-1 do begin

		printf,1, 'circle('+string(ras[i],format='(F11.6)')+','+string(decs[i],format='(F11.6)')+','+string(sizes[i],format='(F7.3)')+'") # text={'+str(ids[i])+'} color='+colors[i],format='(a)'
		
		if(KEYWORD_SET(casa)) then begin
			;printf,2,'circle[['+str(ras[i],/f)+'deg,'+str(decs[i],/f)+'deg], '+str(sizes[i])+'arcsec]'
			printf,2,'circle[['+string(ras[i],format='(F11.6)')+'deg,'+string(decs[i],format='(F11.6)')+'deg], '+string(sizes[i],format='(F7.3)')+'arcsec]'
		endif
	endfor

	close,1
	print, 'Region file written: ', filename1
	
	if(KEYWORD_SET(casa)) then begin
		close,2
		print, 'Clean file written: ', filename2
	endif
end

; write masks for outlier fields
pro ds9regFlank,path, ras,decs, ids=ids, colors=colors, sizes=sizes, casa=casa
	; size of the outlier image is currently fixed to 128x128pix regardless of the size input

	num=N_ELEMENTS(ras)

	if(not KEYWORD_SET(colors)) then begin
		colors=strarr(num)
		colors[*]='green'
	endif else if (size(colors))[0] eq 0 then begin
		val=colors
		colors=strarr(num)
		colors[*]=val
	endif

	if(not KEYWORD_SET(ids)) then begin
		ids=strarr(num)
		ids[*]=''
		for i = 0L, num-1 do begin
			ids[i]='out'+str(i)
		endfor
	endif else if (size(ids))[0] eq 0 then begin
		val=ids
		ids=strarr(num)
		ids[*]=val
	endif

	if(not KEYWORD_SET(sizes)) then begin
		sizes=strarr(num)
		sizes[*]=5
	endif else if (size(sizes))[0] eq 0 then begin
		val=sizes
		sizes=fltarr(num)
		sizes[*]=val
	endif

	filename1=path
	openw,1,filename1
	printf,1,'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
	printf,1,'fk5'

	if(KEYWORD_SET(casa)) then begin
		filename2=filename1+'.outlier'
		openw,2,filename2
		printf,2,'# Outlier info' ;CASA magic word
	endif

	for i = 0L, N_ELEMENTS(ras)-1 do begin
		printf,1, 'box('+string(ras[i],format='(F13.6)')+','+string(decs[i],format='(F13.6)')+','+string(sizes[i],format='(F7.3)')+'"'+','+string(sizes[i],format='(F7.3)')+'"'+',0) # text={'+str(ids[i])+'} color='+colors[i],format='(a)'

		if(KEYWORD_SET(casa)) then begin
			;printf,2,'circle[['+str(ras[i],/f)+'deg,'+str(decs[i],/f)+'deg], '+str(sizes[i])+'arcsec]'
			;printf,2,'circle[['+string(ras[i],format='(F13.6)')+'deg,'+string(decs[i],format='(F13.6)')+'deg], '+string(sizes[i],format='(F7.3)')+'arcsec]'

			printf,2,'imagename="'+str(ids[i])+'"'
			printf,2,'imsize=[128,128]'
			printf,2,'phasecenter="J2000 '+string(ras[i],format='(F13.6)') + 'deg ' +string(decs[i],format='(F13.6)')+'deg"'
		endif
	endfor

	close,1
	print, 'Region file written: ', filename1

	if(KEYWORD_SET(casa)) then begin
		close,2
		print, 'Outlier file written: ', filename2
	endif

end

; perform a 2D Gaussian imfit in CASA and parse the output
pro casafit, dir, file, ra, dec, s, params, paramnames, print=printres
	;ra, dec in degs, s in arcsec
	; must set this
	casa="path-to-casa-executable"
	
	params=[0]
	paramnames=['']
	
	cd,dir
	
	scriptfile='casa_imfit.py'
	logfile='casa_imfit.log'
	
	pytxt="imfit('"+file+"', region='centerbox[["+str(ra)+"deg, "+str(dec)+"deg], ["+str(s*2)+"arcsec"+","+str(s*2)+"arcsec"+"]]', logfile='"+logfile+"')"
	
	if(KEYWORD_SET(printres)) then print,pytxt
	
	openw,1,scriptfile
	printf,1,pytxt,format='(a)'
	close,1
	
	;remove existing content
	openw,1,logfile
	close,1
	
	;script=casa + " --nologger --nologfile -c " + scriptfile
	script=casa + " --nologger --nologfile -c " + scriptfile
	if(KEYWORD_SET(printres)) then print,script
	spawn, script,result,/STDERR
	

	result=readtxt(logfile)
	
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)
	r=""
	for i = 0L, N_ELEMENTS(result)-1 do begin
		r+=result[i]+newline
	endfor
	result=r
	
	if(KEYWORD_SET(printres)) then print,result
	
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*"
	
	rx=STREGEX(result, "--- Integrated:"+num+pm+num+" ([a-zA-Z]+)"+newline,/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
	
		fac=1.
		unit=rx[3]
		
		if unit eq 'uJy' then fac=1
		if unit eq 'mJy' then fac=1e3
		if unit eq 'Jy' then fac=1e6
		
		params=[params,double(rx[1]*fac)]
		params=[params,double(rx[2])*fac]
		paramnames=[paramnames, 'flux']
		paramnames=[paramnames, 'fluxerr']
	endif
	
	rx=STREGEX(result, "--- Peak:"+num+pm+num+" ([a-zA-Z]+)",/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
		fac=1.
		unit=rx[3]
		
		if unit eq 'uJy' then fac=1
		if unit eq 'mJy' then fac=1e3
		if unit eq 'Jy' then fac=1e6
		
		params=[params,double(rx[1]*fac)]
		params=[params,double(rx[2])*fac]
		paramnames=[paramnames, 'peak']
		paramnames=[paramnames, 'peakerr']
	endif
	
	rx=STREGEX(result, "Position ---"+newline+spc+"--- ra:"+num+pm+num + " s \("+num+" arcsec\)" + newline+$
		"[ ]*--- dec:"+num+pm+num+" arcsec",/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
		ra=rx[1]
		dec=rx[4]
		
		decpart=strsplit(dec,'.',/EXTRACT)
		if(N_ELEMENTS(decpart) eq 4) then dec=decpart[0]+':'+decpart[1]+':'+decpart[2]+'.'+decpart[3]
		if(N_ELEMENTS(decpart) eq 3) then  dec=decpart[0]+':'+decpart[1]+':'+decpart[2]
		
		stringad, (ra + ' ' + dec), ra,dec
		params=[params,double(ra)]
		params=[params,double(dec)]
		paramnames=[paramnames, 'ra']
		paramnames=[paramnames, 'dec']
	endif

	rx=STREGEX(result, "--- major axis FWHM:" + num + pm + num + " ([a-zA-Z]+)"+newline,/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
		fac=1.
		unit=rx[3]
		
		if unit eq 'marcsec' then fac=1e-3
		if unit eq 'arcsec' then fac=1.
		
		
		params=[params,double(rx[1]*fac)]
		params=[params,double(rx[2])*fac]
		paramnames=[paramnames, 'bmaj']
		paramnames=[paramnames, 'bmajerr']
	endif
	
	rx=STREGEX(result, "--- minor axis FWHM:" + num + pm + num + " ([a-zA-Z]+)"+newline,/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
		fac=1.
		unit=rx[3]
		
		if unit eq 'marcsec' then fac=1e-3
		if unit eq 'arcsec' then fac=1.

		params=[params,double(rx[1]*fac)]
		params=[params,double(rx[2])*fac]
		paramnames=[paramnames, 'bmin']
		paramnames=[paramnames, 'bminerr']
	endif
	
	rx=STREGEX(result, "--- position angle:" + num + pm + num + " deg"+newline,/SUBEXPR,/EXTRACT)
	if(rx[0] ne '') then begin
		params=[params,double(rx[1])]
		params=[params,double(rx[2])]
		paramnames=[paramnames, 'bpa']
		paramnames=[paramnames, 'bpaerr']
	endif
end

; automate a function from miriad: import
function mirImport, dir, file, print=pr, check=check
	cd, dir
	
	if(KEYWORD_SET(check))then begin
		ext=strmid(file,strlen(file)-5,5)
		
		if(ext eq '.fits' or ext eq '.FITS') then begin
			;check if file is already imported
			if(FILE_TEST(file+'.mir',/DIRECTORY) eq 0) then begin
				print,'Importing fits to Miriad.'
				;mirImport, dir, file
			endif else begin
				if(KEYWORD_SET(pr)) then print, 'Miriad image exists.'
				return, file+'.mir'
			endelse
		endif
		ext=strmid(file,strlen(file)-4,4)
		if(ext eq '.mir') then begin
			return, file
		endif
	endif
	
	print,dir + '/' + file
	
	script="fits in="+file+" out="+file+'.mir'+ " op=xyin"
	spawn, script, result,/STDERR
	
	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result
	
	return, file+'.mir'
end

; automate a function from miriad: adxy, coordinate conversion
pro mirAdxy, dir, file, ra, dec, px,py, value, nearpx, nearpy, print=pr,import=import

	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)
	
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)
	
	script="impos in="+mirfile+" type=absdeg coord="+str(ra,/f)+','+str(dec,/f)+""
	spawn, script, result,/STDERR
	
	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r
	
	if(KEYWORD_SET(pr)) then print, script
	if(KEYWORD_SET(pr)) then print, result
	
	rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---SIN =" + num + newline $
		+"Axis 2: DEC--SIN =" + num $
		,/SUBEXPR,/EXTRACT)
		
	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---TAN =" + num + newline $
		+"Axis 2: DEC--TAN =" + num $
		,/SUBEXPR,/EXTRACT)
		
	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---NCP =" + num + newline $
		+"Axis 2: DEC--NCP =" + num $
		,/SUBEXPR,/EXTRACT)

	px=-99
	py=-99
	if(rx[1] ne '')then px=double(rx[1])
	if(rx[2] ne '')then py=double(rx[2])
	
	nearpx=-99
	nearpy=-99
	value=-99
	rx=STREGEX(result, 'Nearest pixel = '+num+','+num+'\.'+spc+'Value ='+num+ '([ a-zA-Z]*)',/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then nearpx=long(rx[1])
	if(rx[2] ne '')then nearpy=long(rx[2])
	if(rx[3] ne '')then value=double(rx[3])
end

; automate a function from miriad: adxy, coordinate conversion
pro mirAdxyPath, path, ra, dec, px,py, value, nearpx, nearpy, print=pr,import=import
	file=FILE_BASENAME(path)
	dir=FILE_DIRNAME(path)

	cd, dir
	mirfile=file

	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)

	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	script="impos in="+mirfile+" type=absdeg coord="+str(ra,/f)+','+str(dec,/f)+""
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r

	if(KEYWORD_SET(pr)) then print, script
	if(KEYWORD_SET(pr)) then print, result

	rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---SIN =" + num + newline $
		+"Axis 2: DEC--SIN =" + num $
		,/SUBEXPR,/EXTRACT)

	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---TAN =" + num + newline $
		+"Axis 2: DEC--TAN =" + num $
		,/SUBEXPR,/EXTRACT)

	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result, 'Absolute pixels'+newline $
		+"Axis 1: RA---NCP =" + num + newline $
		+"Axis 2: DEC--NCP =" + num $
		,/SUBEXPR,/EXTRACT)

	px=-99
	py=-99
	if(rx[1] ne '')then px=double(rx[1])
	if(rx[2] ne '')then py=double(rx[2])

	nearpx=-99
	nearpy=-99
	value=-99
	rx=STREGEX(result, 'Nearest pixel = '+num+','+num+'\.'+spc+'Value ='+num+ '([ a-zA-Z]*)',/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then nearpx=long(rx[1])
	if(rx[2] ne '')then nearpy=long(rx[2])
	if(rx[3] ne '')then value=double(rx[3])
end

; automate a function from miriad: xyad, coordinate conversion
pro mirXyad, dir, file, px,py, ra,dec, value, nearpx, nearpy, print=pr,import=import

	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)
		
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)
	
	script="impos in="+mirfile+" type=abspix coord="+str(px,/f)+','+str(py,/f)+""
	
	spawn, script, result,/STDERR
	
	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r
	
	if(KEYWORD_SET(pr)) then print, script
	if(KEYWORD_SET(pr)) then print, result
	
	rx=STREGEX(result, 'World coordinates'+newline $
		+"Axis 1: RA---SIN =" + num + newline $
		+"Axis 2: DEC--SIN =" + num $
		,/SUBEXPR,/EXTRACT)

	ra=-99
	dec=-99
	if(rx[1] ne '' and rx[2] ne '') then stringad, rx[1]+' '+rx[2],ra,dec

	nearpx=-99
	nearpy=-99
	value=-99
	rx=STREGEX(result, 'Nearest pixel = '+num+','+num+'\.'+spc+'Value ='+num+ ' ([a-zA-Z/]+)',/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then nearpx=long(rx[1])
	if(rx[2] ne '')then nearpy=long(rx[2])
	if(rx[3] ne '')then value=double(rx[3])
end

; automate a function from miriad: statistics, compute rms
pro mirRMS, dir, file, ra, dec, s, rms, freq, print=pr,import=import ;s in arcmin

	cd, dir

	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)
	
	
	mirAdxy,dir,mirfile,ra+s/60.,dec-s/60.,px0,py0
	mirAdxy,dir,mirfile,ra-s/60.,dec+s/60.,px1,py1
	;script="imfit in="+mirfile+" object=gaussian region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	script="imlist in=" + mirfile + " options=statistics region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"

	spawn, script, result,/STDERR
	;print,result
	
	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result
	
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space	
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r

	rx=STREGEX(result, "plane"+spc+"Frequency"+spc+"Total Flux"+spc+"Maximum"+spc+"Minimum"+spc+"Average"+spc+"rms"+newline $
		+num+num+num+num+num+num+num,/SUBEXPR,/EXTRACT)

		if(rx[2] ne '')then freq=double(rx[2]) else freq=-99 ;GHz
		if(rx[7] ne '')then rms=double(rx[7]) else rms=-99 
	if(rms eq -1) then begin
		rx=STREGEX(result, spc+"Total Flux"+spc+"Maximum"+spc+"Minimum"+spc+"Average"+spc+"rms"+newline $
			+num+num+num+num+num,/SUBEXPR,/EXTRACT)
		if(rx[5] ne '')then rms=double(rx[5]) else rms=-99 
	endif
	
	;print,'rms',rms
	;print,'freq',freq
	
	;print,rx
	
	;Axis:     3     FREQ
	;plane   Frequency   Total Flux   Maximum     Minimum     Average       rms
	;1   2.05205      5.8008E-03  1.4230E-03 -1.1656E-04  3.6926E-07  2.6937E-05
	;if(KEYWORD_SET(pr)) then print,script
	;if(KEYWORD_SET(pr)) then print,result
	;print, result
end

; automate a function from miriad: statistics, compute rms
pro mirRMSpath, path, ra, dec, s, rms, freq, print=pr,import=import ;s in arcmin
	file=FILE_BASENAME(path)
	dir=FILE_DIRNAME(path)
	
	cd, dir

	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)


	mirAdxy,dir,mirfile,ra+s/60.,dec-s/60.,px0,py0
	mirAdxy,dir,mirfile,ra-s/60.,dec+s/60.,px1,py1
	;script="imfit in="+mirfile+" object=gaussian region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	script="imlist in=" + mirfile + " options=statistics region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"

	spawn, script, result,/STDERR
	;print,result

	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result

	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r

	rx=STREGEX(result, "plane"+spc+"Frequency"+spc+"Total Flux"+spc+"Maximum"+spc+"Minimum"+spc+"Average"+spc+"rms"+newline $
		+num+num+num+num+num+num+num,/SUBEXPR,/EXTRACT)

	if(rx[2] ne '')then freq=double(rx[2]) else freq=-99 ;GHz
	if(rx[7] ne '')then rms=double(rx[7]) else rms=-99 
	if(rms eq -1) then begin
		rx=STREGEX(result, spc+"Total Flux"+spc+"Maximum"+spc+"Minimum"+spc+"Average"+spc+"rms"+newline $
			+num+num+num+num+num,/SUBEXPR,/EXTRACT)
		if(rx[5] ne '')then rms=double(rx[5])else rms=-99 
	endif

	;print,'rms',rms
	;print,'freq',freq

	;print,rx

	;Axis:     3     FREQ
	;plane   Frequency   Total Flux   Maximum     Minimum     Average       rms
	;1   2.05205      5.8008E-03  1.4230E-03 -1.1656E-04  3.6926E-07  2.6937E-05

	;if(KEYWORD_SET(pr)) then print,script
	;if(KEYWORD_SET(pr)) then print,result

	;print, result
end

; automate a function from miriad: imfit, fit a 2D Gaussian to a source
pro mirImfit, dir, file, ra, dec, s, params, paramnames, print=pr, import=import

	params=[]
	paramnames=[]

	;if input file is in fits format (must have fits extension)
	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)

	cd, dir
	;print,'mirfile',mirfile
	
	;must escape .^$*+?()[{\|
	
	;regex blocks
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)
	
	;get radec in pixel coord
	mirAdxy,dir,mirfile,ra+s/3600.,dec-s/3600.,px0,py0
	mirAdxy,dir,mirfile,ra-s/3600.,dec+s/3600.,px1,py1
	
	;print,px0,py0,px1,py1
	
	script="imfit in="+mirfile+" object=gaussian region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r
	
	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result

	;rx=STREGEX(result, 'Fatal Error')
	;if rx gt 0 then return
	
	;rx=STREGEX(result, 'Nothing to fit')
	;if rx gt 0 then return
	
	rx=STREGEX(result, '  Peak value:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99 
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'peak']
	paramnames=[paramnames, 'peakerr']
	
	rx=STREGEX(result, '  Total integrated flux:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'flux']
	paramnames=[paramnames, 'fluxerr']
	
	rx=STREGEX(result, '  Major axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmaj']
	paramnames=[paramnames, 'bmajerr']
	
	rx=STREGEX(result, '  Minor axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmin']
	paramnames=[paramnames, 'bminerr']
	
	rx=STREGEX(result, '  Position angle \(degrees\): '+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bpa']
	paramnames=[paramnames, 'bpaerr']
	
	rx=STREGEX(result, '  Right Ascension:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valra=rx[1] else valra='-99'
	rx=STREGEX(result, '  Declination:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valdec=rx[1] else valdec='-99'
	if(not(valra eq '-99' and valdec eq '-99')) then $
		stringad, (valra +' '+valdec), valra,valdec
	params=[params,valra]
	params=[params,valdec]
	paramnames=[paramnames, 'ra']
	paramnames=[paramnames, 'dec']
	
	valdebpa=-99
	valdebmaj=-99
	valdebmin=-99
	
	rx=STREGEX(result, 'Deconvolution appears to produce a point source',/SUBEXPR,/EXTRACT)
	
	if(rx[0] eq '') then begin
	
		rx=STREGEX(result, 'Deconvolved Major, minor axes \(arcsec\):' $
			+num + num,/SUBEXPR,/EXTRACT)
		
		
		if(rx[1] ne '')then valdebmaj=double(rx[1]) else valdebmaj='-99'
		if(rx[2] ne '')then valdebmin=double(rx[2]) else valdebmin='-99'
		
		rx=STREGEX(result, 'Deconvolved Position angle \(degrees\):' $
			+num,/SUBEXPR,/EXTRACT)
		if(rx[1] ne '')then valdebpa=double(rx[1]) else valdebpa='-99'
				
	endif else begin
		valdebmaj=0
		valdebmin=0
		valdebpa=0
	endelse
	
	params=[params,valdebmaj]
	params=[params,valdebmin]
	params=[params,valdebpa]
	paramnames=[paramnames, 'bmajdeconv']
	paramnames=[paramnames, 'bmindeconv']
	paramnames=[paramnames, 'bpadeconv']

	;for i = 0L, N_ELEMENTS(params)-1 do begin
	;	print,paramnames[i],params[i]
	;endfor
	
	;dict=DICTIONARY(paramnames,params)
end

; automate a function from miriad
pro mirImfitPathClip, path, ra, dec, s, params, paramnames, clip, print=pr, import=import
	
	;result in Jy
	;clip in Jy
	
	params=[]
	paramnames=[]

	file=FILE_BASENAME(path)
	dir=FILE_DIRNAME(path)

	;if input file is in fits format (must have fits extension)
	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)


	cd, dir
	;print,'mirfile',mirfile

	;must escape .^$*+?()[{\|

	;regex blocks
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	;get radec in pixel coord
	mirAdxy,dir,mirfile,ra+s/3600.,dec-s/3600.,px0,py0
	mirAdxy,dir,mirfile,ra-s/3600.,dec+s/3600.,px1,py1

	;print,px0,py0,px1,py1

	script="imfit in="+mirfile+" object=gaussian clip="+str(clip)+" region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r


	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result


	;rx=STREGEX(result, 'Fatal Error')
	;if rx gt 0 then return

	;rx=STREGEX(result, 'Nothing to fit')
	;if rx gt 0 then return

	rx=STREGEX(result, '  Peak value:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99 
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'peak']
	paramnames=[paramnames, 'peakerr']

	rx=STREGEX(result, '  Total integrated flux:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'flux']
	paramnames=[paramnames, 'fluxerr']

	rx=STREGEX(result, '  Major axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmaj']
	paramnames=[paramnames, 'bmajerr']

	rx=STREGEX(result, '  Minor axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmin']
	paramnames=[paramnames, 'bminerr']

	rx=STREGEX(result, '  Position angle \(degrees\): '+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bpa']
	paramnames=[paramnames, 'bpaerr']

	rx=STREGEX(result, '  Right Ascension:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valra=rx[1] else valra='-99'
	rx=STREGEX(result, '  Declination:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valdec=rx[1] else valdec='-99'
	if(not(valra eq '-99' and valdec eq '-99')) then $
		stringad, (valra +' '+valdec), valra,valdec
	params=[params,valra]
	params=[params,valdec]
	paramnames=[paramnames, 'ra']
	paramnames=[paramnames, 'dec']

	valdebpa=-99
	valdebmaj=-99
	valdebmin=-99

	rx=STREGEX(result, 'Deconvolution appears to produce a point source',/SUBEXPR,/EXTRACT)

	if(rx[0] eq '') then begin

		rx=STREGEX(result, 'Deconvolved Major, minor axes \(arcsec\):' $
			+num + num,/SUBEXPR,/EXTRACT)


		if(rx[1] ne '')then valdebmaj=double(rx[1]) else valdebmaj='-99'
		if(rx[2] ne '')then valdebmin=double(rx[2]) else valdebmin='-99'

		rx=STREGEX(result, 'Deconvolved Position angle \(degrees\):' $
			+num,/SUBEXPR,/EXTRACT)
		if(rx[1] ne '')then valdebpa=double(rx[1]) else valdebpa='-99'

	endif else begin
		valdebmaj=0
		valdebmin=0
		valdebpa=0
	endelse

	params=[params,valdebmaj]
	params=[params,valdebmin]
	params=[params,valdebpa]
	paramnames=[paramnames, 'bmajdeconv']
	paramnames=[paramnames, 'bmindeconv']
	paramnames=[paramnames, 'bpadeconv']

	;for i = 0L, N_ELEMENTS(params)-1 do begin
	;	print,paramnames[i],params[i]
	;endfor

	;dict=DICTIONARY(paramnames,params)
end

; automate a function from miriad
pro mirImfitPath, path, ra, dec, s, params, paramnames, print=pr, import=import

	params=[]
	paramnames=[]

	file=FILE_BASENAME(path)
	dir=FILE_DIRNAME(path)
	
	;if input file is in fits format (must have fits extension)
	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)


	cd, dir
	;print,'mirfile',mirfile

	;must escape .^$*+?()[{\|

	;regex blocks
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	;get radec in pixel coord
	mirAdxy,dir,mirfile,ra+s/3600.,dec-s/3600.,px0,py0
	mirAdxy,dir,mirfile,ra-s/3600.,dec+s/3600.,px1,py1

	;print,px0,py0,px1,py1

	script="imfit in="+mirfile+" object=gaussian region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r


	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result


	;rx=STREGEX(result, 'Fatal Error')
	;if rx gt 0 then return

	;rx=STREGEX(result, 'Nothing to fit')
	;if rx gt 0 then return

	rx=STREGEX(result, '  Peak value:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'peak']
	paramnames=[paramnames, 'peakerr']

	rx=STREGEX(result, '  Total integrated flux:'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'flux']
	paramnames=[paramnames, 'fluxerr']

	rx=STREGEX(result, '  Major axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmaj']
	paramnames=[paramnames, 'bmajerr']

	rx=STREGEX(result, '  Minor axis \(arcsec\):'+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bmin']
	paramnames=[paramnames, 'bminerr']

	rx=STREGEX(result, '  Position angle \(degrees\): '+num+"("+pm+num+")?",/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then val=double(rx[1]) else val=-99
	if(rx[3] ne '' and rx[3] ne '*******')then err=double(rx[3]) else err=-99
	params=[params,val]
	params=[params,err]
	paramnames=[paramnames, 'bpa']
	paramnames=[paramnames, 'bpaerr']

	rx=STREGEX(result, '  Right Ascension:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valra=rx[1] else valra='-99'
	rx=STREGEX(result, '  Declination:'+num,/SUBEXPR,/EXTRACT)
	if(rx[1] ne '')then valdec=rx[1] else valdec='-99'
	if(not(valra eq '-99' and valdec eq '-99')) then $
		stringad, (valra +' '+valdec), valra,valdec
	params=[params,valra]
	params=[params,valdec]
	paramnames=[paramnames, 'ra']
	paramnames=[paramnames, 'dec']

	valdebpa=-99
	valdebmaj=-99
	valdebmin=-99

	rx=STREGEX(result, 'Deconvolution appears to produce a point source',/SUBEXPR,/EXTRACT)

	if(rx[0] eq '') then begin

		rx=STREGEX(result, 'Deconvolved Major, minor axes \(arcsec\):' $
			+num + num,/SUBEXPR,/EXTRACT)


		if(rx[1] ne '')then valdebmaj=double(rx[1]) else valdebmaj='-99'
		if(rx[2] ne '')then valdebmin=double(rx[2]) else valdebmin='-99'

		rx=STREGEX(result, 'Deconvolved Position angle \(degrees\):' $
			+num,/SUBEXPR,/EXTRACT)
		if(rx[1] ne '')then valdebpa=double(rx[1]) else valdebpa='-99'

	endif else begin
		valdebmaj=0
		valdebmin=0
		valdebpa=0
	endelse

	params=[params,valdebmaj]
	params=[params,valdebmin]
	params=[params,valdebpa]
	paramnames=[paramnames, 'bmajdeconv']
	paramnames=[paramnames, 'bmindeconv']
	paramnames=[paramnames, 'bpadeconv']

	;for i = 0L, N_ELEMENTS(params)-1 do begin
	;	print,paramnames[i],params[i]
	;endfor

	;dict=DICTIONARY(paramnames,params)
end

; +
; In miriad: Fit a 2D parabola to a source to find a maximum
; Usage:
; mirMaxfit,dir,file,ra,dec,s,params,pnames,/import;,/print
; raPeak=params[where(pnames eq 'ra')]
; decPeak=params[where(pnames eq 'dec')]
; -
pro mirMaxfit, dir, file, ra, dec, s, params, paramnames, print=pr, import=import
	params=[]
	paramnames=[]
	
	;if input file is in fits format (must have fits extension)
	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)
	
	
	cd, dir
	;print,'mirfile',mirfile
	
	;must escape .^$*+?()[{\|
	
	;regex blocks
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)
	
	;get radec in pixel coord
	mirAdxy,dir,mirfile,ra+s/3600.,dec-s/3600.,px0,py0
	mirAdxy,dir,mirfile,ra-s/3600.,dec+s/3600.,px1,py1
	
	;print,'px0,py0,px1,py1',px0,py0,px1,py1
	
	script="maxfit in="+mirfile+" region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r
	
	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result

	rx=STREGEX(result, 'Fatal Error')
	;if rx gt 0 then return
	
	rx=STREGEX(result, 'Nothing to fit')
	;if rx gt 0 then return

	;Peak pixel   : (5050,4202,1) = 5.1668E-05
	;Fitted pixel : (5048.78,4198.71,1) = 7.6129E-05
	
	rx=STREGEX(result, 'Peak pixel'+spc+': \('+num+','+num+'(,'+num+')?\) = '+num,/SUBEXPR,/EXTRACT)
	if(rx[5] ne '')then val=double(rx[5]) else val=-99
	params=[params,val]
	paramnames=[paramnames, 'peakmax']

	rx=STREGEX(result, 'Fitted pixel'+spc+': \('+num+','+num+'(,'+num+')?\) = '+num,/SUBEXPR,/EXTRACT)
	if(rx[5] ne '')then val=double(rx[5]) else val=-99
	params=[params,val]
	paramnames=[paramnames, 'peak']

	rx=STREGEX(result,"Axis 1: Fitted RA---SIN  =" + num + newline $
		+spc+"Axis 2: Fitted DEC--SIN  =" + num $
		,/SUBEXPR,/EXTRACT)
	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result,"Axis 1: Fitted RA---TAN  =" + num + newline $
		+spc+"Axis 2: Fitted DEC--TAN  =" + num $
		,/SUBEXPR,/EXTRACT)

	valra=-99
	valdec=-99
	if(rx[1] ne '' and rx[2] ne '') then stringad, rx[1]+' '+rx[2],valra,valdec
	params=[params,valra]
	paramnames=[paramnames, 'ra']
	params=[params,valdec]
	paramnames=[paramnames, 'dec']

	;print,'valra,valdec',valra,valdec

	if(0) then begin
		ramax=-1
		decmax=-1
		if(rx[1] ne '' and rx[2] ne '')then begin
		
			if(strmid(rx[1],0,1) ne '*' and strmid(rx[2],0,1) ne '*') then $
				mirXyad,dir,mirfile,double(rx[1]),double(rx[2]),ramax,decmax
		endif
		params=[params,ramax]
		paramnames=[paramnames, 'ramax']
		params=[params,decmax]
		paramnames=[paramnames, 'decmax']

		ramax=-1
		decmax=-1
		if(rx[1] ne '' and rx[2] ne '')then begin
			if(strmid(rx[1],0,1) ne '*' and strmid(rx[2],0,1) ne '*') then $
				mirXyad,dir,mirfile,double(rx[1]),double(rx[2]),ramax,decmax
		endif
		params=[params,ramax]
		paramnames=[paramnames, 'ra']
		params=[params,decmax]
		paramnames=[paramnames, 'dec']
	endif
	
	;for i = 0L, N_ELEMENTS(params)-1 do begin
	;	print,paramnames[i],params[i]
	;endfor
	;print,params
	;mirXyad,dir,mirfile,5048.78d,4198.71d,rafit,decfit
	;print,rafit,decfit
end

; automate a function from miriad
pro mirMaxfitPath, path, ra, dec, s, params, paramnames, print=pr, import=import
	;result in Jy

	params=[]
	paramnames=[]

	file=FILE_BASENAME(path)
	dir=FILE_DIRNAME(path)

	;if input file is in fits format (must have fits extension)
	mirfile=file
	if(KEYWORD_SET(import)) then mirfile=mirImport(dir, file, /check)

	cd, dir
	;print,'mirfile',mirfile
	;must escape .^$*+?()[{\|

	;regex blocks
	num="[ ]*([:\*0-9E\.\+\-]+)" ;number
	pm="[ ]*\+/\-" ;plus-minus sign
	spc="[ ]*" ;space
	if (!D.NAME eq 'WIN') then newline = string([13B, 10B]) else newline = string(10B)

	;get radec in pixel coord
	mirAdxy,dir,mirfile,ra+s/3600.,dec-s/3600.,px0,py0
	mirAdxy,dir,mirfile,ra-s/3600.,dec+s/3600.,px1,py1

	;print,'px0,py0,px1,py1',px0,py0,px1,py1

	script="maxfit in="+mirfile+" region='abspixel,box("+str(px0)+","+str(py0)+","+str(px1)+","+str(py1)+")'"
	spawn, script, result,/STDERR

	r=""
	for i = 0L, N_ELEMENTS(result)-1 do r+=result[i]+newline
	result=r

	if(KEYWORD_SET(pr)) then print,script
	if(KEYWORD_SET(pr)) then print,result

	rx=STREGEX(result, 'Fatal Error')
	;if rx gt 0 then return

	rx=STREGEX(result, 'Nothing to fit')
	;if rx gt 0 then return

	;Peak pixel   : (5050,4202,1) = 5.1668E-05
	;Fitted pixel : (5048.78,4198.71,1) = 7.6129E-05

	rx=STREGEX(result, 'Peak pixel'+spc+': \('+num+','+num+'(,'+num+')?\) = '+num,/SUBEXPR,/EXTRACT)
	if(rx[5] ne '')then val=double(rx[5]) else val=-99
	params=[params,val]
	;paramnames=[paramnames, 'peakmax']
	paramnames=[paramnames, 'peakpix']

	rx=STREGEX(result, 'Fitted pixel'+spc+': \('+num+','+num+'(,'+num+')?\) = '+num,/SUBEXPR,/EXTRACT)
	if(rx[5] ne '')then val=double(rx[5]) else val=-99
	params=[params,val]
	paramnames=[paramnames, 'peak']

	rx=STREGEX(result,"Axis 1: Fitted RA---SIN  =" + num + newline $
		+spc+"Axis 2: Fitted DEC--SIN  =" + num $
		,/SUBEXPR,/EXTRACT)
	if(rx[1] eq '' or rx[2] eq '') then $
		rx=STREGEX(result,"Axis 1: Fitted RA---TAN  =" + num + newline $
		+spc+"Axis 2: Fitted DEC--TAN  =" + num $
		,/SUBEXPR,/EXTRACT)

	valra=-99
	valdec=-99
	if(rx[1] ne '' and rx[2] ne '') then stringad, rx[1]+' '+rx[2],valra,valdec
	params=[params,valra]
	paramnames=[paramnames, 'ra']
	params=[params,valdec]
	paramnames=[paramnames, 'dec']

	;print,'valra,valdec',valra,valdec

	if(0) then begin
		ramax=-1
		decmax=-1
		if(rx[1] ne '' and rx[2] ne '')then begin

			if(strmid(rx[1],0,1) ne '*' and strmid(rx[2],0,1) ne '*') then $
				mirXyad,dir,mirfile,double(rx[1]),double(rx[2]),ramax,decmax
		endif
		params=[params,ramax]
		paramnames=[paramnames, 'ramax']
		params=[params,decmax]
		paramnames=[paramnames, 'decmax']

		ramax=-1
		decmax=-1
		if(rx[1] ne '' and rx[2] ne '')then begin
			if(strmid(rx[1],0,1) ne '*' and strmid(rx[2],0,1) ne '*') then $
				mirXyad,dir,mirfile,double(rx[1]),double(rx[2]),ramax,decmax
		endif
		params=[params,ramax]
		paramnames=[paramnames, 'ra']
		params=[params,decmax]
		paramnames=[paramnames, 'dec']
	endif

	;for i = 0L, N_ELEMENTS(params)-1 do begin
	;	print,paramnames[i],params[i]
	;endfor
	;print,params
	;mirXyad,dir,mirfile,5048.78d,4198.71d,rafit,decfit
	;print,rafit,decfit
end

