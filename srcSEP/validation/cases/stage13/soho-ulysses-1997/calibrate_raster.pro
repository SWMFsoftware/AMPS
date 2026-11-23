; SolarSoft CDS/NIS preprocessing for the archived May 27 rasters.
; Run in a CDS-enabled SolarSoft IDL session with the pointing, wave and
; VDS/NIS radiometry databases. Preserve the session log and database hashes.
pro calibrate_raster, input_path, output_prefix
  compile_opt idl2
  if file_test(output_prefix + '.sav') or file_test(output_prefix + '.wcs.csv') then $
    message, 'Calibration output already exists'

  ; READCDSFITS loads the instrument housekeeping and nominal wavelength
  ; calibration, then calls UPD_CDS_POINT with the supplied options. /ADJUST
  ; requests the documented helium-to-photospheric limb correction. The
  ; geometry downstream uses a photospheric radius from the start. Use the
  ; documented contiguous detector array format so the uncalibrated and
  ; calibrated pixels are both present in the saved QLDS product.
  error_text = ''
  qlds = readcdsfits(input_path, /batch, /adjust, /oldsys, status=error_text)
  if error_text ne '' then message, error_text
  if size(qlds, /type) ne 8 then message, 'READCDSFITS did not return a raster'

  ; GT_SOLAR_XY applies the NIS slit/band alignment for the specific Mg IX
  ; window. It accepts a zero-based window index. Geometry for a different
  ; band or line must be extracted from its own native window.
  index = where(strupcase(strtrim(qlds.detdesc.label, 2)) eq 'MG_9_368_07', count)
  if count ne 1 then message, 'Expected one Mg IX 368.07 A science window'
  gt_solar_xy, qlds, index[0], x, y, /nopadding, errmsg=error_text
  if error_text ne '' then message, error_text
  dimensions = size(x, /dimensions)
  if n_elements(dimensions) ne 3 then message, 'Expected wavelength and two sky axes'
  if dimensions[0] ne 21 or dimensions[1] ne 30 or dimensions[2] ne 143 or $
    n_elements(y) ne n_elements(x) then message, 'Unexpected CDS raster shape'
  if min(finite(x)) eq 0 or min(finite(y)) eq 0 then $
    message, 'Corrected WCS contains invalid positions'

  ; Keep the original detector values and descriptors, including missing and
  ; saturated samples, for auditing before and after SolarSoft despiking.
  ; NIS_CALIB delegates CCD debias, flat, detector throughput and time
  ; corrections to VDS_CALIB, and applies the official date-dependent
  ; telescope response. Radiance is per spectral pixel and still needs a
  ; line-profile fit with its own quality mask and covariance.
  original_adc = qlds.detdata
  original_windows = qlds.detdesc
  cds_new_spike, qlds
  nis_calib, qlds, /ergs, /steradians, errmsg=error_text
  if error_text ne '' then message, error_text
  units = strtrim(qlds.detdesc[index[0]].units, 2)
  if units ne 'ERGS/CM^2/SEC/STERAD' then message, 'Unexpected NIS radiance units'

  ; Publish only after the instrument preprocessing and unit checks succeed.
  ; A paired package contains the corrected positions and the whole QLDS so
  ; another analyst can inspect the measured signal and fitted provenance.
  runtime_version = !version.release
  save, qlds, original_adc, original_windows, runtime_version, input_path, x, y, $
    filename=output_prefix + '.sav.partial'
  openw, unit, output_prefix + '.wcs.csv.partial', /get_lun
  printf, unit, 'x_index_zero,y_index_zero,x_arcsec,y_arcsec'
  ; The correction varies slightly along the wavelength axis because the
  ; slit image is tilted. The magnetic screen uses the nominal central Mg IX
  ; sample. A line fit can use each spectral pixel's own WCS from the save.
  center_wave_index = dimensions[0] / 2L
  for ix = 0L, dimensions[1] - 1 do begin
    for iy = 0L, dimensions[2] - 1 do begin
      printf, unit, ix, iy, x[center_wave_index,ix,iy], y[center_wave_index,ix,iy], $
        format='(I0,",",I0,",",F14.6,",",F14.6)'
    endfor
  endfor
  free_lun, unit
  file_move, output_prefix + '.sav.partial', output_prefix + '.sav'
  file_move, output_prefix + '.wcs.csv.partial', output_prefix + '.wcs.csv'
end
