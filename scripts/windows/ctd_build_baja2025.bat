@echo off
REM Convenience wrapper for the baja2025 cruise. New cruises should use
REM ctd_build.bat directly:  ctd_build.bat ^<cruise_id^>
call "%~dp0ctd_build.bat" baja2025 --xmlcon cruises\baja2025\calibration.xmlcon --hex-dir cruises\baja2025\hex %*
