#!/bin/bash
# Convenience wrapper for the baja2025 cruise. New cruises should use
# ctd_build.sh directly:  ./ctd_build.sh <cruise_id>
exec "$(dirname "$0")/ctd_build.sh" baja2025 "$@"
