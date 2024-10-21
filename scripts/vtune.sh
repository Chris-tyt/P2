vtune -collect threading -result-dir ./vtune_results ./sph.x

vtune -report hotspots -result-dir ./vtune_results -format text -report-output ./hotspots_report_omp.txt