#!/bin/bash

> processlist_mergir_missing

# Loop over year
for year in {2000..2013}; do
  # Loop over 12 months
  for mon in {01..12}; do
    echo run_mergir_missing_data.sh ${year} ${mon} >> processlist_mergir_missing
  done
done

