#!/bin/bash

for chunk in {0..9}
do
    echo "Running chunk $chunk..."
    python -m analysis.cryptic_mosquito_population_structure.Ae_notoscriptus.scripts.check_reads_for_linkage.main $chunk &
done

wait
echo "All chunks completed"