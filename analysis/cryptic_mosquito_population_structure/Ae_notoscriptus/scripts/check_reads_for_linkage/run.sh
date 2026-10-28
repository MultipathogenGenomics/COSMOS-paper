#!/bin/bash

for chunk in {0..9}
do
    echo "Running chunk $chunk..."
    python -m analysis.Ae_notoscriptus.check_reads_for_linkage.main $chunk &
done

wait
echo "All chunks completed"