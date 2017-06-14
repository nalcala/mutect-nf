#!/bin/bash
cd ~/mutect-nf-nf/
git config --global user.email "robitaillea@students.iarc.fr"
git add dag.png
git add dag.html
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/mutect-nf/trigger/848d00e1-4e32-498a-aaaa-4e17029cdf3c/
