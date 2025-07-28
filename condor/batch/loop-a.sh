#!/bin/bash

# List of volumes
volumes=(
  Dewar
  1mKplate 10mKplate HSplate MCplate HEplate Stillplate Potplate IVCplate NSplate NuclearStage
  110PillarUnion 10MCPillarUnion HSMCPillarUnion NS1PillarUnion PillarNSHS
  Topplate
  Floor Ceiling PbBlocks AlFrame AlHSs PbShot
  RadShield IVCplate PlasticShield BaffleUnion
)

volumes=(
  HSplate MCplate HEplate Stillplate Potplate IVCplate NSplate NuclearStage
  110PillarUnion 10MCPillarUnion HSMCPillarUnion NS1PillarUnion PillarNSHS
  Topplate
  Floor Ceiling PbBlocks AlFrame AlHSs PbShot
  RadShield IVCplate PlasticShield BaffleUnion
)



# List of isotopes
isotopes=(U238 Ra226 Pb210 Th232 Th228 U235 Cs137 K40 Co60 Mn54)

# Loop over each volume and each isotope
for volume in "${volumes[@]}"; do
    for isotope in "${isotopes[@]}"; do
	echo ""
	echo "Processing: $volume, $isotope"
	./count-a.sh $volume-$isotope  | grep  -e 'Root files completed' -e 'Jobs completed'
    done
done
