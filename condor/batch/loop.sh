#!/bin/bash

# List of volumes
volumes=(
  1mKplate 10mKplate HSplate MCplate HEplate Stillplate Potplate IVCplate NSplate NuclearStage
  110PillarUnion 10MCPillarUnion HSMCPillarUnion NS1PillarUnion PillarNSHS
  Dewar
  Topplate
  Floor Ceiling PbBlocks AlFrame AlHSs PbShot
  RadShield IVCplate PlasticShield BaffleUnion
  SampleMagnet1 MagMount1 SampleMagnet2 MagMount2 SM1PillarUnion SM2PillarUnion
  Cell BaseScrewUnion Holder HeatEx Silver Shell
  SourceMount Window Cover
  BoloWall Lhe3detector BaseScrewUnion Base CopperFoil
)

volumes=(
  CopperFoil
)

# List of isotopes
isotopes=(U238 Ra226 Pb210 Th232 Th228 U235 Cs137 K40 Co60 Mn54)

# Loop over each volume and each isotope
for volume in "${volumes[@]}"; do
    for isotope in "${isotopes[@]}"; do
	echo ""
	echo "Processing: $volume, $isotope"
	./count.sh $volume-$isotope  | grep  -e 'Root files completed' -e 'Jobs completed'
    done
done
