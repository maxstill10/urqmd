#!/bin/bash

for ((i=12; i<25; i++)); do
	mkdir /gpfs01/star/scratch/mmorozov/urqmd_3.9_Zheni/log_cms/subdata${i}
	mkdir /star/data01/pwg/mmorozov/urqmd/oo/6gev/subdata${i}

	OutFile="run_urqmd_$i.xml"

	cat > "$OutFile" <<EOF	
<?xml version="1.0" encoding="utf-8" ?>
<job nProcesses="10000">
<shell>singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif</shell>
    <command>
      setenv NODEBUG yes
      stardev

        mkdir out-f13
        setenv SEED `head -c 3 /dev/urandom | xxd -p | sed 's;^;0x;g' | awk --non-decimal-data '{printf "%d\n", $0}'`
        sed -i "s;RNDNUM;\${SEED};g" inputfile_oo_6.urqmd

        # UrQMD environment variables
        setenv ftn09 "inputfile_oo_6.urqmd"
        setenv ftn13 "\${JOBID}.f13"
        # urqmd2mc environment variables
        setenv MCDST `dirname /star/u/mmorozov/urqmd/urqmd_Zheni/McDst/urqmd2mc;`
        setenv LD_LIBRARY_PATH "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/:\$LD_LIBRARY_PATH"

        time ./urqmd.x86_64
        setenv outfile `echo \$ftn13 | sed "s;\.f13;.mcDst.root;g"`
        time ./urqmd2mc "\$ftn13" 1000
  </command>

  <SandBox>
    <Package>
        <File>file:/star/u/mmorozov/urqmd/urqmd_Zheni/urqmd.x86_64</File>
        <File>file:/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/urqmd2mc</File>
        <File>file:/star/u/mmorozov/urqmd/inputfile_oo_6.urqmd</File>
        <File>file:/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/libMcDst.so</File>
    </Package>
  </SandBox>

  <stdout URL="file:/gpfs01/star/scratch/mmorozov/urqmd_3.9_Zheni/log_cms/subdata${i}/\$JOBID.out" />
  <stderr URL="file:/gpfs01/star/scratch/mmorozov/urqmd_3.9_Zheni/log_cms/subdata${i}/\$JOBID.err" />
  <output fromScratch="./\$JOBID.mcDst.root" toURL="file:/star/data01/pwg/mmorozov/urqmd/oo/6gev/subdata${i}"/>

  <Generator>
      <Location>/gpfs01/star/scratch/mmorozov/urqmd_3.9_Zheni/log_cms/subdata${i}/</Location>
  </Generator>
</job>

EOF
	star-submit-beta $OutFile
done
