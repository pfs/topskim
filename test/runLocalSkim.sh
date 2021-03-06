#!/bin/bash

#
# example
# sh /afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src/topskim/test/runLocalSkim.sh \
#    -i /store/himc/HINPbPbWinter16DR/Pythia6_bJet170_pp502_Hydjet_MB/AODSIM/75X_mcRun2_HeavyIon_v13-v1/50000/00A14813-B70E-E611-8C4E-02163E014702.root \
#    -o hiftest.root \
#    -s /tmp/psilva
#

while getopts "i:o:s:h" opt; do
  case $opt in
      i)
	  input=$OPTARG
	  ;;
      o)
	  outfile=$OPTARG
	  ;;
      s)
	  outdir=$OPTARG
	  ;;
      h)
	  echo "Example"
	  echo "sh /afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src/topskim/test/runLocalSkim.sh \ "
	  echo "    -i /store/himc/HINPbPbWinter16DR/Pythia6_bJet170_pp502_Hydjet_MB/AODSIM/75X_mcRun2_HeavyIon_v13-v1/50000/00A14813-B70E-E611-8C4E-02163E014702.root\ "
	  echo "    -o hiftest.root -s /tmp/psilva"
	  exit -1;
	  ;;
  esac
done

export X509_USER_PROXY=~/private/cur_proxy
WORKDIR=`pwd`

echo "Working directory is ${WORKDIR}"
cd /afs/cern.ch/user/p/psilva/scratch0/CMSSW_7_5_8_patch3/src
eval `scram r -sh`
cd ${WORKDIR}
echo "Forestizing ${input}"
cmsRun ${CMSSW_BASE}/src/topskim/test/runOverPbPb_MIX_75X.py inputFile=${input} outFilename=hif_${outfile}
echo "Making jet tree out of HIForest"
root -b -q ${CMSSW_BASE}/src/topskim/makeJetsSkim.C+\(\"${outfile}\",\"hif_${outfile}\"\) && .q;
echo "Moving result to ${outdir}"
rm hif_${outfile}
xrdcp ${outfile} root://eoscms//eos/cms/${outdir}/${outfile}