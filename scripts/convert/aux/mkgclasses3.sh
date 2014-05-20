#! /usr/bin/env bash

do_mkclass="yes"

for fbase in g3numi g4numi flugg g4minerva
do
  tname="xyzzy"
  if [ "$fbase" = "g3numi" ] ; then tname="h10" ; fi
  if [ "$fbase" = "flugg"  ] ; then tname="h10" ; fi
  if [ "$fbase" = "g4numi" ] ; then tname="nudata" ; fi
  if [ "$fbase" = "g4minerva" ] ; then tname="nudata" ; fi

  if [ "yes" = "${do_mkclass}" ] ; then
    cat > make_generic.C <<EOF
#include <string>
#include <vector>
#include <map>
void make_generic()
{
  TFile* g3file = TFile::Open("generic_${fbase}.root","READONLY");
  gROOT->ProcessLine("${tname}->MakeClass();");
}
EOF
    if [ ! -f generic_${fbase}.root ]; then
      echo "skip $fbase for lack of generic_${fbase}.root"
      continue
    fi
    root -b -q make_generic.C
    rm -f make_generic.C
  fi

  if [ ! -d ${fbase} ]; then mkdir ${fbase}; fi
  
  fnameh=${fbase}/${fbase}.h
  fnameC=${fbase}/${fbase}.C
  rm -f ${fnameh} ${fnameC}
  touch ${fnameh} ${fnameC}

  cat ${tname}.h | sed -e "s/${tname}/${fbase}/g" \
                       -e "s%Cut(Long64_t entry)%Cut(Long64_t /* entry */)%" \
                   >> $fnameh
  cat ${tname}.C | sed -e "s/${tname}/${fbase}/g" >> $fnameC

  rm -f ${tname}.h ${tname}.C

done


# End-of-Script
