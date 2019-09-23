base=/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar;

a=(#datacards_2019-09-19_elPt25muPt20_ttbarBRFix
    datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis
   # datacards_2019-09-20_elPt25muPt25_ttbarBRFix
    datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis
)

for i in ${a[@]}; do
    #copy over
    cp -r ${base}/${i} ./

    #combine cards
    a=(`ls ${i}`)
    for j in ${a[@]}; do
        ana=${i}/${j}
        cd ${ana}
        if [[ ${i} == *"jetAnalysis"* ]]; then
            echo "Combining cards for jet analysis @ ${ana}" 
            combineCards.py ee0b=ee_0b.card.txt ee1b=ee_1b.card.txt ee2b=ee_2b.card.txt mm0b=mm_0b.card.txt mm1b=mm_1b.card.txt mm2b=mm_2b.card.txt > sameFlavor.card.txt
            combineCards.py em0b=em_0b.card.txt em1b=em_1b.card.txt em2b=em_2b.card.txt > em.card.txt
            combineCards.py ee0b=ee_0b.card.txt ee1b=ee_1b.card.txt ee2b=ee_2b.card.txt mm0b=mm_0b.card.txt mm1b=mm_1b.card.txt mm2b=mm_2b.card.txt em0b=em_0b.card.txt em1b=em_1b.card.txt em2b=em_2b.card.txt > allFlavors.card.txt        
        else
            echo "Combining cards for inclusive analysis @ ${ana}"
            combineCards.py ee=ee.card.txt mm=mm.card.txt > sameFlavor.card.txt       
            combineCards.py ee=ee.card.txt mm=mm.card.txt em=em.card.txt > allFlavors.card.txt       
        fi
        cd -
    done

    a=(`find ${i} | grep .txt`)
    for j in ${a[@]}; do

        if [[ ${j} == *".tmp"* || ${j} == *".log"*  ]]; then
            continue
        fi
        if [[ ${j} == *"em.card"* || ${j} == *"sameFlavor.card"* || ${j} == *"allFlavors.card"*  ]]; then

            echo ${j};
            sed -i 's/zg_lnN/* autoMCStats 0 1 1\nmuzg rateParam * zg 1\n#zg_lnN/g' ${j};
            sed -i 's/kmax/kmax *\n#kmax/g' ${j};        
            sh runUnblind.sh ${j}  > ${j}.log &
        fi
    done    
done


