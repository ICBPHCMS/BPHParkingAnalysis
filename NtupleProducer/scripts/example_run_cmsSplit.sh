##with input a DASpath in case of cmsRun 
python ~/cmsSplit.py --anType cmsRun --cfg test94X_NANO_template.py  --tag BToKJPsiee  --dasquery  --das /Bu_KJPsi_ee_Pythia/tstreble-BuToKJPsiee_Pythia_MINIAODSIM_18_06_05-393c5b9eb90ffa47a3c9a0f5562ec979/USER --filesperjob 2  --storeArea YOUR_OUTPUT_FOLDER


##with input txt in case of simple script
python ~/cmsSplit.py --anType script --cfg scripts/config_data_runB.sh --tag BToKee_2018B_BPH1 --listquery -i list_nanoDATA_BToKee_2018B_BPH1.txt --filesperjob 1 --storeArea YOUR_OUTPUT_FOLDER


##with input txt in case of simple script and JOBID needed as input
python ~/cmsSplit.py --anType scriptAndJOBID --cfg ../scripts/config_runAnalyis.sh --tag BToKmumu_2018B_BPH4_NN_BkgR_OS --listquery -i list_ntuDATA_BToKmumu_2018B_BPH4.txt --filesperjob 20 --storeArea YOUR_OUTPUT_FOLDER
