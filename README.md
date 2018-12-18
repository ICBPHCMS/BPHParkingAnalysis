# BPHParkingAnalysis
Analysis package for analysis of BPH parked data

### Instructions for 10_2_0
```
cmsrel CMSSW_10_2_0
cd CMSSW_10_2_0/src/
cmsenv

git clone https://github.com/ICBPHCMS/BPHParkingAnalysis.git
scram b
```

### Producing ntuples
```
BToKpipiNtupleProducer mc --input test94X_NANO.root --output test.root
```

### Evaluating BMX neutral network discriminant (only Kmumu so far)
```
BToKmumuNtupleProducer mc --input test94X_NANO.root --output test.root --addNNBMX BPHParkingAnalysis/NtupleProducer/data/weights_50.pb
```
