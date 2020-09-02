import SM 

qcd_htbins = [300, 500, 700, 1000, 1500, 2000, 'Inf']
qcd_htbins = map(str, qcd_htbins)
vjets_htbins = [100, 200, 400, 600, 800, 1200, 2500, 'Inf']
vjets_htbins = map(str, vjets_htbins)

variations = ['ISRRenorm', 'FSRRenorm', 'ISRFact', 'FSRFact', 'MPI']
_variations = ['Nominal']
_variations += [x + 'Up' for x in variations]
_variations += [x + 'Down' for x in variations]
higgs_datasets = [
    'ZH_HToBB_ZToNuNu_M125',
    'ZH_HToBB_ZToLL_M125',
    'ggZH_HToBB_ZToNuNu_M125',
    'ggZH_HToBB_ZToLL_M125',
    #'ggZH_HToBB_ZToQQ_M125',
    'WminusH_HToBB_WToLNu_M125',
    'WplusH_HToBB_WToLNu_M125',
]


processes = {}
'''

for bins in zip(qcd_htbins[:-1],qcd_htbins[1:]):
    name = 'QCD_HT%sto%s'%bins 
    xsec = SM.processes[name+'_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'][-1]
    for v in _variations:
        vname = name + '_' + v
        processes[vname] = (vname, 'MC', xsec)
for bins in zip(vjets_htbins[:-1],vjets_htbins[1:]):
    name = 'WJetsToLNu_HT-%sTo%s'%bins 
    xsec = SM.processes[name+'_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'][-1]
    for v in _variations:
        vname = name + '_' + v
        processes[vname] = (vname, 'MC', xsec)
    name = 'ZJetsToNuNu_HT-%sTo%s'%bins 
    xsec = SM.processes[name+'_13TeV-madgraph'][-1]
    for v in _variations:
        vname = name + '_' + v
        processes[vname] = (vname, 'MC', xsec)

for dataset in higgs_datasets:
    xsec = SM.processes[dataset+'_13TeV_powheg_pythia8'][-1];
    for v in _variations:
        vname = dataset + '_' + v
        processes[vname] = (vname, 'MC', xsec)
'''
