from ROOT import TLorentzVector

class Lepton:

    """a wrapper for a lepton object (p4, id, charge and isolation variables)"""

    def __init__(self,pdgId,pt,eta,phi,m,charge):
        self.pdgId=pdgId
        self.p4=TLorentzVector(0,0,0,0)
        self.p4.SetPtEtaPhiM(pt,eta,phi,m)
        self.isEE = True if abs(eta)>1.4442 else False
        self.charge=charge
        self.chiso=0
        self.nhiso=0
        self.phiso=0
        self.iso=0
        self.miniiso=0
        self.rho={}

    def addIsoComponents(self,chiso,nhiso,phiso,miniiso=0):
        self.chiso=chiso
        self.nhiso=nhiso
        self.phiso=phiso
        self.iso=chiso+nhiso+phiso
        self.miniiso=miniiso

    def addRho(self,key,val):
        self.rho[key]=val

    def getIsolation(self, doMiniIso=False):
        if not 'rho' in self.rho : return -1

        rhoVal=self.rho['rho']
        isoVal=0

        if doMiniIso:
            if abs(self.pdgId)==11:
                ue=-0.000097*((rhoVal+12.3874)**2)+0.10686*(rhoVal+12.3874)
            else:
                ue=0.000165*((rhoVal-5.0441)**2)+0.211759*(rhoVal-5.0441)
            isoVal=(self.miniiso*self.p4.Pt()-ue)/self.p4.Pt()

        else:
            if abs(self.pdgId)==11:
                ue=0.001020*((rhoVal+23.421)**2)-0.1482*(rhoVal+23.421)            
            else:
                ue=0.001285*((rhoVal+15.829)**2)+0.294*(rhoVal+15.829)
                
            isoVal=(self.chiso+self.nhiso+self.phiso-ue)/self.p4.Pt()

        return isoVal

    def isIsolated(self,doMiniIso=False):

        isoThr=999.
        if abs(self.pdgId)==11:
            isoThr=0.3 if doMiniIso else 0.3
        else:
            isoThr=0.0 if doMiniIso else 0.4
        return self.getIsolation(doMiniIso)<isoThr

class Dilepton:

    """a wrapper for a dilepton object"""

    def __init__(self,l1,l2,isOF,isSS,isZ,isIso,isMixed=False):
        self.l1=l1 if l1.p4.Pt()>l2.p4.Pt() else l2
        self.l2=l2 if l1.p4.Pt()>l2.p4.Pt() else l1
        self.p4=self.l1.p4+self.l2.p4
        self.flavour=abs(l1.pdgId*l2.pdgId)
        self.isOF=isOF
        self.isSS=isSS
        self.isZ=isZ
        self.isIso=isIso
        self.isMixed=isMixed

def dileptonBuilder(lepColl,doMiniIsolation=False):

    """takes the leading lepton pair and sets some quality flags"""

    if len(lepColl)<2: raise ValueError('Less than 2 leptons')

    isOF = True if abs(lepColl[0].pdgId*lepColl[1].pdgId)==143                             else False
    isSS = True if lepColl[0].charge*lepColl[1].charge>0                                   else False
    isZ  = True if abs((lepColl[0].p4+lepColl[1].p4).M()-91.)<15 and not isSS and not isOF else False
    isIso= True if lepColl[0].isIsolated(doMiniIsolation) and lepColl[1].isIsolated(doMiniIsolation) else False

    return Dilepton(lepColl[0],lepColl[1],isOF,isSS,isZ,isIso)


def getLeptons(t,pdgIdColl=[13,11]):
    """ get all the electrons in the event """

    lepColl=[]
    for il in range(t.nlep):
        absid=abs(t.lep_pdgId[il])
        if not absid in pdgIdColl: continue
        mass=0.511e-3 if absid==11 else 0.105658
        lepColl.append( Lepton(t.lep_pdgId[il],
                               t.lep_pt[il],
                               t.lep_eta[il],
                               t.lep_phi[il],
                               mass,
                               t.lep_charge[il]) )
        lepColl[-1].addIsoComponents(t.lep_chiso[il],t.lep_nhiso[il],t.lep_phiso[il],t.lep_miniiso[il]*t.lep_pt[il])

        #this is hardcoded, as in the rho tree producer...
        eta=t.lep_eta[il]
        eta_idx=None
        if eta>-3.0 and eta<=-2.1 : eta_idx=1
        if eta>-2.1 and eta<=-1.3 : eta_idx=2
        if eta>-1.3 and eta<=1.3  : eta_idx=3
        if eta>1.3  and eta<=2.1  : eta_idx=4
        if eta>2.1  and eta<=3.0  : eta_idx=5
        
        try:
            for key in ['lep_rho']:
                lepColl[-1].addRho(key,getattr(t,key)[il])
            for key in ['chrho','nhrho','phorho']:
                lepColl[-1].addRho(key,getattr(t,key)[0])
                lepColl[-1].addRho(key+'_restr',getattr(t,key)[2 if abs(eta)<1.4442 else 1])
            lepColl[-1].addRho('rho',t.rho[eta_idx] if eta_idx else 0)
        except:
            pass

    return lepColl


def getDilepton(t,pdgIdColl=[13,11],doMiniIsolation=False):

    """get the dilepton in the event"""
    
    return dileptonBuilder( getLeptons(t,pdgIdColl),doMiniIsolation )
