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
        self.rho={}

    def addIsoComponents(self,chiso,nhiso,phiso):
        self.chiso=chiso
        self.nhiso=nhiso
        self.phiso=phiso
        self.iso=chiso+nhiso+phiso

    def addRho(self,key,val):
        self.rho[key]=val

    def getIsolation(self):
        if not 'rho' in self.rho : return -1

        rhoVal=self.rho['rho']
        if abs(self.pdgId)==11:
            ue=0.0011*((rhoVal+142.4)**2)-0.14*(rhoVal+142.4)
        else:
            ue=0.0013*((rhoVal+15.8)**2)+0.29*(rhoVal+15.8)

        return (self.chiso+self.nhiso+self.phiso-ue)/self.p4.Pt()

    def isIsolated(self):
        if abs(self.pdgId)==11:
            return self.getIsolation()<0.16
        else:
            return self.getIsolation()<0.26

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

def dileptonBuilder(lepColl):

    """takes the leading lepton pair and sets some quality flags"""

    if len(lepColl)<2: raise ValueError('Less than 2 leptons')

    isOF = True if abs(lepColl[0].pdgId*lepColl[1].pdgId)==143                             else False
    isSS = True if lepColl[0].charge*lepColl[1].charge>0                                   else False
    isZ  = True if abs((lepColl[0].p4+lepColl[1].p4).M()-91.)<15 and not isSS and not isOF else False
    isIso= True if lepColl[0].isIsolated() and lepColl[1].isIsolated()                     else False

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
        lepColl[-1].addIsoComponents(t.lep_chiso[il],t.lep_nhiso[il],t.lep_phiso[il])

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


def getDilepton(t,pdgIdColl=[13,11]):

    """get the dilepton in the event"""
    
    return dileptonBuilder( getLeptons(t,pdgIdColl) )
