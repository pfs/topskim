from ROOT import TLorentzVector


LEPTONBRANCHES=['pt','eta','phi',
                'idflags','d0','d0err','dz',
                'phiso','chiso','nhiso','rho',
                'pdgId','charge',
                'isofull','isofull20','isofull25','isofull30','miniiso'] 

DILEPTONBRANCHES=['llpt','lleta','llphi','llm','dphi','deta','sumeta','apt','bdt','bdtrarity']

class PhysicsObject:

    """a container for a physics object, mostly a passive object waiting for properties to be added from the source"""

    def __init__(self,tag='lepton'):
        self.tag=tag
        self.p4=TLorentzVector(0,0,0,0)
        self.matched=False

    def addProperty(self,name,val):
        setattr(self,name,val)

    def buildP4(self):
        try:
            self.p4.SetPtEtaPhiM(self.pt,self.eta,self.phi,self.m)
        except Exception as e:
            print '<'*50
            print 'Unable to set p4'
            print e
            print '<'*50

    def setGlobalEventProperties(self,propDict):
        for key in propDict:
            setattr(self,key,propDict[key])

    def getIsolation(self, isoType):
        return 0

    def isIsolated(self,isoType):
        return True

class DileptonObject:

    """a wrapper for a dilepton object"""

    def __init__(self,l1,l2,isOF,isSS,isZ):
        self.l1=l1 if l1.p4.Pt()>l2.p4.Pt() else l2
        self.l2=l2 if l1.p4.Pt()>l2.p4.Pt() else l1
        self.p4=self.l1.p4+self.l2.p4
        self.flavour=abs(l1.pdgId*l2.pdgId)
        self.isOF=isOF
        self.isSS=isSS
        self.isZ=isZ

def getLeptons(t,pdgIdColl=[13,11]):

    """ get all the electrons in the event """

    lepColl=[]

    for il in range(t.nlep):
        absid=abs(t.lep_pdgId[il])
        if not absid in pdgIdColl: continue
        mass=0.511e-3 if absid==11 else 0.105658
        lepColl.append( PhysicsObject() )

        for name in LEPTONBRANCHES:
            lepColl[-1].addProperty(name,getattr(t,'lep_'+name)[il])
        lepColl[-1].addProperty('m',mass)            
        lepColl[-1].matched=bool(t.lep_matched[il])
        lepColl[-1].buildP4()
        lepColl[-1].setGlobalEventProperties({'cenbin':t.cenbin,
                                              'ncoll':t.ncoll})
        

    return lepColl


def dileptonBuilder(lepColl):

    """takes the leading lepton pair and sets some quality flags"""

    if len(lepColl)<2: raise ValueError('Less than 2 leptons')

    isOF = True if abs(lepColl[0].pdgId*lepColl[1].pdgId)==143                             else False
    isSS = True if lepColl[0].charge*lepColl[1].charge>0                                   else False
    isZ  = True if abs((lepColl[0].p4+lepColl[1].p4).M()-91.)<15 and not isSS and not isOF else False

    return DileptonObject(lepColl[0],lepColl[1],isOF,isSS,isZ)

def getDilepton(t,pdgIdColl=[13,11]):

    """get the dilepton in the event"""
    
    return dileptonBuilder( getLeptons(t,pdgIdColl) )
