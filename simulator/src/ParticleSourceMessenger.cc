#include <G4Geantino.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithABool.hh>
#include <G4Tokenizer.hh>
#include <G4ios.hh>
#include <G4SystemOfUnits.hh>

#include "ParticleSource.hh"
#include "ParticleSourceMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** commands concerning the particle source supplied to UI
 *
 */
ParticleSourceMessenger::ParticleSourceMessenger(ParticleSource *particleSource)
    : m_particleSource(particleSource)
{
    G4UIparameter *param;

    m_particleTable = G4ParticleTable::GetParticleTable();

    // create directory
    m_uiDirectory = new G4UIdirectory("/G4QSim/gun/");
    m_uiDirectory->SetGuidance("Particle Source control commands.");

    // list available particles
    m_listCmd = new G4UIcmdWithoutParameter("/G4QSim/gun/list", this);
    m_listCmd->SetGuidance("List available particles.");
    m_listCmd->SetGuidance("Invoke G4ParticleTable.");

    m_typeCmd = new G4UIcmdWithAString("/G4QSim/gun/type", this);
    m_typeCmd->SetGuidance("Set source type.");
    m_typeCmd->SetGuidance("Can be [`betaDecay', `doubleBetaDecay'].");
    m_typeCmd->SetParameterName("gunType", true, true);
    m_typeCmd->SetDefaultValue("betaDecay");
    m_typeCmd->SetCandidates("betaDecay doubleBetaDecay");

    m_QValueCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/gun/QValue", this);
    m_QValueCmd->SetGuidance("Q Value of beta decay.");
    m_QValueCmd->SetParameterName("QValue", true, true);
    m_QValueCmd->SetDefaultUnit("MeV");

    m_EHistCmd = new G4UIcmdWithAString("/G4QSim/gun/EHist", this);
    m_EHistCmd->SetGuidance("Load beta energy spectrum from a file.");
    m_EHistCmd->SetGuidance("E [MeV] versus probability.");
    m_EHistCmd->SetParameterName("EHistFName", true, true);
    m_EHistCmd->SetDefaultValue("EHist.dat");

    m_dBDEventsCmd = new G4UIcmdWithAString("/G4QSim/gun/dBDEvents", this);
    m_dBDEventsCmd->SetGuidance("Load double beta decay from a file.");
    m_dBDEventsCmd->SetGuidance("The file is generated by DECAY0.");
    m_dBDEventsCmd->SetParameterName("dBDEventsFName", true, true);
    m_dBDEventsCmd->SetDefaultValue("decay0.dat");

    // father ion
    m_fIonPDirCmd = new G4UIcmdWith3Vector("/G4QSim/gun/fIonPDir", this);
    m_fIonPDirCmd->SetGuidance("Set momentum direction of father ion.");
    m_fIonPDirCmd->SetParameterName("PX", "PY", "PZ", true, true);

    m_fIonPosCmd = new G4UIcmdWith3VectorAndUnit("/G4QSim/gun/fIonPos", this);
    m_fIonPosCmd->SetGuidance("Set starting position of father ion.");
    m_fIonPosCmd->SetParameterName("X", "Y", "Z", true, true);
    m_fIonPosCmd->SetDefaultUnit("cm");

    m_fIonEkCmd = new G4UIcmdWithADoubleAndUnit("/G4QSim/gun/fIonEk", this);
    m_fIonEkCmd->SetGuidance("Set father ion kinetic energy.");
    m_fIonEkCmd->SetParameterName("fIonEk", true, true);
    m_fIonEkCmd->SetDefaultUnit("keV");
    m_fIonEkCmd->SetDefaultValue(0.0);

    m_fIonCmd = new G4UIcommand("/G4QSim/gun/fIon", this);
    m_fIonCmd->SetGuidance("Set properties of father ion.");
    m_fIonCmd->SetGuidance("[usage] /gun/fIon Z A Q E");
    m_fIonCmd->SetGuidance("        Z:(int) AtomicNumber");
    m_fIonCmd->SetGuidance("        A:(int) AtomicMass");
    m_fIonCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
    m_fIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

    param = new G4UIparameter("Z", 'i', false);
    param->SetDefaultValue("1");
    m_fIonCmd->SetParameter(param);
    param = new G4UIparameter("A", 'i', false);
    param->SetDefaultValue("1");
    m_fIonCmd->SetParameter(param);
    param = new G4UIparameter("Q", 'i', true);
    param->SetDefaultValue("0");
    m_fIonCmd->SetParameter(param);
    param = new G4UIparameter("E", 'd', true);
    param->SetDefaultValue("0.0");
    m_fIonCmd->SetParameter(param);

    // daughter ion
    m_dIonCmd = new G4UIcommand("/G4QSim/gun/dIon", this);
    m_dIonCmd->SetGuidance("Set properties of daughter ion.");
    m_dIonCmd->SetGuidance("[usage] /gun/dIon Z A Q E");
    m_dIonCmd->SetGuidance("        Z:(int) AtomicNumber");
    m_dIonCmd->SetGuidance("        A:(int) AtomicMass");
    m_dIonCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
    m_dIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

    param = new G4UIparameter("Z", 'i', false);
    param->SetDefaultValue("1");
    m_dIonCmd->SetParameter(param);
    param = new G4UIparameter("A", 'i', false);
    param->SetDefaultValue("1");
    m_dIonCmd->SetParameter(param);
    param = new G4UIparameter("Q", 'i', true);
    param->SetDefaultValue("0");
    m_dIonCmd->SetParameter(param);
    param = new G4UIparameter("E", 'd', true);
    param->SetDefaultValue("0.0");
    m_dIonCmd->SetParameter(param);

    // beta-neutrino correlation coefficient
    m_beNu_aCmd = new G4UIcmdWithADouble("/G4QSim/gun/beNu_a", this);
    m_beNu_aCmd->SetGuidance("Set beta-neutrino correlation coefficient `a'.");
    m_beNu_aCmd->SetParameterName("beNu_a", true, true);
    m_beNu_aCmd->SetDefaultValue(-1.0/3.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParticleSourceMessenger::~ParticleSourceMessenger()
{
    delete m_listCmd;
    delete m_typeCmd;
    delete m_QValueCmd;
    delete m_EHistCmd;
    delete m_dBDEventsCmd;
    delete m_fIonPDirCmd;
    delete m_fIonPosCmd;
    delete m_fIonEkCmd;
    delete m_fIonCmd;
    delete m_dIonCmd;
    delete m_beNu_aCmd;

    delete m_uiDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticleSourceMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
    if(command == m_listCmd)
        m_particleTable->DumpTable();
    else if(command == m_typeCmd)
        m_particleSource->SetGunType(newValues);
    else if(command == m_QValueCmd)
        m_particleSource->SetQValue(m_QValueCmd->GetNewDoubleValue(newValues));
    else if(command == m_EHistCmd)
        m_particleSource->LoadEHist(newValues);
    else if(command == m_dBDEventsCmd)
        m_particleSource->LoadDBDEvents(newValues);
    else if(command == m_fIonPDirCmd)
        m_particleSource->SetFIonPDir(m_fIonPDirCmd->GetNew3VectorValue(newValues));
    else if(command == m_fIonPosCmd)
        m_particleSource->SetFIonPos(m_fIonPosCmd->GetNew3VectorValue(newValues));
    else if(command == m_fIonEkCmd)
        m_particleSource->SetFIonEk(m_fIonEkCmd->GetNewDoubleValue(newValues));
    else if(command == m_fIonCmd) {
        ParseIonValues(newValues);
        m_particleSource->SetFIon(m_ionDef, m_ionCharge);
    } else if(command == m_dIonCmd) {
        ParseIonValues(newValues);
        m_particleSource->SetDIon(m_ionDef, m_ionCharge);
    } else if(command == m_beNu_aCmd) {
        m_particleSource->SetBeNu_a(m_beNu_aCmd->GetNewDoubleValue(newValues));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition *ParticleSourceMessenger::ParseIonValues(G4String newValues)
{
    G4Tokenizer next(newValues);
    int atomicNumber, atomicMass, ionCharge;
    G4double ionExciteEnergy=0.0;
    // check argument
    atomicNumber = StoI(next());
    atomicMass = StoI(next());
    G4String sQ = next();

    if(sQ.isNull()) {
        ionCharge = atomicNumber;
    } else {
        ionCharge = StoI(sQ);
        sQ = next();
        if(sQ.isNull()) {
            ionExciteEnergy = 0.0;
        } else {
            ionExciteEnergy = StoD(sQ) * keV;
        }
    }
    m_ionDef = G4IonTable::GetIonTable()->GetIon(atomicNumber, atomicMass, ionExciteEnergy);
    if(m_ionDef == 0) {
        G4cout << "Ion with Z=" << atomicNumber;
        G4cout << " A=" << atomicMass << "is not be defined." << G4endl;
    }
    m_ionCharge = ionCharge * eplus;
    return m_ionDef;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
