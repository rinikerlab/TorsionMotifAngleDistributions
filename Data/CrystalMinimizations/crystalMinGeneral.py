import openmm
from openmm import app
from openmm.app import ForceField as Forcefield
from openmm.app import Modeller
from openmm.app import CutoffPeriodic
from openmm.unit import nanometer,picosecond,kelvin
from openmm.app import PDBFile
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm.app import CutoffPeriodic
from openff.toolkit.topology import Molecule
import argparse
from rdkit import Chem

debug = False

parser = argparse.ArgumentParser(description="PDB filename")
parser.add_argument("-t","--template",type=str,help="Filename for template for crystal minimization")
parser.add_argument("-f","--filename",type=str,help="Filename for pdb for crystal minimization")
parser.add_argument("-o","--outname",type=str,help="Filename for output pdb")
args = parser.parse_args()

mol = Chem.MolFromMolFile(args.template,removeHs=False)
Chem.AssignStereochemistryFrom3D(mol) 
frags = Chem.GetMolFrags(mol,asMols=True)
for frag in frags:
    if Chem.MolToSmiles(frag) != "[H]O[H]":
        mol = frag
smiles = Chem.MolToSmiles(mol)
if debug:
    print("Molecule SMILES: ")
    print(smiles)

molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
molecule.name = "UNK"

smirnoff = SMIRNOFFTemplateGenerator(molecules=[molecule],forcefield="openff-2.0.0.offxml")
forcefield = Forcefield("amber/tip3p_standard.xml")
forcefield.registerTemplateGenerator(smirnoff.generator)

pdbreader = PDBFile(args.filename)
coordinates = pdbreader.positions
omm_topo = pdbreader.topology
modeller = Modeller(omm_topo, coordinates)
modeller_topo = modeller.getTopology()
modeller_positions = modeller.getPositions()

if debug:
    print("Any unmatched residues: ")
    unmatched_residues = forcefield.getUnmatchedResidues(modeller_topo)
    print(unmatched_residues)

print("Periodic Box Vector: ")
print(modeller_topo.getPeriodicBoxVectors())
print("...done")

# Create the system
print("System creation..")
system = forcefield.createSystem(modeller_topo,nonbondedMethod=CutoffPeriodic,nonbondedCutoff=0.2*nanometer)
print("...done")
print("PBC applied?")
print(system.usesPeriodicBoundaryConditions())
print("...done")

# Create the integrator
integrator = openmm.LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picosecond)

print("Constructing simulation context...")
simulation = app.Simulation(topology=modeller_topo,system=system,integrator=integrator)
simulation.context.setPositions(modeller_positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)
print("...done")        

# Minimize energy
print('Minimizing...')
simulation.minimizeEnergy()

# Save the minimized structure
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(args.outname, 'w'))
print(f'Minimization complete. Output written to {args.outname}')
