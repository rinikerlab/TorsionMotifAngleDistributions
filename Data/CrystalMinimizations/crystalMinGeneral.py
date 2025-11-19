import openmm
from openmm import app
from openmm import MonteCarloBarostat,AndersenThermostat
from openmm import amd
from openmm.app import ForceField as Forcefield
from openmm.app import Modeller
from openmm.app import PME,CutoffPeriodic
from openmm.vec3 import Vec3
from openmm.unit import nanometer,picosecond,femtosecond,kelvin,bar,kilojoule_per_mole,angstrom
from openmm.app import PDBFile
from sys import stdout
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology 
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm.app import CutoffPeriodic, HBonds
from openff import toolkit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
import argparse
from rdkit import Chem
from rdkit.Chem import rdmolops
import numpy as np

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
print('Minimization complete. Output written to {args.outname}')
