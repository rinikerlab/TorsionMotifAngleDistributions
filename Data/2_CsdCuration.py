# Script to populate the modified CSD lwreg instance
# also calculates the cache names, used for the am1bcc charge calculations

import pickle
from rdkit import Chem
import lwreg
import argparse
import hashlib
import psycopg

def hashStringSha256(input_string):
    sha256 = hashlib.sha256()
    sha256.update(input_string.encode('utf-8'))
    hash_value = sha256.hexdigest()
    return hash_value

def _needsHs(mol):
    for atom in mol.GetAtoms():
         nHNbrs = 0
         for nbri in atom.GetNeighbors():
              if nbri.GetAtomicNum() == 1:
                  nHNbrs+=1
         if atom.GetTotalNumHs(False) > nHNbrs:
            return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Curation on the initial CSD instance")
    parser.add_argument("-H","--hostname",type=str,help="The hostname of the postgresql server")
    parser.add_argument("-u","--username",type=str,help="The username for the postgresql server")
    parser.add_argument("-d","--dbname",type=str,help="The name of the database on postgresql server")
    parser.add_argument("-f","--file",type=str,help="Pkl file with bulk register tuples or failure reasons")
    args = parser.parse_args()

    with open(args.file,"rb") as f:
        allTuples = pickle.load(f)
    # filter out entries that are not tuples of integers, this filters out the RegistrationFailureReasons 
    allTuples = [x for x in allTuples if isinstance(x, tuple) and len(x) == 2 and all(isinstance(i, int) for i in x)]

    badAtoms = Chem.MolFromSmarts('[!#1;!#6;!#7;!#8;!#9;!#15;!#16;!#17;!#35;!#53]')
    atleastOneCC = Chem.MolFromSmarts('[#6]~[#6]')
    mols = []
    filteredTuples = []
    config = lwreg.configure_from_database(dbname=args.dbname,dbtype="postgresql",user=args.username,host=args.hostname,lwregSchema="csd202403")
    for molregnoConfid in allTuples:
        molBlock = lwreg.retrieve(config,id=tuple(molregnoConfid))
        mol = Chem.MolFromMolBlock(molBlock[tuple(molregnoConfid)][0],removeHs=False,sanitize=False)
        matches = mol.GetSubstructMatches(badAtoms)
        if not matches:
            Chem.AssignStereochemistryFrom3D(mol)
            molCopy = Chem.GetMolFrags(mol,asMols=True,sanitizeFrags=False)
            for mol in molCopy:
                if _needsHs(mol):
                    mol = Chem.AddHs(mol,addCoords=True)
                matches2 = mol.GetSubstructMatches(atleastOneCC)
                if matches2:
                    mols.append(tuple([molregnoConfid,mol]))
                else:
                    filteredTuples.append(molregnoConfid)
        else:
            filteredTuples.append(molregnoConfid)
    allMols = [x[1] for x in mols]
    modifiedConfig = lwreg.configure_from_database(dbname=args.dbname,dbtype="postgresql",user=args.username,host=args.hostname,lwregSchema="csd202403modified")
    ans = lwreg.bulk_register(modifiedConfig,allMols)

    # populate the csd202403modifiedrdk.mols table for the rdkit cartridge
    cn = psycopg.connect(dbname=args.dbname, user=args.username, host=args.hostname)
    curs = cn.cursor()
    curs.execute('drop table if exists csd202403modifiedrdk.mols cascade')
    curs.execute('select molregno, mol_from_ctab(molblock::cstring,false) m into csd202403modifiedrdk.mols from csd202403modified.molblocks')
    curs.execute('create index molidx on csd202403modifiedrdk.mols using gist(m)')
    cn.commit()
    curs.close()
    cn.close()
