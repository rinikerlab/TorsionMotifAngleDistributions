from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
import argparse
from GNNImplicitSolvent import minimize_mol
import psycopg
import numpy as np
from lwreg import utils
import lwreg
import time
import hashlib

CONFIG = None

def SetLoginInfoAndConfig(hostname, username, dbname):
    global CONFIG
    CONFIG = {}
    CONFIG['dbtype'] ='postgresql'
    CONFIG['host'] = hostname
    CONFIG['username'] = username
    CONFIG['dbname'] = dbname
    CONFIG['schema'] = 'csd202403modified'

# function to handle the connection and the context handler 
# function is just gonna return the connection
def GetConnection(retries=10, delay=10):
    attempt = 0
    while attempt < retries:
        try:
            cn = psycopg.connect(dbname=CONFIG['dbname'], user=CONFIG['username'], host=CONFIG['host'])
            return cn
        except psycopg.OperationalError as e:
            print(f"Connection failed: {e}")
            attempt += 1
            if attempt < retries:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print("All retry attempts failed.")
                raise
    return None

def GNNConformerGenerationAndTorsionExtraction(mol,molregno,cacheid,solvent,cachedir):
    if solvent not in ["vac","tip3p","dmso","methanol","hexane","chloroform"]:
        print("Solvent not in allowed entries.")
        raise
    paulsDict = {
            "vac": "vac",
            "tip3p": "tip3p",
            "dmso": "DMSO",
            "methanol": "Methanol",
            "hexane": "Hexan",
            "chloroform": "Chloroform"
            }
    schemaMetadata = solvent+"metadata"

    cn = GetConnection()
    curs = cn.cursor()
    sql = f"select exists(select 1 from gnn.csd{solvent} where molregno = %s)"
    curs.execute(sql,[molregno])
    exists = curs.fetchone()[0]
    if exists:
        return
    # skipping the general torsions atm
    sql = "select distinct dihedral from extractions.et4csd202403modified where molregno=%s and entrynr<422;"
    curs.execute(sql,(molregno,))
    dihedrals = curs.fetchall()
    cn.close()
    dihedrals = [tuple(d[0]) for d in dihedrals]

    ps = AllChem.ETKDGv3()
    ps.useBasicKnowledge = True
    ps.useExpTorsionAnglePrefs = False
    ps.useSmallRingTorsions = False
    ps.useMacrocycleTorsions = False
    ps.randomSeed = 42
    ps.pruneRmsThresh = -1
    allCoords = []
    AllChem.EmbedMultipleConfs(mol,numConfs=100,params=ps)
    output = minimize_mol(mol,paulsDict[solvent],"./Models/ProductionRun_seed_1612_49_ckpt.pt",strides=10,cache=f"{cachedir}/{cacheid}.cache",return_traj=True)
    mol,traj,ene = output
    eneNorm = ene - np.min(ene)
    cids = [x.GetId() for x in mol.GetConformers()]
    allCoords = []
    allTorsionDist = {d: [] for d in dihedrals}
    for i in cids:
        conf = mol.GetConformer(cids[i])
        for dihedral in dihedrals:
            tmp = rdMolTransforms.GetDihedralDeg(conf, dihedral[0], dihedral[1], dihedral[2], dihedral[3])
            if tmp < 0:
                tmp += 360
            allTorsionDist[dihedral].append(tmp)
        allCoords.append(traj[i].xyz[0])

    dihedralsForDb = [np.array(d).tolist() for d in dihedrals]
    torsionsForDb = []
    for key in allTorsionDist:
        torsionsForDb.append(np.array(allTorsionDist[key]).tolist())
    datToInsert = [(molregno,d,tv) for d,tv in zip(dihedralsForDb,torsionsForDb)]

    cn = GetConnection()
    curs = cn.cursor()
    coordsForDb = np.array(allCoords).tolist()
    eneForDb = ene.tolist()
    eneNormForDb = eneNorm.tolist()
    sql = f"""insert into gnn.csd{solvent} (molregno, dihedral, torsionvalues) 
    values (%s, %s, %s)
    on conflict (molregno, dihedral)
    do update set torsionvalues = excluded.torsionvalues;
    """
    curs.executemany(sql, datToInsert)
    cn.commit()
    sql =  f"""insert into gnn.csd{schemaMetadata} (molregno, coords, energies, energiesnormalized)
    values (%s, %s, %s, %s)
    on conflict (molregno)
    do update set coords = excluded.coords, energies = excluded.energies, energiesnormalized = excluded.energiesnormalized"""
    curs.execute(sql,(molregno,coordsForDb,eneForDb,eneNormForDb))
    cn.commit()
    cn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Conformer Generation with KDG and then minimization")
    parser.add_argument("-H","--host",type=str,help="Host of the database")
    parser.add_argument("-u","--username",type=str,help="Username for the database")
    parser.add_argument("-d","--dbname",type=str,help="Database name")
    parser.add_argument("-f","--file",type=str,help="File with molregnos, cacheid")
    parser.add_argument("-s","--solvent",type=str,help="solvent to minimize in")
    parser.add_argument("-cd","--cachedir",type=str,help="Path to cache directory")
    args = parser.parse_args()
    SetLoginInfoAndConfig(args.host, args.username, args.dbname)
    molregnosAndCacheid = np.loadtxt(args.file,delimiter=",",dtype=str)
    molregnosAndCacheidDict = {}
    for x in molregnosAndCacheid:
        molregnosAndCacheidDict[int(x[0])]=x[1]
    molregnos = molregnosAndCacheid[:,0]
    molregnos = [int(x) for x in molregnos]
    cn = GetConnection()
    curs = cn.cursor()
    placeholders = ', '.join(['%s'] * len(molregnos))
    sql = f"select molregno,molblock from csd202403modified.molblocks where molregno in ({placeholders})"
    curs.execute(sql,molregnos)
    dat = curs.fetchall()
    datDict = dict(dat)
    cn.close()
    for row in molregnosAndCacheid:
        molregno = int(row[0])
        cacheid = str(row[1])
        print(molregno)
        mol = Chem.MolFromMolBlock(datDict[molregno],removeHs=False)
        Chem.AssignStereochemistryFrom3D(mol)
        try:
            GNNConformerGenerationAndTorsionExtraction(mol,molregno,cacheid,args.solvent,args.cachedir)
        except:
            cn = GetConnection()
            curs = cn.cursor()
            schemaMetadata = args.solvent+"metadata"
            sql =  f"""insert into gnn.csd{schemaMetadata} (molregno, error)
            values (%s, %s)
            on conflict (molregno)
            do update set error = excluded.error"""
            curs.execute(sql,(molregno,True))
            cn.commit()
            cn.close()
            continue
