import psycopg
import numpy as np
from rdkit import Chem

CONFIG_MULTISOLVENT = None

def SetLoginInfoAndConfigMultisolvent():
    global CONFIG_MULTISOLVENT
    CONFIG_MULTISOLVENT = {}
    CONFIG_MULTISOLVENT['dbtype'] ='postgresql'
    CONFIG_MULTISOLVENT['host'] = 'ADD HOSTNAME'
    CONFIG_MULTISOLVENT['username'] = 'ADD USERNAME'
    CONFIG_MULTISOLVENT['dbname'] = 'ADD DBNAME'

SetLoginInfoAndConfigMultisolvent()

def GetSmarts(entrynr):
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select distinct alteredsmarts from extractions.et4hierarchy where entrynr=%s;"
    curs.execute(sql, (entrynr,))
    rows = curs.fetchall()
    ans = [row for row in rows]
    cn.rollback()
    cn.close()
    return ans[0][0] if ans else None

def GetHierarchyclass(entrynr):
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select distinct hierarchyclass from extractions.et4hierarchy where entrynr=%s;"
    curs.execute(sql, (entrynr,))
    rows = curs.fetchall()
    ans = [row for row in rows]
    cn.rollback()
    cn.close()
    return ans[0][0] if ans else None

def RetrieveRawDataCsdModified(entrynr):
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select alteredsmarts from extractions.et4hierarchy where entrynr=%s;"
    curs.execute(sql, (entrynr,))
    rows = curs.fetchall()
    cn.rollback()
    sql = """select csd.torsionvalue, csd.molregno 
    from (select molregno, dihedral from extractions.et4csd202403modified where entrynr=%(entrynr)s and hierarchy is true) 
    as subquery 
    join csd202403modified.profiles as csd
    on subquery.molregno=csd.molregno and subquery.dihedral=csd.dihedral;"""
    curs.execute(sql, {'entrynr':entrynr})
    rows = curs.fetchall()
    cn.rollback()
    data = []
    for row in rows:
        tmpData = np.array(row[0],dtype=np.float64)
        data.append(tmpData)
    data = np.array(data)
    data = np.deg2rad(data)
    cn.close()
    return data

def RetrieveRawDataDash(entrynr, threshold):
    environment = 'tip3p'
    envMetadata = environment + 'metadata'
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select distinct alteredsmarts from extractions.et4hierarchy where entrynr=%s;"
    curs.execute(sql, (entrynr,))
    rows = curs.fetchall()
    cn.rollback()
    sql = f"""SELECT {environment}.torsionvalues, {envMetadata}.energiesnormalized
        FROM (
            SELECT molregno, dihedral 
            FROM extractions.et4 
            WHERE entrynr=%(entrynr)s AND hierarchy IS TRUE
        ) AS subquery
        JOIN gnn.{environment} AS {environment} 
        ON subquery.molregno = {environment}.molregno AND subquery.dihedral = {environment}.dihedral
        JOIN gnn.{envMetadata} AS {envMetadata}
        ON subquery.molregno = {envMetadata}.molregno;"""
    curs.execute(sql, {'entrynr':entrynr})
    rows = curs.fetchall()
    cn.rollback()
    data = []
    es = []
    for row in rows:
        tmpData = np.array(row[0],dtype=np.float64)
        tmpEs = np.array(row[1],dtype=np.float64)
        data.append(tmpData)
        es.append(tmpEs)
    if not data:
        return np.array([])
    data = np.concatenate(data)
    es = np.concatenate(es)
    datEs = np.column_stack((data,es))
    # print(datEs.shape)
    data = datEs[datEs[:,1]<threshold][:,0]
    data = np.deg2rad(data)
    cn.close()
    return data

def RetrieveRawDataGNN(entrynr, threshold, environment):
    allowedEnvs = ['vac', 'tip3p', 'hexane']
    if environment not in allowedEnvs:
        raise ValueError(f"Environment must be one of {allowedEnvs}")
    envMetadata = environment + 'metadata'
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select distinct alteredsmarts from extractions.et4hierarchy where entrynr=%s;"
    curs.execute(sql, (entrynr,))
    rows = curs.fetchall()
    cn.rollback()
    sql = f"""SELECT {environment}.torsionvalues, {envMetadata}.energiesnormalized
        FROM (
            SELECT molregno, dihedral 
            FROM extractions.et4csd202403modified 
            WHERE entrynr=%(entrynr)s AND hierarchy IS TRUE
        ) AS subquery
        JOIN gnn.csd{environment} AS {environment} 
        ON subquery.molregno = {environment}.molregno AND subquery.dihedral = {environment}.dihedral
        JOIN gnn.csd{envMetadata} AS {envMetadata}
        ON subquery.molregno = {envMetadata}.molregno;"""
    curs.execute(sql, {'entrynr':entrynr})
    rows = curs.fetchall()
    cn.rollback()
    data = []
    es = []
    for row in rows:
        tmpData = np.array(row[0],dtype=np.float64)
        tmpEs = np.array(row[1],dtype=np.float64)
        data.append(tmpData)
        es.append(tmpEs)
    if not data:
        return np.array([])
    data = np.concatenate(data)
    es = np.concatenate(es)
    datEs = np.column_stack((data,es))
    # print(datEs.shape)
    data = datEs[datEs[:,1]<threshold][:,0]
    data = np.deg2rad(data)
    cn.close()
    return data

def RetrieveEnsembleDataVac(molregno):
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select molregno,molblock from csd202403modified.molblocks where molregno=%s"
    curs.execute(sql, (molregno,))
    rows = curs.fetchall()
    cn.rollback()
    mol = Chem.MolFromMolBlock(rows[0][1],removeHs=False,sanitize=False)
    conf = Chem.Conformer(mol.GetConformer(0))
    mol.RemoveAllConformers()
    sql = "select coords from gnn.csdvacmetadata where molregno=%s"
    curs.execute(sql, (molregno,))
    rows = curs.fetchall()
    cn.rollback()
    coords = np.array(rows,dtype=np.float64)[0][0]
    for i in range(np.shape(coords)[0]):
        for j in range(np.shape(coords)[1]):
            conf.SetAtomPosition(j,coords[i][j]*10)
        mol.AddConformer(conf, assignId=True)
    return mol

def RetrieveEnsembleDataTip3p(molregno):
    cn = psycopg.connect(dbname=CONFIG_MULTISOLVENT['dbname'], user=CONFIG_MULTISOLVENT['username'], host=CONFIG_MULTISOLVENT['host'])
    curs = cn.cursor()
    sql = "select molregno,molblock from csd202403modified.molblocks where molregno=%s"
    curs.execute(sql, (molregno,))
    rows = curs.fetchall()
    cn.rollback()
    mol = Chem.MolFromMolBlock(rows[0][1],removeHs=False,sanitize=False)
    conf = Chem.Conformer(mol.GetConformer(0))
    mol.RemoveAllConformers()
    sql = "select coords from gnn.csdtip3pmetadata where molregno=%s"
    curs.execute(sql, (molregno,))
    rows = curs.fetchall()
    cn.rollback()
    coords = np.array(rows,dtype=np.float64)[0][0]
    for i in range(np.shape(coords)[0]):
        for j in range(np.shape(coords)[1]):
            conf.SetAtomPosition(j,coords[i][j]*10)
        mol.AddConformer(conf, assignId=True)
    return mol