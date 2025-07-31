from rdkit import Chem
import psycopg
import argparse

def ExtractDihedralsMatchingMotifs(dbname, hostname, username):
    cn = psycopg.connect(dbname=dbname,user=username,host=hostname)
    curs = cn.cursor()
    sql = "select entrynr, alteredsmarts from extractions.et4hierarchy where ignored is not true and entrynr<422 order by entrynr;"
    curs.execute(sql)
    results = curs.fetchall()
    entrynrs = [r[0] for r in results]
    alteredSmarts = [r[1] for r in results]

    for entrynr, smi in zip(entrynrs,alteredSmarts):
        toRegister = []
        curs.execute("select * from csd202403modifiedrdk.mols where m @> %s::qmol",(smi,))
        d = curs.fetchall()
        sma = Chem.MolFromSmarts(smi)
        atomMap = {}
        for atom in sma.GetAtoms():
            atomMap[atom.GetAtomMapNum()] = atom.GetIdx()
        for m in d:
            sql = "select molblock from csd202403modified.molblocks where molregno=%s"
            curs.execute(sql,(m[0],))
            molblock = curs.fetchall()
            mol = Chem.MolFromMolBlock(molblock[0][0],removeHs=False)
            hits = mol.GetSubstructMatches(sma)
            for hit in hits:
                tmp = []
                tmp2 = []
                tmp = [hit[atomMap[1]],hit[atomMap[2]],hit[atomMap[3]],hit[atomMap[4]]]
                if hit[atomMap[2]] < hit[atomMap[3]]:
                    tmp2 = [hit[atomMap[2]],hit[atomMap[3]]]
                else:
                    tmp2 = [hit[atomMap[3]],hit[atomMap[2]]]
                toRegister.append([int(entrynr),m[0],tmp,tmp2])
        sql = "insert into extractions.et4csd202403modified (entrynr,molregno,dihedral,hierarchy,bond) values (%s,%s,%s,false,%s)"
        curs.executemany(sql,toRegister)
        cn.commit()

    # matches with hierarchy switched on 
    sql = "SELECT DISTINCT ON (bond) entrynr, molregno, dihedral, bond FROM extractions.et4csd202403modified ORDER BY bond, entrynr ASC"
    curs.execute(sql)
    results = curs.fetchall()
    sql = "insert into extractions.et4csd202403modified (entrynr,molregno,dihedral,hierarchy,bond) values (%s,%s,%s,true,%s)"
    curs.executemany(sql,results)
    cn.commit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract dihedrals matching motifs from CSD")
    parser.add_argument("-H","--hostname",type=str,help="The hostname of the postgresql server")
    parser.add_argument("-d","--dbname",type=str,help="The database name")
    parser.add_argument("-u","--username",type=str,help="The username for the postgresql server")
    args = parser.parse_args()
    ExtractDihedralsMatchingMotifs(args.dbname, args.hostname, args.username)