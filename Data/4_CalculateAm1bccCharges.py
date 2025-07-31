# Calculation of AM1BCC charges for CSD molecules
import argparse
from ForceField.Forcefield import OpenFF_forcefield
import hashlib
import psycopg

def hashStringSha256(input_string):
    sha256 = hashlib.sha256()
    sha256.update(input_string.encode('utf-8'))
    hash_value = sha256.hexdigest()
    return hash_value

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculation of am1bcc charges")
    parser.add_argument("-H","--hostname",type=str,help="The hostname of the postgresql server")
    parser.add_argument("-u","--username",type=str,help="The username for the postgresql server")
    parser.add_argument("-d","--dbname",type=str,help="The database name")
    parser.add_argument("-cd","--cachedir",type=str,help="Path to cache directory")
    args = parser.parse_args()
    cn = psycopg.connect(dbname=args.dbname, user=args.username, host=args.hostname)
    curs = cn.cursor()
    sql = "select molregno, canonical_smiles from csd202403modified.hashes"
    curs.execute(sql)
    dat = curs.fetchall()
    cn.close()
    molregnos = [x[0] for x in dat]
    allSmiles = [x[1] for x in dat]
    cacheids = []
    for smiles in allSmiles:
        hashname = hashStringSha256(smiles)[-3:]
        cacheids.append(hashname)
        OpenFF_forcefield(f"{smiles}_in_v",cache=f"{args.cachedir}/{hashname}.cache")
    cn = psycopg.connect(dbname=args.dbname, user=args.username, host=args.hostname)
    curs = cn.cursor()
    sql = "insert into csd202403modified.caches (molregno, cache) values (%s,%s)"
    curs.executemany(sql, zip(molregnos, cacheids))
    cn.commit()
    curs.close()
    cn.close()