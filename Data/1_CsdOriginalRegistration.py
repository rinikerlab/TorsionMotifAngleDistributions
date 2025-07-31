import lwreg
from lwreg import utils
import argparse
import pickle

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Registration of the initial CSD instance")
    parser.add_argument("-H","--hostname",type=str,help="The hostname of the postgresql server")
    parser.add_argument("-u","--username",type=str,help="The username for the postgresql server")
    parser.add_argument("-d","--dbname",type=str,help="The name of the database to register the CSD data")
    parser.add_argument("-f","--file",type=str,help="Filename of the sdf with the CSD data")
    args = parser.parse_args()
    orgConfig = lwreg.configure_from_database(dbname=args.dbname,dbtype="postgresql",user=args.username,host=args.hostname,lwregSchema="csd202403")
    ans = utils.bulk_register(config=orgConfig, sdfile=args.file)
    # dump as to a logfile in pkl format
    with open("csd202403registration.pkl","wb") as f:
        pickle.dump(ans,f)