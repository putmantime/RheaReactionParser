import glob
import os
import urllib.request
import tarfile
import pprint


class ParseRheaReactions(object):
    """
    retrieves reaction data from http://www.rhea-db.org/ ftp site and uses mappings to construct full reaction with
    rhea id, ec number, and chebi ids for all reaction constituents
    """
    def __init__(self):
        self.endpoint = 'ftp://ftp.ebi.ac.uk/pub/databases/rhea/'
        self.reaction_data = self.get_rhea_reaction_rd_data()
        self.r2c = self.rhea2chebi()
        self.ecData = self.get_rhea_ec_tsv_data()
        self.rhea_chebi_ec = self.rhea_chebi2ec()

    def get_rhea_reaction_rd_data(self):
        #         retrieve rhea id with linked chebi ids
        #         returns unarchived collection of rd files for each rhea id
        rd_path = 'ctfiles/rhea-rd.tar.gz'
        reactions_tar = urllib.request.urlretrieve(self.endpoint + rd_path, 'rd.tar.gz')
        reactions = tarfile.open(reactions_tar[0], 'r')
        reactions.extractall()
        reactions.close()

    def get_rhea_ec_tsv_data(self):
        #         retrieve ec rhea mapping
        #         return tempfile of ec 2 rhea mapping
        ec_path = 'tsv/ec-rhea-dir.tsv'
        ec2rhea = urllib.request.urlretrieve(self.endpoint + ec_path)
        return ec2rhea[0]

    def get_chebi_names_list(self):
        #        retrieve tsv of chebi name id pairs
        chebi_list_path = 'tsv/chebiId_name.tsv'
        chebiid2name = urllib.request.urlretrieve(self.endpoint + chebi_list_path)
        cheb_dict = {}
        with open(chebiid2name[0]) as cnl:
            for line in cnl:
                line = line.strip()
                cnl_list = line.split("\t")
                cheb_dict[cnl_list[0]] = cnl_list[-1]
        return cheb_dict

    @staticmethod
    def rhea2chebi():
        #       get chebi ids for reaction from each rhea rd file
        rxn_all = []
        for filename in glob.glob(os.path.join('rd/*.rd')):
            rxn = {}
            rheaid1 = filename.split('.')[0]
            rheaid = rheaid1.split('/')[1]
            with open(filename) as fp:
                rxn['rhea_id'] = rheaid
                reaction = []
                for line in fp:
                    line = line.strip()
                    if 'CHEBI:' in line:
                        reaction.append(line)
                rxn['chebi_id'] = reaction
                rxn_all.append(rxn)
        return rxn_all

    def rhea_chebi2ec(self):
        #      map ec to rxn with rhea ids
        rxn_all = []
        for rxn in self.r2c:
            with open(self.ecData) as ec:
                for line in ec:
                    line = line.strip()
                    eclist = line.split("\t")
                    ec = eclist[0]
                    rhea = eclist[1]
                    if rxn['rhea_id'] == rhea:
                        rxn['ecnumber'] = ec
                        rxn_all.append(rxn)
        return rxn_all




