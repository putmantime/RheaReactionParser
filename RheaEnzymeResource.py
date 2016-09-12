import Bio.ExPASy.Enzyme as bee
import re
import pprint
import requests
import urllib.request
import sys
import pandas as pd
import tempfile
import glob
import os
import urllib.request
import tarfile
from pymongo import MongoClient
from datetime import datetime


client = MongoClient()
reactions_db = client['reactions']
expasy_coll = reactions_db.expasy
rhea_coll = reactions_db.rhea


class ParseRheaReactions(object):
    """
    retrieves reaction data from http://www.rhea-db.org/ ftp site and uses mappings to construct full reaction with
    rhea id, ec number, and chebi ids for all reaction constituents
    """
    def __init__(self):
        self.endpoint = 'ftp://ftp.ebi.ac.uk/pub/databases/rhea/'
        self.chebi_names = self.get_chebi_names_list()
        self.rhea2ec = self.get_rhea_ec_tsv_data()
        self.r2c = self.rhea2chebi()

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
        ec_dict = {}
        with open(ec2rhea[0]) as ec:
            ec_dict = {}
            for line in ec:
                line = line.strip()
                ec_list = line.split("\t")
                ec_dict[ec_list[1]] = ec_list[0]
            return ec_dict

    def get_chebi_names_list(self):
        #        retrieve tsv of chebi name id pairs
        chebi_list_path = 'tsv/chebiId_name.tsv'
        chebiid2name = urllib.request.urlretrieve(self.endpoint + chebi_list_path)
        chebi_dict ={}
        with open(chebiid2name[0]) as cnl:
            for line in cnl:
                line = line.strip()
                cnl_list = line.split("\t")
                chebi_dict[cnl_list[0]] = cnl_list[-1]

        return chebi_dict

    def rhea2chebi(self):
        #       get chebi ids for reaction from each rhea rd file
        rxn_all = []
        count = 0
        for filename in glob.glob(os.path.join('rd/*.rd')):

            rxn = {}
            rxn['timestamp'] = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.now())
            rheaid1 = filename.split('.')[0]
            rheaid = rheaid1.split('/')[1]
            with open(filename) as fp:
                rxn['_id'] = rheaid
                reaction = {}
                for line in fp:
                    line = line.strip()
                    # extract chebi id and link name from self.chebi_names
                    if 'CHEBI:' in line:
                        try:
                            reaction[line] = self.chebi_names[line]
                        except Exception as e:
                            reaction[line] = 'none'
                rxn['chebi_id'] = reaction
                try:
                    rxn['ecnumber'] = self.rhea2ec[rxn['_id']]
                except Exception as e:
                    next
                rxn_all.append(rxn)
                # find and replace entire doc if exists else insert new if does not exist
                rhea_coll.replace_one({'_id': rxn['_id']}, replacement=rxn, upsert=True)

        return rxn_all


def link_compound2chebi(compound, apik):
    def replace_strings(compound_name):
        """
        takes list of compound names and replaces with names formatted for NCBO annotator
        """
        names = {'H(2)O': 'water',
                 'CO(2)': 'carbon dioxide',
                 'An alcohol': 'alcohol',
                 'A secondary alcohol': 'secondary alcohol',
                 'O(2)': 'dioxygen',
                 'H(+)': 'hydron',
                 'a ketone': 'ketone',
                }

        if compound_name in names:
            return names[compound_name]

        else:
            return compound_name

    compound = replace_strings(compound)
    """
    used NCBO Annotator from BioPortal to return ChEBI IDS
    for substrates and products of reactions from Expasy enzyme

    """
    url = 'http://data.bioontology.org/annotator'
    params = dict(apikey=apik, text=compound, ontologies='CHEBI', longest_only='true',
                  include='properties', exlude_numbers='false', exclude_synonyms='false', mappins='all')
    tm_results = requests.get(url=url, params=params).json()

    for i in tm_results:
        prefLabel = i['annotatedClass']['properties']['http://data.bioontology.org/metadata/def/prefLabel']
        text_input = i['annotations'][0]['text']
        if prefLabel[0].lower() == text_input.lower():
            final_chebi = i['annotatedClass']['@id'].split("/")[-1]
            return final_chebi.split("_")[-1]


class AnnotateExpasy(object):
    def __init__(self, apikey):
        self.apikey = apikey
        self.reactions = self.get_parse_expasy_reactions()

    def get_parse_expasy_reactions(self):
        reaction_list = []
        count = 0
        url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"
        print("Retrieving enzyme records from Expasy Enzyme")
        enzyme = urllib.request.urlretrieve(url)
        enzyme_p = bee.parse(open(enzyme[0], 'r'))
        # use BioPython to parse enzyme records
        for record in enzyme_p:
            enz_rec = {
                'reaction(s)': {},
                '_id': record['ID'],
                'description': record['DE'].strip("."),
                'timestamp': '{:%Y-%m-%d %H:%M:%S}'.format(datetime.now())
            }
            count += 1

            # parse reactions to remove leading numbers and populate dictionary
            reactions = [x for x in record['CA'].split(".") if x]
            reactions = [x[3:].strip() if len(reactions) > 1 else x for x in reactions]
            rcount = 0
            # parse the constituents out of each side of reaction
            for rxn in reactions:
                print("reaction_" + str(count))
                rcount += 1
                constituents = rxn.split('=')
                if len(constituents) == 2:
                    left = constituents[0]
                    right = constituents[1]
                    # split each side of reaction on '+' not '(+)'
                    r = re.compile(r'(?:[^\+(]|\([^)]*\))+')
                    left = [x.strip() for x in r.findall(left)]
                    left = {x: link_compound2chebi(x, self.apikey) for x in left}
                    right = [x.strip() for x in r.findall(right)]
                    right = {x: link_compound2chebi(x, self.apikey) for x in right}
                    enz_rec['reaction(s)']['rxn_' + str(rcount)] = {'reaction': rxn,
                                                                    'left': left,
                                                                    'right': right
                                                                     }
                    print(enz_rec['reaction(s)']['rxn_' + str(rcount)])
            print('rxn_' + str(count))
            expasy_coll.replace_one({'_id': enz_rec['_id']}, replacement=enz_rec, upsert=True)
            reaction_list.append(enz_rec)

        return reaction_list
#
api_key = '0a99c359-d2a2-483a-8dca-148c3bb4e8c1'
AnnotateExpasy(apikey=api_key)