#!/usr/bin/env python

#
# Usage: process-otu-data.py INPUT-FILE OUTPUT-FILE
#
# This script processes files containing data on relative abundance of OTUs in
# one or more samples.
#
# Input should be a CSV file with rows representing different OTUs and columns
# representing the sampling site(s). The first column should be labeled "OTU"
# and have the names of each OTU. One of the columns should be labeled "control"
# and have the relative abundance data for the elution buffer. All other columns
# can be labeled freely and should have the relative abundance data for each
# sampling site. Output is another CSV file with columns for each taxonomic rank
# followed by columns for each sampling site. Each row represents the organism
# that an OTU is attributed to and each numerical value represents the relative
# abundance of the OTU converted to a proportion after pruning.
#
# Pruning is done according to relative abundance in the elution buffer. If an
# OTU has less than PRUNING_VALUE times as many reads in a sample as in the
# buffer, its presence is attributed to contamination and the number of hits is
# adjusted to zero. PRUNING_VALUE has been set to 1.5 by default, but may be
# changed.
#
# Open Tree of Life API code is based on https://github.com/brunoasm/TaxReformer
#

# Adjustable pruning value
PRUNING_VALUE = 1.5


import csv, requests, sys, time, warnings
from requests.exceptions import ConnectionError, SSLError
from numpy import nan


# convert Open Tree of Life taxonomy results into dictionary
list2dict = lambda taxlist : {x.split(':')[0] : x.split(':')[1] for x in taxlist}


# query Open Tree of Life's taxonomic resolution services
def otl_tnrs(query, do_approximate = True, wait_time = 600):
    contact_otl = True

    while contact_otl:
        try:
            res = requests.post('https://api.opentreeoflife.org/v3/tnrs/match_names',
                    json = {'names' : [query, query], 'do_approximate_matching' : do_approximate})
        except (SSLError, ConnectionError):
            sys.stderr.write(time.ctime() + ': error connecting to Open Tree of Life, retrying in ' + str(wait_time) + ' seconds')
            time.sleep(wait_time)
            continue

        if res.status_code == 200:
            contact_otl = False
        else:
            sys.stderr.write(time.ctime() + ': error with Open Tree of Life response, retrying in ' + str(wait_time) + ' seconds')
            time.sleep(wait_time)
            continue
    return res


# query Open Tree of Life's taxonomy
def otl_taxon(query, wait_time = 600, ncbi = False):
    contact_otl = True

    while contact_otl:
        try:
            if ncbi:
                res = requests.post('https://api.opentreeoflife.org/v3/taxonomy/taxon_info',
                        json = {'source_id' : 'ncbi:' + str(query), 'include_lineage' : True})
            else:
                res = requests.post('https://api.opentreeoflife.org/v3/taxonomy/taxon_info',
                        json = {"ott_id" : query, "include_lineage" : True})
        except (SSLError, ConnectionError):
            sys.stderr.write(time.ctime() + ': error connecting to Open Tree of Life, retrying in ' + str(wait_time) + ' seconds')
            time.sleep(wait_time)
            continue

        if res.status_code == 200:
            contact_otl = False
        elif res.status_code == 400:
            contact_otl = False
            sys.stderr.write(res.json()['message'])
            sys.stderr.write('skipping')
            return None
        else:
            sys.stderr.write(time.ctime() + ': error contacting Open Tree of Life, retrying in ' + str(wait_time) + ' seconds')
            time.sleep(wait_time)
    return res


# get taxonomy up to order from a genus name
def taxonomy_OTT(ott_id = None):
    res = otl_taxon(ott_id, wait_time = 3600)

    # taxonomy dictionary
    out_dict = {('tax_' + higher['rank']) : higher['name'] for higher in res.json()['lineage']}
    out_dict['tax_higher_source'] = 'OTT'
    out_dict['rank'] = res.json()['rank']

    # remove unnecessary ranks
    for rank in ['tax_no rank']:
        out_dict.pop(rank,0)

    # add OTT id and accepted name
    out_dict['tax_ott_id'] = ott_id
    out_dict['tax_ott_accepted_name'] = res.json()['name']

    # add NCBI id
    try:
        out_dict.update({'tax_ncbi_id' : list2dict(res.json()['tax_sources'])['ncbi']})
    except KeyError:
        pass

    # add genus information for species or subspecies queries
    if res.json()['rank'] in ['species', 'subspecies']:
        try:
            genus_tax = [tax for tax in res.json()['lineage'] if tax['rank'] == 'genus'][0]
        except IndexError:
            out_dict['tax_cg_ott_id'] = nan
            if out_dict['rank'] in ['species', 'subspecies']:
                out_dict['cg'] = res.json()['unique_name'].split()[0]
        else:
            out_dict['tax_cg_ott_id'] = genus_tax['ott_id']
            out_dict['cg'] = out_dict['tax_genus']
            del out_dict['tax_genus']
            try:
                out_dict.update({'tax_cg_ncbi_id' : list2dict(genus_tax['tax_sources'])['ncbi']})
            except KeyError:
                pass

    # update species ids for subspecies queries
    try:
        species_tax = [tax for tax in res.json()['lineage'] if tax['rank'] == 'species'][0]
        out_dict['tax_cs_ott_id'] = species_tax['ott_id']
    except:
        pass

    try:
        del out_dict['tax_species']
    except KeyError:
        pass

    if out_dict['rank'] not in ['species', 'subspecies', 'genus', 'subgenus']:
        out_dict['cg'] = ''

    return out_dict


# get taxonomy information from a query
def otl_checkname(query):
    outdict = None

    res = otl_tnrs(query, do_approximate = False)

    if res.json()['results']:
        if not res.json()['results'][0]['matches']:
            return { }

        result = res.json()['results'][0]['matches'][0]
        outdict = {'current_name' : result['taxon']['name'], 'id' : result['taxon']['ott_id'], 'name_source' : 'OTT'}
        outdict['higher_taxonomy'] = taxonomy_OTT(result['taxon']['ott_id'])

        if result['taxon']['rank'] == 'species':
            outdict['level'] = 'species'
            try:
                outdict['ncbi_id'] = list2dict(result['taxon']['tax_sources'])['ncbi']
            except:
                pass
        elif result['taxon']['rank'] == 'genus':
            outdict['level'] = 'genus'
        else:
            warnings.warn(query + ': found in ott, but not as genus or species')
            return None

    return outdict


# read data from csv file
def read_data(file):
    try:
        with open(file, newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            header = [item for item in next(csvreader)]
            header_lower = [item.lower() for item in header]
            csvrows = list(csvreader)
    except FileNotFoundError:
        sys.stderr.write('file \'' + file + '\' not found\n')
        sys.exit(1)

    csvcols = list(zip(*csvrows))

    # first column should be 'OTU'
    if header_lower[0] != 'otu':
        sys.stderr.write('no OTU column found\n')
        sys.exit(1)

    species = [item.lower() for item in csvcols[0]]

    # find 'control' column
    try:
        control_col = header_lower.index('control')
    except ValueError:
        sys.stderr.write('no control column found\n')
        sys.exit(1)

    controls = [int(item) for item in csvcols[control_col]]

    # find all experimental columns
    data = [[int(row[i]) for i in range(1, len(row)) if i != control_col] for row in csvrows]

    return header, species, controls, data


# process data
def process_data(species, controls, data):
    species_processed = []
    data_processed = []

    # prune counts below cutoff
    for i in range(len(data)):
        data_processed += [[]]
        for j in range(len(data[i])):
            if controls[i] and data[i][j] / controls[i] <= PRUNING_VALUE:
                data_processed[i] += [0]
            else:
                data_processed[i] += [data[i][j]]

    # get sums of pruned counts
    sums = [sum(i) for i in list(zip(*data_processed))]

    # get proportions
    for i in range(len(data_processed)):
        for j in range(len(data_processed[i])):
            data_processed[i][j] = data_processed[i][j] / sums[j]

    # get taxonomy data
    for sp in species:
        name = [sp]
        species_processed += [[]]

        taxa = otl_checkname(sp)

        while not taxa and ' ' in name[0]:
            name = name[0].rsplit(' ', 1) + name[1:]
            taxa = otl_checkname(name[0])

        if taxa:
            for key in ['tax_domain', 'tax_kingdom', 'tax_phylum', 'tax_class', 'tax_order', 'tax_family']:
                try:
                    species_processed[-1] += [taxa['higher_taxonomy'][key]]
                except KeyError:
                    species_processed[-1] += ['-']

            species_processed[-1] += [taxa['current_name']]

            if len(name) > 1:
                species_processed[-1][-1] += ' ' + ' '.join(name[1:])

        else:
            species_processed[-1] += ['-', '-', '-', '-', '-', '-', sp]

    return species_processed, data_processed


# write processed data to csv file
def write_data(file, header, species, data):
    with open(file, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus, species'] + header[1:-1])

        combined = [species[i] + data[i] for i in range(len(species))]
        combined.sort()

        csvwriter.writerows(combined)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.stderr.write('no input file given\n')
        sys.exit(1)

    if len(sys.argv) == 2:
        sys.stderr.write('no output file give\n')
        sys.exit(1)

    header, species, controls, data = read_data(sys.argv[1])

    print('processing...')

    species_processed, data_processed = process_data(species, controls, data)
    write_data(sys.argv[2], header, species_processed, data_processed)

    print('done.')
