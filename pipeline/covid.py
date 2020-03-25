"""covid.py

Utilities to process COVID data.

== Converting into LIF

$ python3 covid.py --convert METADATA_FILE DATA_DIR OUT_DIR N?

Converts Covid JSON files in DATA_DIR while using the COVID metadata in
METADATA_FILE. The DATA_DIR could be the comm_use_subset or any of the other
subsets of the data and METADATA_FILE is whatever COVID metadata you downloaded,
in this case it was all_sources_metadata_2020-03-13.csv. LIF files are written
to OUT_DIR. The optional last argument can be used to restrict processing to
some number of files, default is to process all of them.

== Creating a relations file

$ python3 covid.py --create-relations METADATA_FILE PROCESSING_RESULTS OUT_FILE

Creates a file with a Python object containing all relations structured around
reified relations (relation-object pairs, see the docstring in reify_relations()
for more details). Requires the COVID metadata file and the Harvard processing
results.


== Importing Harvard processing results

$ python3 covid.py --import METADATA_FILE PROCESSING_RESULTS LIF_DIR HAR_DIR N?

For each file in LIF_DIR a file will be created in HAR_DIR which has the reified
relations from the Harvard data as metadata (because there were no offsets).


"""

import os
import sys
import csv
import json
import random
import textwrap

from io import StringIO
from collections import Counter

from lif import LIF, View, Text, Annotation


RELTYPES = {'Activation': ('activates', 'activator'),
            'Inhibition': ('inhibits','inhibitor'),  
            'IncreaseAmount': ('increases', 'increaser'),
            'DecreaseAmount': ('decreases', 'decreaser')}


def convert_into_lif(metadata_file, data_dir, out_dir, n=99999):
    """Load Covid metadata and convert Covid JSON files fron data_dir into LIF files
    and save them in out_dir."""
    print('Loading metadata...')
    metadata = Metadata(metadata_file)
    for fname in os.listdir(data_dir)[:n]:
        infile = os.path.join(data_dir, fname)
        outfile = os.path.join(out_dir, fname)
        Converter(infile, outfile, metadata).convert()


def create_relations_file(metadata_file, results_file, out_file):
    """Create a file with reified relations and write it to out_file. This requires
    the COVID metadata and the Harvard processing results."""
    print('Loading metadata...')
    metadata = Metadata(metadata_file)
    print('Loading Harvard processing results...')
    results = HarvardResults(results_file)
    rels = results.collect_relations(metadata)
    idx = index_by_fname(rels)
    reified_rels = reify_relations(rels)
    for rt in reified_rels:
        print(rt, len(reified_rels[rt]))
    json.dump(reified_rels, open(out_file, 'w'), indent=True)


def translate_reltype_into_action(reltype):
    return RELTYPES.get(reltype, (None, None))[0]

def translate_reltype_into_role(reltype):
    return RELTYPES.get(reltype, (None, None))[1]


class Metadata(object):

    """Reads the CSV file with the COVID metadata and stores each line as a list of
    fields in the data instance variable."""

    # Some counts (first column for 3/13 download, second for 3/20 download):
    #
    # doi        26,358  40,751
    # pmcid      27,338  23,320
    # pubmed_id  16,731  22,944

    def __init__(self, csv_file):
        self.fh = open(csv_file)
        self.reader = csv.reader(self.fh)
        self.data = []
        for fields in self.reader:
            self.data.append(fields)
        self.sha2pmid = {}
        self.sha2year = {}
        self.pmid2sha = {}
        for fields in self.data:
            sha = fields[0]
            pmcid = fields[4]
            pmid = fields[5]
            year = fields[8]
            if len(year) > 4 and year[4] == ' ' and year[:4].isdigit():
                year = year[:4]
            self.sha2pmid[sha] = pmid
            self.sha2year[sha] = year
            self.pmid2sha[pmid] = sha

    def __getitem__(self, i):
        return self.data[i]
    
    def __len__(self):
        return len(self.data)

    def get_pmid(self, sha):
        return self.sha2pmid.get(sha)

    def get_year(self, sha):
        return self.sha2year.get(sha)

    def get_sha(self, pmid):
        return self.pmid2sha.get(pmid)

    def count_identifiers(self):
        """Return counts for doi, pmcid and pubmed_id identifiers."""
        doi = 0
        pmcid = 0
        pubmed_id = 0
        for l in self.data:
            if l[3]: doi += 1
            if l[4]: pmcid += 1
            if l[5]: pubmed_id += 1
        print(doi, pmcid, pubmed_id)


class CovidData(object):

    """Class whose only goal is to create files named class-activators.txt,
    class-decreasers.txt, class-increasers.txt and class-inhibitors.txt. Those
    files have lists with all acitivators etcetera."""
    
    # TODO: this class has lost most functionality to MetaData and the remaining
    # code here probably should be in HarvardResults or be discarded.
    
    def __init__(self, results):
        self.results = results

    def collect_relations(self, relation_type):
        rels = []
        for sample in self.results.types.get(relation_type, []):
            if 'subj' in sample and 'obj' in sample:
                rels.append((relation_type, sample['subj']['name'], sample['obj']['name']))
        return rels

    def write_relation_classes(self):
        """Write classes of inhibitors, activators, increasers and decreasers to a
        couple of files in the working directory."""
        classes = [('Inhibition', 'inhibitor', 'inhibitee'),
                   ('Activation', 'activator', 'activatee'),
                   ('IncreaseAmount', 'increaser', 'increasee'),
                   ('DecreaseAmount', 'decreaser', 'decreasee')]
        for rel, agent, patient in classes:
            fh = open('class-%ss.txt' % agent, 'w')
            rels = self.collect_relations(rel)
            counter = Counter([r[2] for r in rels])
            args = {}
            for rel in rels:
                args.setdefault(rel[2], []).append(rel[1])
            for key in args.keys():
                args[key] = Counter(args[key])
            for ent, count in counter.most_common():
                if count >= 25:
                    fh.write("%s %s %ss\n" % (count, ent, agent))
                    for ent2, count2 in args[ent].most_common():
                        if count2 > 0:
                            fh.write('  %s %s\n' % (count2, ent2))
                    fh.write('\n')
            fh.close()


class CovidDoc(object):
    
    MINIMUM_TEXT_SIZE = 1000
    
    def __init__(self, fname, covid_data):
        """Collect data from the filename. Note how some information comes from
        the COVID meta data."""
        # TODO: process authors at this spot, like with the sections
        # TODO: process the body text a bit further too (no duplicate headers)
        self.json = json.load(open(fname))
        self.id = self.json['paper_id']
        self.pmid = covid_data.get_pmid(self.id)
        self.year = covid_data.get_year(self.id)
        self.title = self.json['metadata']['title']
        self.authors = self.json['metadata']['authors']
        self.abstract = []
        self.body_text = []
        for a in self.json['abstract']:
            self.abstract.append(a['text'])
        for s in self.json['body_text']:
            self.body_text.append((s['section'], s['text']))
        
    def is_complete(self):
        """A document is complete if it has a PMID a year and a title and enough
        text to work with."""
        return self.pmid and self.year and self.title and self.has_enough_text()
    
    def has_enough_text(self):
        abstract_size = sum([len(a) for a in self.abstract])
        body_size = sum([len(t) for h, t in self.body_text])
        return abstract_size + body_size > CovidDoc.MINIMUM_TEXT_SIZE


class Identifiers(object):

    identifiers = {}

    @classmethod
    def reset(cls):
        cls.identifiers = {}

    @classmethod
    def new_id(cls, tagname):
        cls.identifiers[tagname] = cls.identifiers.get(tagname, 0) + 1
        return "{}{:d}".format(tagname, cls.identifiers[tagname])


class HarvardResults(object):

    # We are not using the database references, but for future reference it uses
    # the following:
    #
    #   HGNC for gene names
    #   UP for Uniprot
    #   FPLX for the FamPlex namespace
    
    def __init__(self, json_file):
        self.fname = json_file
        self.results = json.load(open(self.fname))
        self.types = {}
        self.characterizations = {}
        for result in self.results:
            self.types.setdefault(result['type'], []).append(result)
        self._init_characterization()

    def __getitem__(self, i):
        return self.results[i]

    def __len__(self):
        return len(self.results)

    def print_types(self):
        """Show the list of relation types."""
        c = Counter([result['type'] for result in self.results])
        return c.most_common()
    
    def get_sample(self, result_type):
        """Return a single random sample result for a particular relation type."""
        return self.get_samples(result_type)

    def get_samples(self, result_type, n=1):
        """Return n random sample results for a particular relation type. """
        pool = self.types.get(result_type, [])
        return random.choices(pool, k=n)

    def _init_characterization(self):
        """Show what kind of arguments we have for each relation type."""
        for t in sorted(self.types.keys()):
            count = len(self.types[t])
            keys = []
            for r in self.types[t]:
                keys.extend(r.keys())
            c = Counter(keys)
            for k in ('type', 'matches_hash', 'id', 'belief'):
                if c[k] == count:
                    del(c[k])
            self.characterizations[t] = {'COUNT': count, 'ARGS': c }

    def collect_relations(self, metadata, n=100000):
        """Collect all relations and return them as a list of four-tuples with relation
        type, subject, object and an evidence list. This method requires access
        to an instace of CovidData. Four relation types are collected:
        Activation, Inhibition, IncreaseAmount and DecreaseAmount."""
        target_relations = ('Activation', 'Inhibition',
                            'IncreaseAmount', 'DecreaseAmount')
        relations = []
        for result in self.results[:n]:
            reltype = result['type']
            if reltype not in target_relations:
                continue
            try:
                sub = result['subj']['name']
                obj = result['obj']['name']
            except KeyError:
                # we want both a subject and an object
                continue
            evidence = []
            for e in result['evidence']:
                pmid = e.get('pmid')
                text = e.get('text')
                sha = metadata.get_sha(pmid)
                if sha is not None:
                    evidence.append((pmid, sha, text))
            # don't keep relations that come without evidence
            if evidence:
                relations.append((reltype, sub, obj, evidence))
        return relations

    def print_characterization(self):
        for relation_type in self.characterizations.keys():
            type_characterization = self.characterizations[relation_type]
            print(relation_type, type_characterization['COUNT'])
            for k,v in type_characterization['ARGS'].items():
                print('   ', k,v)
            print()

    def print_samples(self, relation_type, n, covid):
        """Given the Harvard results, print the specified number of samples for a
        particualr relation type."""
        for sample in self.get_samples(relation_type, n):
            print('type  ', sample['type'])
            if 'subj' in sample:
                self._print_arg('subj', sample['subj'])
            if 'obj' in sample:
                self._print_arg('obj ', sample['obj'])
            for e in sample['evidence']:
                pmid = e['pmid']
                sha = covid.pmid2sha.get(pmid)
                print('\n   ', pmid, sha)
                for line in textwrap.wrap(e['text'], width=80):
                    print('   ', line)
                print()

    def _print_arg(self, header, arg):
        print('%s   %-25s' %  (header, arg['name']), end='')
        if 'db_refs' in arg:
            refs = ' '.join(["%s=%s" % (k,v) for k,v in arg['db_refs'].items()])
            print('db_refs: { %s }' % refs)
        else:
            print()


def index_by_fname(relations):
    """Take a list of relations as created by HarvardResults.collect_relations() and
    return the relations indexed on the filename."""
    idx = {}
    for (reltype, sub, obj, evidence) in relations:
        for e in evidence:
            sha = e[1]
            fname = sha + '.json'
            idx.setdefault(sha, []).append((reltype, sub, obj, evidence))
    return idx


def reify_relations(relations):
    """Take a list of relations as created by HarvardResults.collect_relations() and
    return an index of reified reations where the object is folded into the
    relation. At the toplevel the index is keyed by the four relation types and
    at the second level on the reified relations where the value is a dictionary
    keyed on subjects with occurrences as values:

    "TNF-activator": {
       "IL12": [
       {
          "pmid": "21188201",
          "sha": "d3f7afa8b4d0f21b23ccb1135dec12356375f5cc",
          "text": "IL-12 stimulates production of IFNgamma and TNFalpha by T and natural..."
       },
       {
          "pmid": "19325820",
          "sha": "fd955af6faa247ced61613abf842cf3affac6bb2",
          "text": "As poly (ICLC) stimulates IL-12 and CpG ODN stimulates IL-6, -12 and..."
       }]

    """
    rels = { reltype: {} for reltype in RELTYPES }
    for rel in relations:
        reified_rel = "%s-%s" % (rel[2], translate_reltype_into_role(rel[0]))
        rels[rel[0]].setdefault(reified_rel, {})
        ds = [{'pmid': e[0], 'sha': e[1], 'text': e[2]} for e in rel[3]]
        rels[rel[0]][reified_rel].setdefault(rel[1], []).extend(ds)
    return rels


def print_relations(relations):
    """Print relations as created by HarvardResults.collect_relations()."""
    for (reltype, sub, obj, evidence) in relations:
        action = translate_reltype_into_action(reltype)
        print("%-20s  %-12s  %-25s  %3d documents" % (sub, action, obj, len(evidence)))


def print_index(idx, n=25):
    """Print the index that was created by index_by_fname()."""
    c = 0
    for fname in idx:
        c += 1
        if c > n:
            break
        print(fname)
        for rel in idx[fname][:5]:
            print('   ', rel)
        if len(idx[fname]) > 5:
            print('    ...')


class Converter(object):
    
    """Converts the JSON from a COVID file into a LIF document."""

    # TODO: add the directory of the sourcefile to the metadata
    # TODO: (this is to destinguish between the licenses)
    
    def __init__(self, infile, outfile, metadata):
        self.infile = infile
        self.outfile = outfile
        self.doc = CovidDoc(self.infile, metadata)

    def convert(self):
        print('Converting', os.path.basename(self.infile))
        if not self.doc.is_complete():
            print('skipping')
            return
        with open(self.outfile, 'w') as fh:
            self._setup()
            self._collect_metadata()
            self._add_abstract()
            self._add_sections()
            self._finish()
            
    def _setup(self):
        Identifiers.reset()
        self.p = 0
        self.lif = LIF()
        self.text = StringIO()
        self.view = View('docstruct')

    def _collect_metadata(self):
        self.lif.metadata['title'] = self.doc.title
        self.lif.metadata['sha'] = self.doc.id
        self.lif.metadata['pmid'] = self.doc.pmid
        self.lif.metadata['year'] = self.doc.year
        self.lif.metadata['authors'] = []
        for author in self.doc.authors:
            self.lif.metadata['authors'].append("%s %s" % (author['first'], author['last']))

    def _add_docelement_anno(self, docelement_type, p1, p2):
        self.view.add(
            Annotation(
                {'id': Identifiers.new_id('de'),
                 '@type': 'Section',
                 'start': p1,
                 'end': p2,
                 'features': {'section_type': docelement_type}}))
        
    def _add_abstract(self):
        # TODO: would like to add the section header
        # TODO: should make sure that the docelement ends not after the newlines
        abstract_p0 = self.p
        for text_str in self.doc.abstract:
            text_str += u"\n\n"
            chars = len(text_str)
            self.p += chars
            self.text.write(text_str)
        self._add_docelement_anno('Abstract', abstract_p0, self.p)

    def _add_sections(self):
        # TODO: add section header
        previous_header = None
        section_p0 = self.p
        for header_str, text_str in self.doc.body_text:
            text_str += u"\n\n"
            header_str += u"\n\n"
            chars = len(text_str)
            self.p += chars
            if header_str != previous_header:
                # fh.write(header_str)
                previous_header = header_str
            self.text.write(text_str)
            self._add_docelement_anno('Paragraph', section_p0, self.p)
            section_p0 = self.p

    def _finish(self):
        """Gather it all up and write output."""
        self.lif.text = Text(json_obj={'language': 'en', '@value': self.text.getvalue()})
        self.lif.views.append(self.view)
        self.lif.write(self.outfile, pretty=True)


class RelationImporter():

    def __init__(self, metadata_file, results_file, lif_dir, out_dir, n=99999):
        print('Loading metadata...')
        self.metadata = Metadata(metadata_file)
        print('Loading Harvard processing results...')
        self.results = HarvardResults(results_file)
        self.lif_dir = lif_dir
        self.out_dir = out_dir

    def convert(self, n=100000):
        print('Converting files...')
        self.relations = self.results.collect_relations(self.metadata)
        self.reified_rels = reify_relations(self.relations)
        self.filter_relobjs()
        self.invert_filtered_relobjs()
        #self.print_reified_rels_counts()
        #self.print_filtered_relobjs()
        #print(len(self.inverted_rels))
        for fname in os.listdir(self.lif_dir)[:n]:
            infile = os.path.join(self.lif_dir, fname)
            outfile = os.path.join(self.out_dir, fname)
            if outfile.endswith('.json'):
                outfile = outfile[:-4] + 'lif'
            self.convert_file(fname, infile, outfile)

    def convert_file(self, fname, infile, outfile):
        print(infile)
        lif = LIF()
        lif.text.value = None
        for relobj, subj in self.inverted_rels.get(fname, []):
            lif.metadata.setdefault(relobj, []).append(subj)
        lif.write(outfile, pretty=True)

    def filter_relobjs(self):
        self.filtered_rels = { reltype: {} for reltype in RELTYPES }
        for reltype in self.reified_rels:
            for rel_obj in self.reified_rels[reltype]:
                subjects = self.reified_rels[reltype][rel_obj]
                if not is_big(subjects.values(), reltype):
                    continue
                self.filtered_rels[reltype].setdefault(rel_obj, {'size': size_of(subjects.values()), 'data': {}})
                subjects_sorted = reversed(sorted([(len(v), k, v) for k, v in subjects.items()]))
                for (i, subj, val) in list(subjects_sorted)[:8]:
                    self.filtered_rels[reltype][rel_obj]['data'].setdefault(subj, []).extend(val)

    def invert_filtered_relobjs(self):
        """Create the inverted_rels index with the relations indexed on filenames."""
        self.inverted_rels = {}
        for reltype in self.reified_rels:
            for relobj in self.filtered_rels[reltype]:
                for subj in self.filtered_rels[reltype][relobj]['data']:
                    for e in self.filtered_rels[reltype][relobj]['data'][subj]:
                        self.inverted_rels.setdefault(e['sha'] + '.json', []).append((relobj, subj))

    def print_significant_rel_objs(self):
        for reltype in self.filtered_rels:
            print(">>> %s\n" % reltype)
            for relobj in self.filtered_rels[reltype]:
                print(relobj, self.filtered_rels[reltype][relobj]['size'])
                for subj in self.filtered_rels[reltype][relobj]['data']:
                    count = len(self.filtered_rels[reltype][relobj]['data'][subj])
                    print('    ', count, subj)
                print()

    def print_filtered_relobjs(self):
        for reltype in self.filtered_rels:
            print(">>> %s\n" % reltype)
            for relobj in self.filtered_rels[reltype]:
                print(relobj, self.filtered_rels[reltype][relobj]['size'])
                for subj in self.filtered_rels[reltype][relobj]['data']:
                    count = len(self.filtered_rels[reltype][relobj]['data'][subj])
                    print('    ', count, subj)
                print()

    def print_reified_rels_counts(self):
        for rt in self.reified_rels:
            print(rt, len(self.reified_rels[rt]))


def size_of(subj_sentences):
    return sum([len(sents) for sents in subj_sentences])


def is_big(subj_sentences, reltype):
    minimum_size = 25 if reltype.endswith('Amount') else 100
    return size_of(subj_sentences) >= minimum_size



if __name__ == '__main__':

    if sys.argv[1] == '--convert':
        metadata = sys.argv[2]
        data_dir = sys.argv[3]
        out_dir = sys.argv[4]
        n = int(sys.argv[5]) if len(sys.argv) > 5 else None
        convert_into_lif(metadata, data_dir, out_dir, n)

    elif sys.argv[1] == '--create-relations':
        metadata = sys.argv[2]
        results = sys.argv[3]
        out_file = sys.argv[4]
        create_relations_file(metadata, results, out_file)

    elif sys.argv[1] == '--import':
        metadata = sys.argv[2]
        results = sys.argv[3]
        lif_dir = sys.argv[4]
        out_dir = sys.argv[5]
        n = int(sys.argv[6]) if len(sys.argv) > 6 else 100000
        RelationImporter(metadata, results, lif_dir, out_dir).convert(n)

    else:
        pass
