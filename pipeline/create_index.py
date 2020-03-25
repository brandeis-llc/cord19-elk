"""create_index_docs.py

Merge information from different views and export JSON files that can be loaded
into ElasticSearch.

Usage:

$ python create_index_docs.py -d DATA_DIR -f FILELIST (-b BEGIN) (-e END) (--crash)

Directories:

lif   LIF files created from the Covid data
top   topics
har   relations from Harvard results
ela   output

Processes about 40 Covid documents per second.

"""

import os, sys, json
from pprint import pformat
from collections import Counter

from lif import LIF, Annotation
from utils import time_elapsed, elements, print_element, get_options



@time_elapsed
def create_documents(data_dir, filelist, start, end, crash=False):
    print("$ python3 %s\n" % ' '.join(sys.argv))
    ela_dir = os.path.join(data_dir, 'ela')
    if not os.path.exists(ela_dir):
        os.mkdir(ela_dir)
    for n, fname in elements(filelist, start, end):
        print_element(n, fname)
        if crash:
            create_document(data_dir, fname)
        else:
            try:
                create_document(data_dir, fname)
            except Exception as e:
                print('ERROR:', Exception, e)


def create_document(data_dir, fname):
    # the subdir is really the document identifier
    lif_file = os.path.join(data_dir, 'lif', fname)
    top_file = os.path.join(data_dir, 'top', fname[:-4] + 'lif')
    har_file = os.path.join(data_dir, 'har', fname[:-4] + 'lif')
    if not os.path.exists(lif_file):
        print('Skipping...  %s' % fname)
    else:
        doc = Document(fname, data_dir, lif_file, top_file, har_file)
        doc.write(os.path.join(data_dir, 'ela'))


def fix_view(identifier, view):
    annos = view.id['annotations']
    view.id = identifier
    view.annotations = []
    for a in annos:
        view.annotations.append(Annotation(a))


class Document(object):

    def __init__(self, fname, data_dir, lif_file, top_file, har_file):

        """Build a single LIF object with all relevant annotations. The annotations
        themselves are stored in the Annotations object in self.annotations."""
        self.id = fname
        self.fname = fname
        self.data_dir = data_dir
        self.lif = LIF(json_file=lif_file)
        self.top = LIF(json_file=top_file)
        self.har = LIF(json_file=har_file)
        # NOTE: no idea why this was needed
        fix_view('doc', self.lif.views[0])
        fix_view('top', self.top.views[0])
        self.lif.views.append(self.top.views[0])
        #print(self.lif)
        self.annotations = Annotations(self.id, fname, doc=self, text=self.lif.text.value)
        # NOTE: commented out for now because it is not clear we need the text
        # NOTE: we will need it when we need to do full text searches for the demo
        self.annotations.text = self.lif.text.value
        self._collect_authors()
        self._collect_topics()
        self._collect_relations()

    def get_view(self, identifier):
        return self.lif.get_view(identifier)

    def _collect_authors(self):
        """Just get the authors from the metadata and put them in the index."""
        self.annotations.authors = self.lif.metadata['authors']

    def _collect_topics(self):
        """Collect the topics and put them on a list in the index."""
        view = self.get_view("top")
        for annotation in view.annotations:
            if annotation.type.endswith('SemanticTag'):
                topic_name = annotation.features['topic_name']
                self.annotations.topics.append(topic_name)
                for topic_element in topic_name.split():
                    self.annotations.topic_elements.append(topic_element)
        self.annotations.topic_elements = sorted(set(self.annotations.topic_elements))

    def _collect_relations(self):
        added = False
        for relobj, subjs in self.har.metadata['relations'].items():
            for subj in subjs:
                if relobj in self.annotations.relations:
                    #print(relobj, subj)
                    added = True
                    self.annotations.relations[relobj].append(subj)
        #if added:
        #    print(self.annotations.relations)
        
    def write(self, dirname):
        self.annotations.write(os.path.join(dirname, self.fname),
                               self.lif.metadata["year"])

    def pp(self, prefix=''):
        views = ["%s:%d" % (view.id, len(view)) for view in self.lif.views]
        print("%s<Document id=%s '%s'>" % (prefix, self.id, self.fname))
        print("    <Views %s>" % ' '.join(views))
        print("    %s\n" % self.annotations)


class Annotations(object):

    """Object that holds all annotations for a file as well as the text of the
    document. Annotations include (1) metadata like authors and topics, which
    are not associated with offsets, (2) entities and events, which all are
    associated with text positions, and (3) relations, which currently have a
    special status in that they are the only complex annotation."""

    RELATIONS = ['TNF-activator',
                 'CD4-activator',
                 'IFNB1-activator',
                 'apoptotic process-activator',
                 'NFkappaB-activator',
                 'immune response-activator',
                 'cell death-activator 570',
                 'CD8-activator',
                 'Interferon-activator',
                 'inflammatory response-activator',
                 'autophagy-activator',
                 'endocytosis-activator',
                 'IL10-activator',
                 'IFNG-activator',
                 'innate immune response-activator',
                 'IRF3-activator',
                 'Interferon-inhibitor',
                 'replication-inhibitor',
                 'EIF2AK2-inhibitor',
                 'apoptotic process-inhibitor',
                 'TNF-inhibitor',
                 'NFkappaB-inhibitor',
                 'translation-inhibitor',
                 'autophagy-inhibitor',
                 'IFNB1-inhibitor',
                 'cell death-inhibitor',
                 'inflammatory response-inhibitor',
                 'cell population proliferation-inhibitor']

    def __init__(self, docid, fname, doc=None, text=None):
        self.docid = docid
        self.fname = fname
        self.doc = doc
        self.text = text
        self.authors = []
        self.year = None
        self.topics = []
        self.topic_elements = []
        self.relations = {}
        for rel in Annotations.RELATIONS:
            self.relations[rel] = []
        self.text = None

    def write(self, fname, year=None):
        """Writes the document with the search fields to a json file."""
        json_object = {
            "text": self.text,
            "docid": self.docid,
            "docname": self.fname,
            "year": year,
            "author": self.authors,
            "topic": self.topics,
            "topic_element": self.topic_elements,
        }
        for relobj, subj in self.relations.items():
            json_object[relobj] = subj
        with open(fname, 'w', encoding='utf8') as fh:
            fh.write(json.dumps(json_object, sort_keys=True, indent=4))

    def pp(self, indent=''):
        print("%s%s\n" % (indent, self))


if __name__ == '__main__':

    data_dir, filelist, start, end, crash = get_options()
    create_documents(data_dir, filelist, start, end, crash=crash)

