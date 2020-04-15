"""load_index.py

Loading documents into the elastic search demo index.

Usage:

$ python load_index.py INDEX_NAME DIRECTORY

Load JSON documents from DIRECTORY into an index named INDEX_NAME.

Edit the HOST and PORT variables below if you do not need the defaults
(localhost:9200).

"""

import os
import sys
import codecs
import json

from elastic import Index


HOST = 'localhost'
PORT = 9200


def read_documents(document_directory):
    documents = []
    for directory_element in sorted(os.listdir(document_directory)):
        if directory_element.endswith('.json'):
            print(directory_element)
            fname = os.path.join(document_directory, directory_element)
            json_obj = json.loads(codecs.open(fname, encoding='utf8').read())
            documents.append(json_obj)
    return documents


if __name__ == '__main__':

    if len(sys.argv) > 2:
        index_name = sys.argv[1]
        source_directory = sys.argv[2]
    else:
        exit('ERROR: missing arguments\nUsage: python load_index.py INDEX_NAME DIRECTORY\n')
    if len(sys.argv) > 3:
        mapping_fname = sys.argv[3]
    else:
        mapping_fname = None

    docs = read_documents(source_directory)
    idx = Index(index_name, host=HOST, port=PORT)
    if mapping_fname is not None:
        idx.es.indices.delete(index=index_name, ignore=[400, 404])
        mappings = json.load(open(mapping_fname))
        idx.es.indices.create(index_name, body=mappings)
    idx.load(docs)

    print("Loading documents into the index...")
