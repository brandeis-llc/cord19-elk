"""elastic.py

Module with some convenience code for acessing an Elastic Search index.

"""

from pprint import pprint
from collections import Counter
import json

from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch.exceptions import NotFoundError


class Index(object):

    def __init__(self, index_name, host='localhost', port=9200, index_elements=None):
        self.index = index_name
        self.es = Elasticsearch([{'host': host, 'port': port}])
        if index_elements is not None:
            self.load(index_elements)

    def to_bulk_iterable(self, elements):
        for i, element in enumerate(elements):
            docid = element.get('docid')
            identifier = i if docid is None else docid
            print(docid)
            yield {
                "_type":"_doc",
                "_id":identifier,
                "_index": self.index,
                "_source": element } 

    def load(self, elements):
        helpers.bulk(self.es, self.to_bulk_iterable(elements))

    def get(self, message, doc_id, dribble=False):
        print("\n{}".format(message))
        try:
            doc = self.es.get(index=self.index, id=doc_id)
            if dribble:
                pprint(doc['_source'])
            return doc
        except NotFoundError as e:
            print(e)

    def search(self, message, query, dribble=False):
        print("\n{}".format(message))
        result = Result(self.es.search(index=self.index, body=query))
        result.print_sources(dribble)
        return result


class Result(object):

    """Class to wrap an ElasticSearch result."""
    
    def __init__(self, result):
        self.result = result
        self.hits = [Hit(hit) for hit in self.result['hits']['hits']]
        self.total_hits = self.result['hits']['total']['value']
        self.sources = [hit.source for hit in self.hits]

    def write(self):
        fname = "{:04d}.txt".format(nextint())
        with open(fname, 'w', encoding='utf8') as fh:
            fh.write(json.dumps(self.result, sort_keys=True, indent=4))

    def pp(self):
        print("\n    Number of hits: {:d}".format(self.total_hits))
        for hit in self.hits:
            print("    {}  {:.4f}  {}".format(hit.docid, hit.score, hit.docname[:80]))

    def print_sources(self, dribble):
        if dribble:
            sources = self.sources
            print('   Got {:d} hits'.format(self.total_hits))
            for source in self.sources:
                print('   {}'.format(source))


class Hit(object):

    def __init__(self, hit):
        self.hit = hit
        self.id = hit['_id']
        self.score = hit['_score']
        self.source = hit['_source']
        self.docid = self.source.get('docid')
        self.docname = self.source.get('docname')


def nextint(data=Counter()):
    data['count'] += 1
    return data['count']



if __name__ == '__main__':

    entities = [
        {"name": {"first":"john", "last":"doe"}, "age": 27, "interests": ['sports','music']},
        {"name": {"first": "jane", "last": "doe"}, "age": 27, "interests": ["rap music yeah"]},
        {"name": {"first": "june", "last": "doe"}, "age": 65, "interests": ["forestry"]}]

    idx = Index('test_entities', entities)

    queries = [

        ["Matching all elements",
         {'query': {'match_all': {}}}],

        ["Find the young ones",
         {'query': {'match': {'age': 27}}}],

        ["Find the Does",
         {'query': {'match': {'name.last': 'doe'}}}],

        ["Find Jane Doe",
         {'query':
          {'bool':
           {'must': [
               {'match': {'name.first': 'jane'}},
               {'match': {'name.last': 'doe'}}]}}}],

        ["Find Sripi Palaver",
         {'query':
          {'bool':
           {'must': [
               {'match': {'name.first': 'sripi'}},
               {'match': {'name.last':'palaver'}}]}}}],

        ["Find rap or music",
         {'query': {'match': {'interests': 'music rap'}}}],

        ["Find \"rap music\"",
         {'query': {'match_phrase': {'interests': 'rap music'}}}]
        ]

    idx.get("Retrieving document with id=1", 1, dribble=True)
    for message, query in queries:
        idx.search(message, query, dribble=True)
