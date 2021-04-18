#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from Bio import Entrez, Medline
import Utils
from dependencies import ROOTDIR


def search_medlinePID(id_complete, email):
    """Medline search with a list of IDs, so that the server does not block searches"""

    Entrez.email = email
    handle = Entrez.epost(db="pubmed", id=",".join(id_complete))
    result = Entrez.read(handle)
    handle.close()

    try:
        return result
    except Exception as e:
        raise IOError(str(e))
    finally:
        handle.close()


class GetAbtractsFromPubmed:
    def __init__(self, list_PID, _email='pedrosac@staff.uni-marburg.de', _debug=False):

        self.email = _email
        self.process_abstracts(list_PID, _email)

    @staticmethod
    def process_abstracts(set_PID, email, batch_size=50):
        """Function extracting content of abstracts using the Entrez routines and saving all abstracts to datafile
         as recording per row"""

        set_PID = set(set_PID)
        print('\nThe sum of all abstracts is {}'.format(len(set_PID)))
        print('\t ...extracting Abstracts to *.txt-files\n')

        result = search_medlinePID(set_PID, email)
        webenv, query_key = result["WebEnv"], result["QueryKey"]

        abstract_unavailable, idx = [0 for _ in range(2)]
        Utils.HelperFunctions.printProgressBar(idx, len(set_PID),
                                               prefix='Progress:', suffix='Complete', length=50, decimals=2)

        for start in range(0, len(set_PID), batch_size):
            fetch_handle = Entrez.efetch(db="pubmed",
                                         rettype="Medline", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv,
                                         query_key=query_key)

            result = Medline.parse(fetch_handle)
            for record in result:
                idx += 1
                Utils.HelperFunctions.printProgressBar(idx, len(set_PID), prefix='Progress:', suffix='Complete',
                                                       decimals=2, length=50)

                if 'AB' in record:  # text is only saved when "AB" is present
                    save_filename = os.path.join(ROOTDIR, 'data', 'all_abstracts.txt')

                    with open(save_filename, 'a+') as file_object:
                        file_object.seek(0)
                        data = file_object.read(100)
                        if len(data) > 0:
                            file_object.write('\n')
                        file_object.write(record["AB"])
                    file_object.close()

            print('\nabstracts:{}, handlers:{}, none_available:{}'.format(idx, len(set_PID), abstract_unavailable))


if __name__ == '__main__':
    list_PID = list(['29560414',
                     '31654235',
                     '33267627',
                     '32450842',
                     '30827864',
                     '28645854',
                     '27091412',
                     '29606416',
                     '23487540',
                     '30052807',
                     '22173967',
                     '28459950',
                     '28729930',
                     '28754506',
                     '24467817',
                     '24661791',
                     '25339758',
                     '27501132',
                     '24848641',
                     '23056463',
                     '23778146',
                     '22809566',
                     '22526237',
                     '32450842',
                     ])

    GetAbtractsFromPubmed(list_PID, _email='pedrosac@staff.uni-marburg.de', _debug=False)