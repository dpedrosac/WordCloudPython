#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from Bio import Entrez, Medline
import requests
from dependencies import ROOTDIR
from io import BytesIO
from tqdm import tqdm

s = requests.Session()  # Create a session object for making persistent HTTP connections


def search_medlinePIDparallel(id_complete, entered_mail):
    """Medline search with a list of IDs, so that the server does not block searches"""

    Entrez.email = entered_mail
    response = s.post(url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi",
                      data={"db": "pubmed", "id": ",".join(id_complete)})  # Use the session object to make the API call

    fp = BytesIO(response.content)  # Create a file-like object from the response text
    result = Entrez.read(fp)  # Pass the file-like object to the Entrez.read() function

    return result


class GetAbtractsFromPubmed:
    def __init__(self, entered_PMID, entered_email='pedrosac@staff.uni-marburg.de', _debug=False):

        self.email = entered_email
        self.process_abstracts(entered_PMID, entered_email)

    def process_abstracts(self, set_PID, entered_mail, batch_size=4):
        """Function extracting content of abstracts using the Entrez routines and saving all abstracts to datafile
         as recording per row"""

        set_PID = set(set_PID)
        print(f'\nThe sum of all abstracts is {len(set_PID)}')
        print(f'\t ...extracting Abstracts to *.txt-files\n')

        result = search_medlinePIDparallel(set_PID, entered_mail)
        webenv, query_key = result["WebEnv"], result["QueryKey"]

        with tqdm(total=len(set_PID)) as pbar:
            for start in range(0, len(set_PID), batch_size):
                processed_records = self.process_batch(set_PID, start, batch_size, webenv, query_key)
                pbar.update(processed_records)


    def process_batch(self, set_PID, start, batch_size, webenv, query_key):
        """Function processing a batch of records"""

        fetch_handle = Entrez.efetch(db="pubmed",
                                     rettype="Medline", retmode="text",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv,
                                     query_key=query_key)

        result = Medline.parse(fetch_handle)

        processed_abstracts = 0
        for record in result:
            if 'AB' in record:  # text is only saved when "AB" is present
                save_filename = os.path.join(ROOTDIR, 'data', 'all_abstracts.txt')

                with open(save_filename, 'a+') as file_object:
                    file_object.seek(0)
                    data = file_object.read(100)
                    if len(data) > 0:
                        file_object.write('\n')
                    file_object.write(record["AB"])
                file_object.close()

            processed_abstracts += 1

        return processed_abstracts


def main(email):
    """Creates the wordcloud by calling several functions """
    list_PID = list(['22173967', '22526237', '22809566', '23056463', '23487540',
                     '23778146', '24467817', '24661791', '24848641', '25339758',
                     '27091412', '27501132', '28459950', '28645854', '28729930',
                     '28754506', '29560414', '29560414', '29606416', '30052807',
                     '30827864', '31654235', '31654235', '32057084', '32450842',
                     '32450842', '32853481', '32854328', '33267627', '33638213',
                     '33958263', '34209024', '34209024', '34352356', '34404706',
                     '34475849', '34911202', '35259208', '35292120', '35629226',
                     '35732679', '35887553'])

    if not email: # Check if an email address was provided
        print('No email was provided in GetData.py')
        return

    GetAbtractsFromPubmed(list_PID, entered_email=email, _debug=False)  # perform medline search and save to txt-file


if __name__ == '__main__':
    main(email='pedrosac@staff.uni-marburg.de')

# close the session when you are done making API calls
s.close()
