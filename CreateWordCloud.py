#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import nltk
import re
import pandas as pds
import Utils
from dependencies import ROOTDIR


class ProcessResults:
    def __init__(self):
        data_directory = os.path.join(ROOTDIR, 'data')
        self.remove_stopwords_punctuation(data_directory)

    def remove_stopwords_punctuation(self, data_directory):
        import string
        from nltk.stem import WordNetLemmatizer
        nltk.download('stopwords')
        nltk.download('punkt')
        nltk.download('wordnet')
        nltk.download('averaged_perceptron_tagger')

        print('Start removing stop words')
        if not data_directory:
            print('\t ...no directory provided, ending the script!')
            return

        stop_words = set(nltk.corpus.stopwords.words('english'))
        stop_words = stop_words.union(Utils.HelperFunctions.stopwords_general())
        stop_words = stop_words.union(Utils.HelperFunctions.stopwords_medicine())
        punctuation = list(string.punctuation)

        replace_terms = {"patients": ["patients", "patient"],
                         "parents": ["parents", "parent"],
                         "carers": ["carers", "carer"],
                         "paediatric": ["paediatric", "pediatric"],
                         "Patient safety": ["Patient safety", "patient-safety"],
                         "patients' safety": ["patients' safety", "patient-safety"],
                         "parkinson's disease": ["parkinson's' disease", "parkinsons-disease"],
                         "parkinson's disease": ["parkinson syndrome", "parkinsons-disease"],
                         "patient's safety": ["patient's' safety", "patient-safety"],
                         "patient safety": ["patient safety", "patient-safety"],
                         "terminally-ill": ["terminally ill", "terminally-ill"],
                         "follow-up": ["follow up", "follow-up"],
                         "end-of-life": ["end of life", "end-of-life"],
                         "palliative care": ["palliative care", "palliative-care"],
                         "palliative medicine": ["palliative medicine", "palliative-medicine"],
                         "patient care": ["patient care", "patient-care"],
                         "hospice care": ["hospice care", "hospice-care"],
                         "cancer care": ["cancer care", "cancer-care"],
                         "end-of-life care": ["end-of-life care", "end-of-life-care"],
                         "acp": [" acp ", " advanced-care-planning "],
                         "quality of life": ["quality of life", "quality-of-life"],
                         "qol": ["qol", "quality-of-life"],
                         "eol": [" eol ", " end-of-life "],
                         "levodopa": ["ldopa", "l-dopa"],
                         "eolc": [" eolc ", " end-of-life-care "],
                         "gp": [" gp ", " general-practitioners "],
                         "gps": [" gps ", " general-practitioners "],
                         "icu": [" icu ", " intensive-care-unit "],
                         "dnr": ["dnr", "do-not-resuscitate"],
                         "DBS": ["dbs", "deep-brain-stimulation"],
                         "deep brain stimulation": ["brain stimulation", "brain-stimulation", "deep brain",
                                                    "deep-brain"]
                         }

        allfiles = os.listdir(data_directory)
        file_list = [file_id for file_id in allfiles if file_id.startswith("all_abstract")]
        [os.remove(os.path.join(data_directory, file_id)) for file_id in allfiles if file_id.startswith("pp_")]

        print("\n\t...removing all stop words, punctuation and replacing words")
        Utils.HelperFunctions.printProgressBar(0, len(file_list), prefix='Progress:',
                                               suffix='Complete', length=50, decimals=2)

        lemmatizer = WordNetLemmatizer()
        for idx, file_id in enumerate(file_list, start=1):
            Utils.HelperFunctions.printProgressBar(idx, len(file_list), prefix='Progress:',
                                                   suffix='Complete', length=50, decimals=2)

            filename, file_extension = os.path.splitext(file_id)
            save_filename = "pp_" + filename + file_extension
            with open(os.path.join(data_directory, file_id), 'r+') as file_object:
                text = file_object.read()
                for item, values in replace_terms.items():
                    if item in text:
                        text = text.replace(values[0], values[1])

                text = re.sub(r'\d+', '', text)
                words = [i.strip("".join(punctuation)) for i in nltk.word_tokenize(text) if i not in punctuation]

                words_clean = [lemmatizer.lemmatize(r, Utils.HelperFunctions.get_wordnet_pos(r)) for r in words
                               if not (r.lower() in stop_words or r in string.punctuation or r.isdigit())]

                appendFile = open(os.path.join(data_directory, save_filename), 'a')
                appendFile.write('\n'.join(words_clean))
                appendFile.close()

    def mostCommonWords(self):
        """this function returns a list of the most common words found in the text and saves them additionally per
        year of the abstract"""

        if not self.data_directory:
            print('No directory was provided, ending the script')
            return

        allfiles = os.listdir(self.data_directory)
        file_list = [list_pp for list_pp in allfiles if list_pp.startswith("pp_")]

        if not file_list:
            print('No preprocessed data was found, running the part where stop words and punctuation is removed again')
            self.remove_stopwords_punctuation()

        print("\nLoading all available tokens from files")
        Utils.HelperFunctions.printProgressBar(0, len(file_list), prefix='Progress:',
                                               suffix='Complete', length=50, decimals=2)
        allWords = []
        words_per_year = []
        for idx, file_id in enumerate(file_list, start=1):
            Utils.HelperFunctions.printProgressBar(idx, len(file_list), prefix='Progress:',
                                                   suffix='Complete', length=50, decimals=2)
            with open(os.path.join(self.data_directory, file_id), 'r+') as text_data:
                words_temp = text_data.read()
                try:
                    words_per_year.append(nltk.tokenize.word_tokenize(words_temp))
                    allWords.extend(nltk.tokenize.word_tokenize(words_temp))
                except LookupError:
                    nltk.download('punkt')
                    words_per_year.append(nltk.tokenize.word_tokenize(words_temp))
                    allWords.extend(nltk.tokenize.word_tokenize(words_temp))

        allWordDist = nltk.FreqDist(w.lower() for w in allWords)
        WordperYearDist = [nltk.FreqDist([x.lower() for x in year_list]) for year_list in words_per_year]

        df_allWords = pds.DataFrame(allWordDist.most_common(len(allWordDist)), columns=['token', 'count'])
        df_WordsPerYear = [pds.DataFrame(year_list.most_common(len(year_list)), columns=['token', 'count'])
                           for idx, year_list in enumerate(WordperYearDist)]

        print("\nSaving the obtained dataframes as csv-files")
        df_allWords.to_csv(os.path.join(self.data_directory, "results_complete.csv"), sep='\t')

        years = [re.findall("\d+", file_id)[0] for file_id in file_list]
        Utils.HelperFunctions.printProgressBar(0, len(years), prefix='Progress:',
                                               suffix='Complete', length=50, decimals=2)

        # Saves the results according to the years of publication
        for idx, publication_date in enumerate(years):
            Utils.HelperFunctions.printProgressBar(idx, len(years), prefix='Progress:',
                                                   suffix='Complete', length=50, decimals=2)
            filename2save = os.path.join(self.data_directory, "results" + publication_date + ".csv")
            df_WordsPerYear[idx].to_csv(filename2save, sep='\t')

    @staticmethod
    def get_wordnet_pos(word):
        from nltk.corpus import wordnet

        """Map POS tag to first character lemmatize() accepts"""
        tag = nltk.pos_tag([word])[0][1][0].upper()
        tag_dict = {"J": wordnet.ADJ,
                    "N": wordnet.NOUN,
                    "V": wordnet.VERB,
                    "R": wordnet.ADV}

        return tag_dict.get(tag, wordnet.NOUN)


class GenerateOutput:
    def __init__(self):
        self.create_wordcloud()

    @staticmethod
    def create_wordcloud():
        from wordcloud import WordCloud
        import matplotlib.pyplot as plt
        with open(os.path.join(ROOTDIR, 'data', 'pp_all_abstracts.txt'), 'r') as file_object:
            text = file_object.read()

        wordcloud = WordCloud(background_color="white", width=2000, height=1000, max_words=200).generate(text)

        def black_color_func(word, font_size, position, orientation, random_state=None, **kwargs):
            return "hsl(0,100%, 1%)"
        wordcloud.recolor(color_func=black_color_func)

        plt.imshow(wordcloud, interpolation="bilinear")
        plt.axis('off')
        plt.show()


if __name__ == '__main__':
    ProcessResults()
    GenerateOutput()
