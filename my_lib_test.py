import pytest
import pandas as pd
import my_lib
import xml.etree.ElementTree as Et

@pytest.fixture
def root():
    tree = Et.parse('drugbank_partial.xml')
    root = tree.getroot()
    return root

@pytest.fixture
def namespace():
    return {'ns': 'http://www.drugbank.ca'}


def test_int_to_db_string():
    assert my_lib.int_to_db_string(123) == "DB00123"
    assert my_lib.int_to_db_string(0) == "DB00000"
    assert my_lib.int_to_db_string(1) == "DB00001"
    assert my_lib.int_to_db_string(2) == "DB00002"
    assert my_lib.int_to_db_string(20000) == "DB20000"


def test_extract_drugs(root, namespace):
    result  = my_lib.extract_drugs(root, namespace)

    assert isinstance(result , pd.DataFrame)
    assert result .shape[0] == 100
    assert result .shape[1] == 6
    assert result .loc['DB00001', 'name'] == 'Lepirudin'
    assert len(result .loc['DB00001', 'food interactions']) == 1


def test_extract_synonyms(root, namespace):
    result  = my_lib.extract_synonyms(root, namespace)
    assert isinstance(result , pd.DataFrame)
    assert result .shape[0] == 100
    assert result .loc['DB00001', 'name'] == 'Lepirudin'
    assert len(result .loc['DB00100', 'synonyms']) == 7


def test_extract_products(root, namespace):
    result = my_lib.extract_products(root, namespace)
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 4584
    assert result .shape[1] == 9

def test_extract_pathway_ids(root, namespace):
    result = my_lib.extract_pathway_ids(root, namespace, my_lib.extract_drugs(root, namespace))
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 102
    assert result.shape[1] == 1


def test_extract_prices(root, namespace):
    result = my_lib.extract_prices(root, namespace)
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 499
    assert result.shape[1] == 3


def test_filter_prices(root, namespace):
    result = my_lib.filter_prices(my_lib.extract_prices(root, namespace))
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 176
    assert result.shape[1] == 4


def test_extract_drug_approval_status(root, namespace):
    result = my_lib.extract_drug_approval_status(root, namespace)
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 100
    assert result.shape[1] == 7

def test_summarise_drug_approval_status(root, namespace):
    result = my_lib.summarise_drug_approval_status(my_lib.extract_drug_approval_status(root, namespace))
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 5
    assert result.shape[1] == 2

def test_extract_drug_interactions(root, namespace):
    result = my_lib.extract_drug_interactions(root, namespace)
    assert isinstance(result , pd.DataFrame)
    assert result.shape[0] == 50688
    assert result.shape[1] == 5