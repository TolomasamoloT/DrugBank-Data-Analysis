import pandas as pd
from re import findall as re_findall
import xml.etree.ElementTree as Et


# a library of functions used in <name here>, in a separate file for testing via pytest


def extract_drugs(root, namespace) -> pd.DataFrame:
    drugs_data = [
        {
            'id': drug.find('ns:drugbank-id', namespace).text,
            'name': drug.find('ns:name', namespace).text,
            'description': drug.find('ns:description', namespace).text,
            'state': drug.find('ns:state', namespace).text,
            'indication': drug.find('ns:indication', namespace).text,
            'mechanism of action': drug.find('ns:mechanism-of-action', namespace).text,
            'food interactions': [
                interaction.text.strip()
                for interaction in drug.find('ns:food-interactions', namespace)
                if interaction.text and interaction.text.strip()
            ] if drug.find('ns:food-interactions', namespace) is not None else None
        }
        for drug in root
    ]

    drugs_df = pd.DataFrame(drugs_data).set_index('id')
    drugs_df = drugs_df.convert_dtypes(convert_string=True)
    return drugs_df


def extract_synonyms(root, namespace) -> pd.DataFrame:
    synonyms_data = [
        {
            'id': drug.find('ns:drugbank-id', namespace).text,
            'name': drug.find('ns:name', namespace).text,  # this might be redundant
            'synonyms': [
                synonym.text.strip()
                for synonym in drug.find('ns:synonyms', namespace)
                if synonym.text and synonym.text.strip()
            ]
        }
        for drug in root
    ]

    synonyms_df = pd.DataFrame(synonyms_data).set_index('id')
    synonyms_df = synonyms_df.convert_dtypes(convert_string=True)
    return synonyms_df


def extract_products(root, ns) -> pd.DataFrame:
    products_data = [
        {
            'id': drug.find('ns:drugbank-id', ns).text,
            'name': product.find('ns:name', ns).text,
            'labeller': product.find('ns:labeller', ns).text,
            'ndc-product-code': product.find('ns:ndc-product-code', ns).text,
            'dosage-form': product.find('ns:dosage-form', ns).text,
            'route': product.find('ns:route', ns).text,
            'strength': product.find('ns:strength', ns).text,
            'country': product.find('ns:country', ns).text,
            'source': product.find('ns:source', ns).text
        }
        for drug in root
        for product in drug.find('ns:products', ns)
    ]

    products_df = pd.DataFrame(products_data).convert_dtypes(convert_string=True)
    return products_df


def extract_pathways(root, ns) -> pd.DataFrame:
    pathways_data = [
        {
            'smpdb-id': pathway[0].text,
            'name': pathway[1].text,
            'category': pathway[2].text,

            # # alternatively:
            # 'smpdb-id': pathway.find('ns:smpdb-id', ns).text,
            # etc ...
        }
        for pathway in root.findall('.//ns:pathway', ns)
    ]

    pathways_df = pd.DataFrame(pathways_data).convert_dtypes(convert_string=True)
    return pathways_df


def append_pathway_drugs(root, ns, pathways: pd.DataFrame) -> pd.DataFrame:
    pathway_drugs = [
        [
            drug.find('ns:name', ns).text
            for drug in pathway.find('ns:drugs', ns).iter()
            if drug.find('ns:name', ns) is not None
        ]
        for pathway in root.findall('.//ns:pathway', ns)
    ]

    pathways = pathways.copy()
    pathways['drugs'] = pathway_drugs

    return pathways


def extract_pathway_ids(root, ns, drugs:pd.DataFrame) -> pd.DataFrame:
    pathway_ids = [
        [
            drug.find('ns:drugbank-id', ns).text
            for drug in pathway.find('ns:drugs', ns).iter()
            if drug.find('ns:drugbank-id', ns) is not None
        ]
        for pathway in root.findall('.//ns:pathway', ns)
    ]

    pathways_id_df = pd.DataFrame()
    pathways_id_df['id'] = pathway_ids
    pathway_id_df = pathways_id_df.explode('id').reset_index(drop=True)
    pathways_id_df = pathway_id_df.convert_dtypes(convert_string=True)
    pathways_id_df = pathways_id_df['id'].value_counts()

    pathways_id_df = pathways_id_df.to_frame()

    drugs = drugs.copy()

    merged = pathways_id_df.merge(drugs, left_index=True, right_index=True, how='outer')
    merged['count'] = merged['count'].fillna(0).astype(int)
    merged = merged[['count']]
    return merged


def explode_pathways(pathways: pd.DataFrame) -> pd.DataFrame:
    boom = pathways.explode('drugs').drop(columns=['smpdb-id', 'name', 'category']).reset_index(drop=True)
    boom.rename(columns={'drugs': 'drug'}, inplace=True)
    return boom


def extract_targets(root, ns) -> pd.DataFrame:
    targets = []

    for drug in root:
        for target in drug.find('ns:targets', ns):
            row = {
                'drug id': drug.find('ns:drugbank-id', ns).text,
                'target id': target.find('ns:id', ns).text,
            }
            poly = target.find('ns:polypeptide', ns)
            if poly is not None:
                row['source'] = poly.get('source')
                row['source id'] = poly.get('id')
                row['polypeptide name'] = poly.find('ns:name', ns).text
                row['gene name'] = poly.find('ns:gene-name', ns).text

                # GenAtlas // just getting the target.gene-name would probably be better since it's the same
                for e_id in poly.find('ns:external-identifiers', ns):
                    if e_id[0].text == 'GenAtlas':
                        row['GenAtlas ID'] = e_id[1].text
                        break

                row['chromosome location'] = poly.find('ns:chromosome-location', ns).text
                row['cellular location'] = poly.find('ns:cellular-location', ns).text

            targets.append(row)

    return pd.DataFrame(targets).convert_dtypes(convert_string=True)


def extract_drug_approval_status(root, ns) -> pd.DataFrame:
    drug_approval_status = [
        {
            'drug id': drug.find('ns:drugbank-id', ns).text,
            'name': drug.find('ns:name', ns).text,
            'approved': 'approved' in statuses,
            'withdrawn': 'withdrawn' in statuses,
            'experimental': 'experimental' in statuses,
            'investigational': 'investigational' in statuses,
            'vet_approved': 'vet_approved' in statuses
        }
        for drug in root
        for statuses in [{t.text.lower() for t in drug.find("ns:groups", ns)}]
    ]

    return pd.DataFrame(drug_approval_status).convert_dtypes(convert_string=True, convert_boolean=True)


def summarise_drug_approval_status(drug_approval_status: pd.DataFrame) -> pd.DataFrame:
    # True -> 1, False -> 0
    numeric = drug_approval_status.copy().iloc[:, 2:].astype('int')

    # sum the columns
    summary = numeric.sum().reset_index()
    summary.columns = ["status", "number of drugs"]

    return summary


def extract_drug_interactions(root, ns) -> pd.DataFrame:
    drug_interactions = [
        {
            'drug name': drug.find("ns:name", ns).text,
            'drug id': drug.find("ns:drugbank-id", ns).text,
            'interacts with': interaction.find("ns:name", ns).text,
            'interactee id': interaction.find("ns:drugbank-id", ns).text,
            'interaction description': interaction.find("ns:description", ns).text,
        }
        for drug in root
        for interaction in drug.find("ns:drug-interactions", ns)
    ]

    return pd.DataFrame(drug_interactions).convert_dtypes(convert_string=True)


def extract_prices(root, ns) -> pd.DataFrame:
    prices = [
        {
            'description': price.find('ns:description', ns).text,
            'cost': price.find('ns:cost', ns).text,
            'unit': price.find('ns:unit', ns).text,
        }
        for price in root.findall('.//ns:price', ns)
    ]

    df = pd.DataFrame(prices).convert_dtypes(convert_string=True)
    df["cost"] = df["cost"].astype(float)
    return df


def filter_prices(prices: pd.DataFrame) -> pd.DataFrame:
    def return_the_only_number(s):
        # finds all ints / floats int string s
        numbers = re_findall(r"\d*\.*\d+", s)

        if len(numbers) == 1:
            return float(numbers[0])
        else:
            return 0

    prices = prices.copy()
    prices['amount'] = prices['description'].apply(return_the_only_number)

    prices_filtered = prices.query('unit not in ["box", "kit"] and amount > 0').reset_index(drop=True)
    prices_filtered.loc[prices_filtered['description'].str.contains(" mcg", na=False), "amount"] /= 1000

    # IDK what a unit of "unit" is and frankly im scared of them so they're gone
    prices_filtered = prices_filtered.query('description.str.contains("unit") == False').reset_index(drop=True)

    return prices_filtered


def int_to_db_string(int_value):
  """Converts an integer to a string with the format 'DB' followed by four zeros
  and then the integer value.

  Args:
    int_value: The integer to be converted.

  Returns:
    A string with the format 'DB00000' + str(int_value).
  """
  return "DB" + str(int_value).zfill(5)


import random
import copy


def generate_random(first_id, last_id, input_file, output_file):
    Et.register_namespace('', "http://www.drugbank.ca")
    tree = Et.parse(input_file)
    root = tree.getroot()
    existing_drugs = [drug for drug in root]

    namespace = "{http://www.drugbank.ca}"

    def get_all_tags():
        first_drug = existing_drugs[0] if existing_drugs else None
        return [child.tag.replace(namespace, "") for child in first_drug] if first_drug else []

    tags = get_all_tags()

    while "drugbank-id" in tags:
        tags.remove("drugbank-id")

    def get_random_subtree(tag, namespace):
        selected_drug = random.choice(existing_drugs)
        subtree = selected_drug.find(f'ns:{tag}', namespace)
        return subtree if subtree is not None else None

        # return copy.deepcopy(subtree) if subtree is not None else None

    namespace = {'ns': 'http://www.drugbank.ca'}

    for i in range(first_id, last_id + 1):  # + 12
        print("generated: ", i)
        new_drug = Et.Element("drug")
        Et.SubElement(new_drug, "drugbank-id", {"primary": "true", "created-by": "tolo"}).text = int_to_db_string(i)

        for tag in tags:
            subtree = get_random_subtree(tag, namespace)
            if subtree is not None:
                new_drug.append(subtree)
            else:
                new_drug.append(Et.Element(tag))  # Tworzy pusty tag XML
        root.append(new_drug)

    print("writing to: ", output_file)
    tree.write(output_file, encoding="UTF-8")