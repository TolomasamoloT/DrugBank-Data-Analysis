from fastapi import FastAPI
from pydantic import BaseModel
import pandas as pd
from my_lib import extract_pathway_ids, extract_drugs
import xml.etree.ElementTree as Et


# Load the XML file
tree = Et.parse('drugbank_partial.xml')
root = tree.getroot()

# Use the namespace to access elements
namespace = {'ns': 'http://www.drugbank.ca'}

drugs_df = extract_drugs(root, namespace)

# Initialize FastAPI app
app = FastAPI()

# Create a DataFrame from the given data
df = extract_pathway_ids(root, namespace, drugs_df)

# Define request model
class DrugRequest(BaseModel):
    drug_id: str

# Define endpoint
@app.post("/get_drug_count/")
async def get_drug_count(request: DrugRequest):
    drug_id = request.drug_id
    if drug_id in df.index:
        return {"count": int(df.loc[drug_id, "count"])}
    return {"error": "Drug not found"}