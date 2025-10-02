import requests
import pandas as pd
def searchProteinByName(protein_names, limit=1):
    """
    Search UniProt by a list of protein names (restricted to mouse) 
    and return a dict of {name: [UniProt IDs]}.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    results_dict = {}

    for name in protein_names:
        query = f'{name} AND (taxonomy_id:10090)'  # restrict to mouse
        params = {
            "query": query,
            "format": "json",
            "size": limit
        }
        response = requests.get(url, params=params)
        response.raise_for_status()
        results = response.json()
        
        ids = [entry["primaryAccession"] for entry in results.get("results", [])]
        results_dict[name] = ids
    
    return results_dict

proteinIds = searchProteinByName(["Pparg","Hpn","Phb2","Tgfbr1","Aldh1a2","Lama5","Tnfaip3","Arg1","Stxbp4","Nr1d1","Zfp36l1","Ptch1","Ctsl","Cnmd",
                                  "Sema5a","Rida","Cldn1","Hes1","Prkdc","Cxadr","Robo1","Vegfa","Xdh","Pdpk1","Lims2","Pold4","Cflar","Stat1",
                                  "Eng","Egfl7","Jag1","Ecm1","Nr4a3","Errfi1","Mtor","Ptn","Wnt7a","Itpr1","Fgfr1","Phip","Smad3","Mst1",
                                  "Sgpp2","Atp7a","Bmp6","Tgfb2","Thbs1","Cdk6","Esrp1","Prkd2","Dicer1","Mmrn2","Aqp11","Zfp36","Ar","Agtr1a",
                                  "Rictor","Gja1","Tacstd2","Gas1","Ceacam2","Fut2","B2m","Cdh3","Erbb4","Erbb2","Serpinb5","Ceacam1","Col4a3","Esrp2","Rian"])
print(proteinIds)

def getProteinDescriptions(uniprotId):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprotId}.json"
    response = requests.get(url)
    response.raise_for_status()
    output = response.json()
    return output['references']#,output['comments'][0]['texts'] #up to 2
#res = getProteinDescriptions(proteinIds["Lama5"][0])#['texts'][0]['value']
#for val in res:
#    print(val['citation']["id"])



rows = []
for name, ids in proteinIds.items():
    for uid in ids:
      try:
        info = getProteinDescriptions(uid)
        listOfCit = [val['citation']['id'] for val in info]
        rows.append({"geneName":name,
                     "references": listOfCit})
      except (KeyError,IndexError):
        rows.append({geneName:name,
                     "references":"N/A"})

proteinTable = pd.DataFrame(rows)
proteinTable.to_csv("RegEpCellProlifDel34Liver.csv")
