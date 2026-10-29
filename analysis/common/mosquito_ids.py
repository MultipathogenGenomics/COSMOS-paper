
import re

def extract_genus_from_target_id(target_id, report_error: bool = False):
    id_eles = re.split(r"[-_\s]+", target_id)
    if "Aedes" in id_eles:
        return "Aedes"
    elif "Anopheles" in id_eles:
        return "Anopheles"
    elif "Culex" in id_eles:
        return "Culex"
    elif "Toxorhynchites" in id_eles:
        return "Toxorhynchites"
    elif "Ochlerotatus" in id_eles:
        return "Ochlerotatus"
    elif "Haemagogus" in id_eles:
        return "Haemagogus"
    elif "Mansonia" in id_eles:
        return "Mansonia"
    elif "Psorophora" in id_eles:
        return "Psorophora"
    elif "Coquillettidia" in id_eles:
        return "Coquillettidia"
    elif "Uranotaenia" in id_eles:
        return "Uranotaenia"
    elif "Tripteroides" in id_eles:
        return "Tripteroides"
    elif "Verrallina" in id_eles:
        return "Verrallina"
    elif "Sabethes" in id_eles:
        return "Sabethes"
    else:
        if report_error:
            print(f"Could not extract genus from target ID: {target_id}")
        return "unknown"
    
def return_species_from_target_id(target_id):
    id_eles = re.split(r"[-_\s]+", target_id)
    genus = extract_genus_from_target_id(target_id)
    genus_idx = id_eles.index(genus) if genus in id_eles else -1
    if genus_idx != -1 and genus_idx + 1 < len(id_eles):
        species = id_eles[genus_idx + 1]
        return genus + " " + species
    else:
        return "unknown"
    
def shorted_genus(species):
    replace_dict = {
        "Aedes": "Ae.",
        "Culex": "Cx.",
        "Anopheles": "An.",
        "Mansonia": "Ma.",
        "Coquillettidia": "Cq.",
        "Ochlerotatus": "Ae.",
        "Haemagogus": "Hg.",
        "Tripteroides": "Tr.",
        "Psorophora": "Ps.",
        "Uranotaenia": "Ur.",
        "Sabethes": "Sa.",
        "Verrallina": "Ve.",
        "Toxorhynchites": "Tx."
    }
    current_genus = species.split(" ")[0]
    shortened_genus = species.replace(current_genus, replace_dict[current_genus])
    return shortened_genus


# Retained from the original compare_kallisto_kma_minimap naming helpers.
def shortened_mosquito_target_id(target_id):
    # split by dashes or underscores
    id_eles = re.split("-|_", target_id)
    id_eles_to_keep = [ele for ele in id_eles if not ele.isnumeric() and ele not in
                       ["mosquito", "NC", "cytochrome", "c", "oxidase", "subunit",
                        "mitochondrial", "mitochondrion", "complete", "genome", "COI",
                        "gene", "partial", "cds", "COX1"]]
    return " ".join(id_eles_to_keep)
