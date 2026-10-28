
def parse_json(file_path: str, remove_nest: bool = False) -> dict:
    import json
    with open(file_path, 'r') as f:
        data = json.load(f)
    if remove_nest:
        if len(list(data.keys())) > 1:
            print("WARNING: more than one top-level key in JSON, taking the first one")
        data = data[list(data.keys())[0]]
    return data