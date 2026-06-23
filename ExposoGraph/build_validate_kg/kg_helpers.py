import json


def load_json(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        result = json.load(file)
    return result
