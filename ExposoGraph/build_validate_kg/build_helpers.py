import json
import re


def parse_static_js_const(file_path: str, var_name: str = "GRAPH_DATA") -> dict:
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
    # Matches: const VAR_NAME = { ... };
    pattern = rf"const\s+{var_name}\s*=\s*({{.*?}});"
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        raise ValueError(f"Variable {var_name} not found or improperly formatted.")
    js_object_string = match.group(1)
    return json.loads(js_object_string)
