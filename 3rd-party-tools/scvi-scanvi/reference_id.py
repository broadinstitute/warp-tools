"""A unique, readable name for the reference a model was trained on, carried inside the model.

reference_id(path, uns) is the reference's uns['title'] (the Allen taxonomy and CELLxGENE schemas
both define it), or else its file name, made file-name safe, plus the first 8 hex digits of the
file's sha256: e.g. 'BICCN-mouse_whole_brain_consensus-1a2b3c4d'. The title says which reference;
the hash tells apart subsets and rebuilds that share a title.

A saved model carries it in reference.txt next to model.pt (write_reference / read_reference).
"""
import hashlib
import os
import re


def reference_id(path, uns):
    name = str(uns.get("title") or "") or os.path.basename(path).removesuffix(".h5ad")
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return f"{name}-{digest.hexdigest()[:8]}"


def write_reference(model_dir, reference):
    with open(os.path.join(model_dir, "reference.txt"), "w") as f:
        f.write(reference + "\n")


def read_reference(model_dir):
    path = os.path.join(model_dir, "reference.txt")
    return open(path).read().strip() if os.path.exists(path) else None
