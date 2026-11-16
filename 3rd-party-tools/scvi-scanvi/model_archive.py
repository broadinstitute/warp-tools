"""The scANVI model archive: one .tar.gz holding a reference's scVI model and its SCANVI models.

  scanvi_model/
    reference.txt             the reference id (below): the model's name in every output
    levels.txt                the label levels, one per line, in order; the first is SCP's default
    scvi/model.pt             the scVI model the SCANVI models were trained from (optional)
    scanvi/<level>/model.pt   one SCANVI model per label level (an obs column of the reference,
                              e.g. class, subclass), all trained from that scVI model

Models are saved without data (save_anndata=False). write_archive() packs saved model directories;
read_archive() unpacks one. A level name becomes an obs column and part of obsm keys and paths, so
it is limited to letters, digits, '_', '.' and '-'.

reference_id(path, uns) is the reference's uns['title'] (the Allen taxonomy and CELLxGENE schemas
both define it), or else its file name, made file-name safe, plus the first 8 hex digits of the
file's sha256: e.g. 'BICCN-mouse_whole_brain_consensus-1a2b3c4d'. The title says which reference;
the hash tells apart subsets and rebuilds that share a title.
"""
import hashlib
import io
import os
import re
import tarfile

ROOT = "scanvi_model"


def reference_id(path, uns):
    name = str(uns.get("title") or "") or os.path.basename(path).removesuffix(".h5ad")
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return f"{name}-{digest.hexdigest()[:8]}"


def check_levels(levels):
    bad = [level for level in levels if not re.fullmatch(r"[A-Za-z0-9_.-]+", level)]
    if not levels or bad or len(set(levels)) != len(levels):
        raise SystemExit(f"ERROR: label levels must be one or more distinct names of letters, digits, "
                         f"'_', '.' and '-'; got {list(levels)}")


def write_archive(out_tar, reference, scvi_dir, scanvi_dirs):
    """Pack saved models into out_tar. scanvi_dirs maps each level, in order, to its saved SCANVI
    model directory; scvi_dir may be None."""
    check_levels(list(scanvi_dirs))
    with tarfile.open(out_tar, "w:gz") as tf:
        if scvi_dir:
            tf.add(os.path.join(scvi_dir, "model.pt"), f"{ROOT}/scvi/model.pt")
        for level, model_dir in scanvi_dirs.items():
            tf.add(os.path.join(model_dir, "model.pt"), f"{ROOT}/scanvi/{level}/model.pt")
        for name, text in (("reference.txt", reference), ("levels.txt", "\n".join(scanvi_dirs))):
            data = (text + "\n").encode()
            info = tarfile.TarInfo(f"{ROOT}/{name}")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def read_archive(model_tar, dest="scanvi_model_in"):
    """Unpack model_tar into dest. Returns (reference id, {level: SCANVI model directory}) in level
    order."""
    with tarfile.open(model_tar) as tf:
        tf.extractall(dest, filter="data")
    root = os.path.join(dest, ROOT)
    if not os.path.isfile(os.path.join(root, "levels.txt")):
        raise SystemExit(f"ERROR: {os.path.basename(model_tar)} is not a scANVI model archive: it has no "
                         f"{ROOT}/levels.txt (see model_archive.py for the layout).")
    with open(os.path.join(root, "levels.txt")) as f:
        levels = f.read().split()
    with open(os.path.join(root, "reference.txt")) as f:
        reference = f.read().strip()
    return reference, {level: os.path.join(root, "scanvi", level) for level in levels}
