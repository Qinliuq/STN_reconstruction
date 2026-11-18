"""The model side (model/) must never import the experimental-data loader."""
import glob
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_model_code_does_not_load_data():
    for f in glob.glob(os.path.join(HERE, 'model', '**', '*.py'), recursive=True):
        src = open(f, encoding='utf-8').read()
        assert not re.search(r'data_io|heldout_bundle|HELDOUT_DATA_DIR', src), f
