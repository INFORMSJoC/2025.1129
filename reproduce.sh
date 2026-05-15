# Reproduce Table 1 of the paper.
uv run src_new/LPT.py

# Reproduce Table 2 of the paper.
uv run src_new/MLPT.py

# Reproduce Table 3 of the paper.
uv run src_new/FFD.py
uv run src_new/FFD_reversed.py
uv run src_new/FFD_bidirectional.py
uv run src_new/COMBINE.py
uv run src_new/COMBINE_bidirectional.py