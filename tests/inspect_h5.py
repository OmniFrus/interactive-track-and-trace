import h5py
import os

for filename in ["hydrodynamic_U.h5", "hydrodynamic_V.h5", "grid.h5", "shore_distance.h5", "shore_distance_on_grid.h5"]:
    filepath = os.path.join("../../data", filename)
    print(f"\nInspecting {filepath}")
    try:
        with h5py.File(filepath, "r") as f:
            def walk(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  {name}: shape={obj.shape}")
            f.visititems(walk)
    except Exception as e:
        print(f"  Failed to open or read {filepath}: {e}")