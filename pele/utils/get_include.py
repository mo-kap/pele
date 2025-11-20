import pathlib

def get_include() -> pathlib.Path:
    import pele
    include_path = pathlib.Path(pele.__file__).parent / "include"
    current_path = pathlib.Path(".").absolute()
    if include_path.is_relative_to(current_path):
        return include_path.relative_to(current_path)
    else:
        return include_path


def get_lammps_include() -> pathlib.Path:
    import lammps
    lammps_package_path = pathlib.Path(lammps.__file__).parent
    pip_include_path = lammps_package_path / "include" / "lammps"
    if pip_include_path.exists():
        return pip_include_path
    local_include_path = lammps_package_path.parent.parent.parent.parent / "include" / "lammps"
    if local_include_path.exists():
        return local_include_path
    return None
