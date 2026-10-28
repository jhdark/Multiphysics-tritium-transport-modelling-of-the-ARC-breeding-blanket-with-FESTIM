from pathlib import Path

from mpi4py import MPI

import dolfinx
import io4dolfinx
from foam2dolfinx import OpenFOAMReader
from dolfinx.io import VTXWriter
from data_utils import ensure_openfoam_case


def read_data_from_foam_file(foam_file, time_value):
    ensure_openfoam_case()

    reader = OpenFOAMReader(filename=foam_file, cell_type=12)

    print("Reading mesh and functions from OpenFOAM...")

    u = reader.create_dolfinx_function_with_point_data(t=time_value, name="U")
    T = reader.create_dolfinx_function_with_point_data(t=time_value, name="T")
    nut = reader.create_dolfinx_function_with_point_data(t=time_value, name="nut")

    tritProd = reader.create_dolfinx_function_with_point_data(
        t=time_value, name="tritProd"
    )
    nuc_heating = reader.create_dolfinx_function_with_point_data(
        t=time_value, name="heatGen"
    )

    mesh = reader.dolfinx_meshes_dict["default"]

    writer1 = VTXWriter(
        MPI.COMM_WORLD, "data/NewOF_symmBCs/CG_fields/velocity.bp", u, "BP5"
    )
    writer1.write(t=0)
    writer2 = VTXWriter(
        MPI.COMM_WORLD, "data/NewOF_symmBCs/CG_fields/temperature.bp", T, "BP5"
    )
    writer2.write(t=0)
    writer3 = VTXWriter(
        MPI.COMM_WORLD, "data/NewOF_symmBCs/CG_fields/nut.bp", nut, "BP5"
    )
    writer3.write(t=0)
    writer4 = VTXWriter(
        MPI.COMM_WORLD, "data/NewOF_symmBCs/CG_fields/t_prod.bp", tritProd, "BP5"
    )
    writer4.write(t=0)
    writer5 = VTXWriter(
        MPI.COMM_WORLD, "data/NewOF_symmBCs/CG_fields/heat_gen.bp", nuc_heating, "BP5"
    )
    writer5.write(t=0)

    functions = [u, T, nut, tritProd, nuc_heating]

    facet_meshtags = reader.create_facet_meshtags()

    return mesh, functions, facet_meshtags


def write_checkpoints(
    checkpoint_file: Path,
    mesh: dolfinx.mesh.Mesh,
    funcs: list[dolfinx.fem.Function],
    facet_meshtags: dolfinx.mesh.MeshTags,
):
    io4dolfinx.write_mesh(checkpoint_file, mesh)

    names = ["U", "T", "tritProd", "nut"]
    for func, name in zip(funcs, names):
        io4dolfinx.write_function(checkpoint_file, func, time=0.0, name=name)

    io4dolfinx.write_meshtags(
        checkpoint_file, mesh, facet_meshtags, meshtag_name="facet_tags"
    )


if __name__ == "__main__":
    comm = MPI.COMM_WORLD

    # foam_file = Path("data/run1_April28_1.4Mcells/run1_OpenFOAM_tritSolve/arc.foam")
    foam_file = Path("data/NewOF_symmBCs/arc.foam")
    time_value = 745
    checkpoint_file = Path("data/NewOF_symmBCs/arc_checkpoint_cell.bp")

    print("Reading data from OpenFOAM...")

    dolfinx_mesh, functions, facet_meshtags = read_data_from_foam_file(
        foam_file, time_value
    )

    # print("Writing checkpoint data")

    # write_checkpoints(
    #     checkpoint_file=checkpoint_file,
    #     mesh=dolfinx_mesh,
    #     funcs=functions,
    #     facet_meshtags=facet_meshtags,
    # )
