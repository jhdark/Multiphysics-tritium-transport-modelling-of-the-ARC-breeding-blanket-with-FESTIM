from dolfinx.mesh import meshtags, exterior_facet_indices
from scipy.spatial import cKDTree
import numpy as np
import ufl
from dolfinx import fem
import festim as F
from dolfinx import cpp as _cpp
from scifem import assemble_scalar


def tag_boundary_patch(dolfinx_mesh, patch_dataset, patch_id, tol=1e-6):
    fdim = dolfinx_mesh.topology.dim - 1
    dolfinx_mesh.topology.create_connectivity(fdim, 0)
    dolfinx_mesh.topology.create_connectivity(0, fdim)
    dolfinx_mesh.topology.create_connectivity(fdim, dolfinx_mesh.topology.dim)

    facet_indices = exterior_facet_indices(dolfinx_mesh.topology)
    x = dolfinx_mesh.geometry.x
    patch_points = patch_dataset.points
    tree = cKDTree(x)
    matched_vertex_indices = tree.query_ball_point(patch_points, tol)
    matched_vertex_indices = list(set(i for sub in matched_vertex_indices for i in sub))

    matched_facets = []
    for facet in facet_indices:
        vertices = dolfinx_mesh.topology.connectivity(fdim, 0).links(facet)
        if all(v in matched_vertex_indices for v in vertices):
            matched_facets.append(facet)

    # print(f"Tagging {len(matched_facets)} facets for patch ID {patch_id}")
    return np.array(matched_facets, dtype=np.int32), np.full(
        len(matched_facets), patch_id, dtype=np.int32
    )


def define_meshtags(cfd_reader):
    OF_multiblock = cfd_reader.reader.read()
    boundary = OF_multiblock["boundary"]
    mesh = cfd_reader.dolfinx_meshes_dict["default"]

    all_facets = np.array([], dtype=np.int32)
    all_tags = np.array([], dtype=np.int32)

    for i, name in enumerate(boundary.keys()):
        if name not in ["inlet", "outlet"]:
            continue
        facets, tags = tag_boundary_patch(mesh, boundary[name], i + 1)
        all_facets = np.concatenate([all_facets, facets])
        all_tags = np.concatenate([all_tags, tags])

        # print(f"Tagging {len(facets)} facets for patch {name} with ID {i + 1}")

    facet_tags = meshtags(
        mesh,
        mesh.topology.dim - 1,
        all_facets,
        all_tags,
    )

    num_cells = mesh.topology.index_map(mesh.topology.dim).size_local
    mesh_cell_indices = np.arange(num_cells, dtype=np.int32)
    tags_volumes = np.full(num_cells, 1, dtype=np.int32)
    volume_meshtags = meshtags(mesh, mesh.topology.dim, mesh_cell_indices, tags_volumes)

    num_cells = (
        mesh.topology.index_map(mesh.topology.dim).size_local
        + mesh.topology.index_map(mesh.topology.dim).num_ghosts
    )

    print(f"number of cells in mesh: {num_cells}")
    return (
        facet_tags,
        volume_meshtags,
    )


def evaluate_stabalisation_term(mesh, u, delta):
    """See more at https://www.comsol.com/blogs/understanding-stabilization-methods"""

    # evaluate Cell size
    tdim = mesh.topology.dim
    num_cells = mesh.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)
    mesh_ = _cpp.mesh.Mesh_float64(
        mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
    )
    h = _cpp.mesh.h(mesh_, tdim, cells)
    V0 = fem.functionspace(mesh, ("DG", 0))
    h_as_function = fem.Function(V0)
    h_as_function.x.array[:] = h

    # Compute magnitude of velocity
    v_mag = ufl.sqrt(ufl.dot(u, u))

    D_art = delta * v_mag * h_as_function

    return D_art


class Advection_flux(F.SurfaceFlux):
    def __init__(self, field, surface, filename, velocity_field):
        super().__init__(field=field, surface=surface, filename=filename)

        self.velocity_field = velocity_field

    @property
    def title(self):
        return f"{self.field.name} advective flux surface {self.surface.id}"

    def compute(self, u, ds: ufl.Measure, entity_maps=None):
        """Computes the value of the flux at the surface

        Args:
            ds (ufl.Measure): surface measure of the model
        """

        # obtain mesh normal from field
        # if case multispecies, solution is an index, use sub_function_space
        if isinstance(u, ufl.indexed.Indexed):
            mesh = self.field.sub_function_space.mesh
        else:
            mesh = u.function_space.mesh
        n = ufl.FacetNormal(mesh)

        surface_flux = assemble_scalar(
            fem.form(
                -self.D * ufl.dot(ufl.grad(u), n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )
        advective_flux = assemble_scalar(
            fem.form(
                u * ufl.inner(self.velocity_field, n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        self.value = surface_flux + advective_flux
        self.data.append(self.value)
