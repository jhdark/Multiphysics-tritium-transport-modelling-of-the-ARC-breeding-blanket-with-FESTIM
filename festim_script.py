import festim as F

from dolfinx.log import set_log_level, LogLevel
from foam2dolfinx import OpenFOAMReader
from dolfinx import fem
import ufl
from helpers import define_meshtags, evaluate_stabalisation_term, Advection_flux
from data_utils import ensure_openfoam_case
from dolfinx.io import VTXWriter
from mpi4py import MPI


def load_openfoam_data(filename, time):
    ensure_openfoam_case()
    print("Reading openfoam data")
    cfd_reader = OpenFOAMReader(f"{filename}", cell_type=12)
    u = cfd_reader.create_dolfinx_function_with_point_data(t=time, name="U")
    T = cfd_reader.create_dolfinx_function_with_point_data(t=time, name="T")
    nut = cfd_reader.create_dolfinx_function_with_point_data(t=time, name="nut")
    tritium_production_rate = cfd_reader.create_dolfinx_function_with_point_data(
        t=time, name="tritProd"
    )

    mesh = cfd_reader.dolfinx_meshes_dict["default"]

    # facet_tags = cfd_reader.create_facet_meshtags(t=time)

    return mesh, u, T, nut, tritium_production_rate, cfd_reader


def build_festim_model(
    mesh,
    facet_meshtags,
    volume_meshtags,
    u,
    T,
    nut,
    tritium_production_rate,
    delta=0.1,
    Sc_no=0.5,
    export_foldername="results",
    steady=True,
):
    print("Building festim model")

    my_model = F.HydrogenTransportProblem()
    my_mesh = F.Mesh(mesh=mesh)
    my_model.facet_meshtags = facet_meshtags
    my_model.volume_meshtags = volume_meshtags
    my_model.mesh = my_mesh

    # define species
    H = F.Species("H", mobile=True)
    my_model.species = [H]

    # define species diffusivities and materials
    D_0 = 6.75e-06
    E_D = 0.615
    D_fickian = D_0 * ufl.exp(-E_D / (F.k_B * T))

    V = fem.functionspace(mesh, ("CG", 1))
    D_flibe_fick = fem.Function(V)
    D_flibe_fick.interpolate(fem.Expression(D_fickian, V.element.interpolation_points))

    # writer1 = VTXWriter(
    #     MPI.COMM_WORLD, "data/NewOF_symmBCs/D_flibe_fick.bp", D_flibe_fick, "BP5"
    # )
    # writer1.write(t=0)

    # eval turbulent diffusion
    turb_Schmidt_no = Sc_no
    D_turb = nut / turb_Schmidt_no

    # eval stablisation term
    D_art = evaluate_stabalisation_term(mesh=mesh, u=u, delta=delta)

    D_expr = D_fickian + D_turb + D_art
    V = fem.functionspace(mesh, ("CG", 1))
    D_flibe = fem.Function(V)
    D_flibe.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    # writer2 = VTXWriter(MPI.COMM_WORLD, "data/NewOF_symmBCs/D_flibe.bp", D_flibe, "BP5")
    # writer2.write(t=0)

    flibe = F.Material(D=D_flibe)

    # define boundary conditions
    inlet = F.SurfaceSubdomain(id=1)
    outlet = F.SurfaceSubdomain(id=2)
    fluid = F.VolumeSubdomain(id=1, material=flibe)
    my_model.subdomains = [inlet, outlet, fluid]

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, species=H, value=0),
    ]

    # define temperature
    my_model.temperature = T

    # define sources

    my_model.sources = [
        F.ParticleSource(volume=fluid, species=H, value=tritium_production_rate),
    ]

    # define advection term
    my_model.advection_terms = [F.AdvectionTerm(velocity=u, species=H, subdomain=fluid)]

    if not steady:
        # define model settings
        dt = F.Stepsize(
            initial_value=10,
            growth_factor=1.05,
            cutback_factor=0.9,
            target_nb_iterations=10,
        )
        final_time = 3600 * 1.5  # 1.5 hours
        my_model.settings = F.Settings(
            atol=1e10, rtol=1e-15, transient=True, stepsize=dt, final_time=final_time
        )

    else:
        my_model.settings = F.Settings(
            atol=1e5, rtol=1e-12, max_iterations=5, transient=False
        )

    # define exports
    outlet_surface_flux = Advection_flux(
        field=H,
        surface=outlet,
        filename=f"{export_foldername}/outlet_surface_flux.csv",
        velocity_field=u,
    )
    inlet_surface_flux = F.SurfaceFlux(
        field=H,
        surface=inlet,
        filename=f"{export_foldername}/inlet_surface_flux.csv",
    )
    total_inventory = F.TotalVolume(
        field=H,
        volume=fluid,
        filename=f"{export_foldername}/total_inventory.csv",
    )
    outlet_concentration = F.AverageSurface(
        field=H,
        surface=outlet,
        filename=f"{export_foldername}/outlet_average_concentration.csv",
    )
    mobile_field_vtx = F.VTXSpeciesExport(
        filename=f"{export_foldername}/mobile_concentration.bp",
        field=H,
    )
    mobile_field_xdmf = F.XDMFExport(
        filename=f"{export_foldername}/mobile_concentration.xdmf",
        field=H,
    )

    my_model.exports = [
        outlet_surface_flux,
        inlet_surface_flux,
        total_inventory,
        outlet_concentration,
        # mobile_field_vtx,
        mobile_field_xdmf,
    ]

    return my_model


if __name__ == "__main__":
    # set_log_level(LogLevel.INFO)

    run1_filename = "data/run1_April28_1.4Mcells/run1_OpenFOAM_tritSolve/arc.foam"
    run1_time = 500

    run2_filename = "data/run2_May28_6Mcells/run2_OFcase/arc.foam"
    run2_time = 44.65

    run3_filename = "data/run3_May28_11.3Mcells/run3_OFcase/arc.foam"
    run3_time = 407

    new_file = "data/NewOF_symmBCs/arc.foam"
    new_time = 745

    # load data
    mesh, u, T, nut, t_prod, cfd_reader = load_openfoam_data(
        filename=new_file, time=new_time
    )

    # convert tritprod from kg/m3/s to T/m3/s
    tritium_at_mass = 3.01604928  # Da or u
    tritium_at_mass = tritium_at_mass * 1.660539e-27  # kg
    t_prod.x.array[:] /= tritium_at_mass  # T/m3/s

    # define meshtags
    facet_meshtags, volume_meshtags = define_meshtags(cfd_reader)

    # build festim model
    my_model = build_festim_model(
        mesh=mesh,
        facet_meshtags=facet_meshtags,
        volume_meshtags=volume_meshtags,
        u=u,
        T=T,
        nut=nut,
        tritium_production_rate=t_prod,
        steady=True,
        export_foldername="results/symmBCs/steady",
        delta=0.1,
    )

    print("Initialising and running festim model")
    my_model.initialise()
    set_log_level(LogLevel.INFO)
    import time

    start = time.time()
    my_model.run()
    end = time.time() - start
    print(f"Simulation took {end} seconds")

    import dolfinx
    import numpy as np
    from mpi4py import MPI
    from petsc4py import PETSc

    def compute_consistent_flux(form: ufl.Form, bc: dolfinx.fem.DirichletBC):
        residual_flux = dolfinx.fem.form(form)
        residual = dolfinx.fem.assemble_vector(residual_flux)
        residual.scatter_reverse(dolfinx.la.InsertMode.add)
        residual.scatter_forward()
        bc_fs = bc.function_space
        res_fs = residual_flux.function_spaces
        assert len(res_fs) == 1
        if not res_fs[0].contains(bc_fs):
            return 0
        dofs, num_owned = bc.dof_indices()
        owned_dofs = dofs[:num_owned]
        local_flux = np.sum(residual.array[owned_dofs])
        return residual_flux.mesh.comm.allreduce(local_flux, op=MPI.SUM)

    # flux_in_new = compute_consistent_flux(F_mobile, bc_inlet)
    flux_in_new = compute_consistent_flux(my_model.formulation, my_model.bc_forms[0])
    print(f"Consistent inlet flux: {flux_in_new:.3e} T/s")

    outlet_facets = facet_meshtags.find(2)
    outlet_dofs = fem.locate_dofs_topological(
        my_model.function_space, my_model.facet_meshtags.dim, outlet_facets
    )
    bc_temp = fem.dirichletbc(
        dolfinx.fem.Constant(my_model.mesh.mesh, PETSc.ScalarType(1.0)),
        outlet_dofs,
        my_model.species[0].sub_function_space,
    )

    flux_out_new = compute_consistent_flux(my_model.formulation, bc_temp)

    import scifem

    total_t_prod = scifem.assemble_scalar(
        my_model.sources[0].value.fenics_object * my_model.dx
    )
    # print(f"Total tritium produced: {total_t_prod:.3e} T/s")
    out_flux = my_model.exports[0].data[-1]
    inlet_flux = my_model.exports[1].data[-1]
    # print(f"Total tritium out (from outlet flux): {out_flux:.3e} T/s")
    # print(f"Total tritium out (from inlet flux): {inlet_flux:.3e} T/s")

    # convert from kg/s to mg/s
    total_t_prod *= 1e6
    out_flux *= 1e6
    inlet_flux *= 1e6

    # convert from 15 deg to 360 deg
    total_t_prod *= 72
    out_flux *= 72
    inlet_flux *= 72

    print(f"Total tritium produced: {total_t_prod:.3e} mg/s")
    print(f"Total tritium out (from outlet flux): {out_flux:.3e} mg/s")
    print(f"Total tritium out (from inlet flux): {inlet_flux:.3e} mg/s")

    # print(f"Consistent inlet flux: {flux_in_new:.3e} T/s")
    # print(f"Consistent outlet flux: {flux_out_new:.3e} T/s")
    total_out = inlet_flux + out_flux
    # print(f"Total tritium out: {total_out:.3e} T/s")

    print(f"Relative difference: {(total_t_prod - total_out) / total_t_prod:.2%}")
