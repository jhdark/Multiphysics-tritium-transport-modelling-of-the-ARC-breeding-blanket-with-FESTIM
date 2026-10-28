from festim_script import build_festim_model, load_openfoam_data, define_meshtags

if __name__ == "__main__":
    run_filename = "data/NewOF_symmBCs/arc.foam"
    run_time = 745

    # load data
    mesh, u, T, nut, t_prod, cfd_reader = load_openfoam_data(
        filename=run_filename, time=run_time
    )
    # convert tritprod from kg/m3/s to T/m3/s
    tritium_at_mass = 3.01604928  # Da or u
    tritium_at_mass = tritium_at_mass * 1.660539e-27  # kg
    t_prod.x.array[:] /= tritium_at_mass  # T/m3/s

    # define meshtags
    facet_meshtags, volume_meshtags = define_meshtags(cfd_reader)

    for delta_value in [10, 1, 0.1, 0.01]:
        print(f"Running with delta = {delta_value}")
        my_model = build_festim_model(
            mesh=mesh,
            facet_meshtags=facet_meshtags,
            volume_meshtags=volume_meshtags,
            u=u,
            T=T,
            nut=nut,
            tritium_production_rate=t_prod,
            delta=delta_value,
            export_foldername=f"results/symmBCs/delta_testing/delta_{delta_value}",
            steady=True,
        )

        print("Initialising and running festim model")
        my_model.initialise()
        my_model.run()
